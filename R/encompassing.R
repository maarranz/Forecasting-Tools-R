# Forecast encompassing -----------------------------------------------------------
#
# Matches encompassing_test() in the Python package Forecasting-Tools.

#' Forecast encompassing tests
#'
#' Does one forecast contain all the useful (linear) information in the other?
#' With \eqn{e_1 = y - \hat y_1}, \eqn{e_2 = y - \hat y_2} and
#' \eqn{d = \hat y_1 - \hat y_2}:
#'
#' * `method = "difference"` (default): \eqn{e_1 = \beta d + u} (direction
#'   `"1_encompasses_2"`) or \eqn{e_2 = \beta d + u} (`"2_encompasses_1"`),
#'   testing \eqn{\beta = 0}. The sign of d is the same in both directions.
#' * `method = "cross_forecast"`: \eqn{e_1 = \delta \hat y_2 + u} or
#'   \eqn{e_2 = \gamma \hat y_1 + u}, testing that coefficient is zero. These are
#'   different orthogonality restrictions, not generally equivalent to the
#'   difference regressions.
#' * `method = "unrestricted"`: \eqn{y = a_0 + a_1 \hat y_1 + a_2 \hat y_2 + u},
#'   testing \eqn{(a_1, a_2) = (1, 0)} (direction 1) or \eqn{(0, 1)}
#'   (direction 2) jointly; the intercept is always unrestricted.
#'
#' The difference and cross-forecast regressions have no intercept unless
#' `include_intercept = TRUE`, which changes the restriction from an
#' uncentered to a centered orthogonality condition.
#'
#' **Inference** is always HAC: Newey-West/Bartlett covariance without
#' small-sample adjustment, two-sided normal tests for one restriction and an
#' asymptotic chi-square(2) for the unrestricted joint test.
#'
#' The two directions are separate hypotheses and must be tested separately.
#' Failure to reject means that the encompassing restriction cannot be
#' rejected, not that it holds; rejection does not mean the competing forecast
#' is more accurate. Identical forecasts make the difference and unrestricted
#' designs unidentified.
#'
#' @inheritParams mgn_test
#' @param method `"difference"`, `"cross_forecast"` or `"unrestricted"`.
#' @param direction `"1_encompasses_2"` or `"2_encompasses_1"`.
#' @param include_intercept Add an unrestricted intercept to the difference or
#'   cross-forecast regression. Must be `FALSE` for `"unrestricted"`, which
#'   always has one.
#' @param maxlags Bartlett bandwidth (non-negative integer below n), or `NULL`
#'   for \eqn{\min(n-1, \lfloor 4 (n/100)^{2/9} \rfloor)}.
#' @return A list of class `forecastingtools_test` with elements `method`,
#'   `direction`, `nobs`, `include_intercept` (whether the fitted regression has
#'   one), `null_hypothesis`, `statistic`, `distribution` (`"normal"` or
#'   `"chi2"`), `df` (`NA` or 2), `df_resid`, `n_restrictions`, `p.value`,
#'   `alpha`, `reject_null`, `cov_type` (`"HAC"`), `maxlags` and `coefficients`
#'   (terms `intercept` when fitted and `forecast_difference`,
#'   `competing_forecast`, or `forecast1` and `forecast2`).
#' @examples
#' actual <- c(10, 20, 30, 40, 50)
#' f1 <- c(9, 22, 27, 39, 48)
#' f2 <- c(12, 19, 31, 37, 51)
#' encompassing_test(actual, f1, f2, maxlags = 1)
#' encompassing_test(actual, f1, f2, method = "unrestricted")$df
#' @export
encompassing_test <- function(actual, forecast1, forecast2, method = "difference",
                              direction = "1_encompasses_2", include_intercept = FALSE,
                              maxlags = NULL, alpha = 0.05) {
  method <- .check_choice(method, c("difference", "cross_forecast", "unrestricted"),
                          "method")
  direction <- .check_choice(direction, c("1_encompasses_2", "2_encompasses_1"),
                             "direction")
  include_intercept <- .check_flag(include_intercept, "include_intercept")
  if (method == "unrestricted" && include_intercept) {
    stop("The unrestricted regression always includes an intercept; leave ",
         "`include_intercept = FALSE`.", call. = FALSE)
  }
  alpha <- .check_unit(alpha, "alpha")
  maxlags <- .check_maxlags(maxlags)
  s <- .three(actual, forecast1, forecast2)
  y <- s$actual; f1 <- s$forecast1; f2 <- s$forecast2
  n <- length(y)
  if (!is.null(maxlags)) .bandwidth(maxlags, n)
  joint <- method == "unrestricted"
  has_intercept <- joint || include_intercept
  first <- direction == "1_encompasses_2"
  if (joint) {
    target <- y
    X <- cbind(1, f1, f2)
    terms <- c("intercept", "forecast1", "forecast2")
    R <- rbind(c(0, 1, 0), c(0, 0, 1))
    r <- if (first) c(1, 0) else c(0, 1)
    null_label <- if (first) "forecast1 = 1 and forecast2 = 0"
                  else "forecast1 = 0 and forecast2 = 1"
  } else {
    target <- y - (if (first) f1 else f2)
    regressor <- if (method == "difference") f1 - f2 else if (first) f2 else f1
    label <- if (method == "difference") "forecast_difference" else "competing_forecast"
    terms <- c(if (has_intercept) "intercept", label)
    X <- if (has_intercept) cbind(1, regressor) else matrix(regressor, ncol = 1)
    R <- matrix(c(rep(0, length(terms) - 1), 1), nrow = 1)
    r <- 0
    null_label <- paste(label, "= 0")
  }
  k <- ncol(X)
  if (n <= k) {
    stop("Encompassing inference requires more observations than regression ",
         "parameters (", k, ").", call. = FALSE)
  }
  if (!all(is.finite(X)) || !all(is.finite(target))) {
    stop("The encompassing design or errors are non-finite; rescale the ",
         "observations.", call. = FALSE)
  }
  if (.rank(X) != k) {
    stop("The encompassing design is rank-deficient: the regressors are not ",
         "identified (for example, identical forecasts).", call. = FALSE)
  }
  fit <- .ols(target, X)
  res_scale <- max(max(abs(target)), max(abs(stats::fitted(fit))))
  if (res_scale == 0 ||
      max(abs(stats::residuals(fit))) / res_scale <= n * .Machine$double.eps) {
    stop("The encompassing residual variance is zero or numerically degenerate.",
         call. = FALSE)
  }
  lags <- .bandwidth(maxlags, n)
  V <- .check_vcov(.vcov(fit, "HAC", lags), k, "Encompassing HAC parameter")
  b <- unname(stats::coef(fit)); se <- sqrt(diag(V))
  if (joint) {
    statistic <- .wald(b, V, R, r)
    p <- stats::pchisq(statistic, 2, lower.tail = FALSE)
  } else {
    statistic <- (b[k] - r) / se[k]
    p <- 2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
  }
  tab <- .coef_table(terms, b, se, stats::qnorm(1 - alpha / 2))
  if (!all(is.finite(c(statistic, p, b, se, tab$conf.low, tab$conf.high)))) {
    stop("Encompassing inference is non-finite; rescale or revise the inputs.",
         call. = FALSE)
  }
  .ft_test(list(
    method = method, direction = direction, nobs = n,
    include_intercept = has_intercept, null_hypothesis = null_label,
    statistic = statistic, distribution = if (joint) "chi2" else "normal",
    df = if (joint) 2L else NA_integer_, df_resid = n - k,
    n_restrictions = if (joint) 2L else 1L, p.value = p, alpha = alpha,
    reject_null = p < alpha, cov_type = "HAC", maxlags = lags, coefficients = tab
  ), paste0("Forecast encompassing test (", method, ", ",
            sub("_encompasses_", " encompasses ", direction), ")"))
}
