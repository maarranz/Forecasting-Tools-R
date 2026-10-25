# Forecast unbiasedness -----------------------------------------------------------
#
# Matches unbiasedness_test() in the Python package Forecasting-Tools.

#' Forecast unbiasedness: mean-error and Mincer-Zarnowitz tests
#'
#' **Mean error** (`method = "mean_error"`): estimates \eqn{e_t = a + u_t}, with
#' \eqn{e_t = y_t - \hat y_t}, and tests \eqn{H_0: a = 0}. The estimate of
#' \eqn{a} is the mean error: positive values indicate underforecasting,
#' negative values overforecasting. At least two observations are required.
#'
#' **Mincer-Zarnowitz** (`method = "mincer_zarnowitz"`): estimates
#' \eqn{y_t = a + b \hat y_t + u_t} and tests the joint restriction
#' \eqn{H_0: a = 0, b = 1}. This is one joint test, not two separate
#' coefficient tests. At least three observations and varying forecasts are
#' required.
#'
#' **Inference.** With `cov_type = "HAC"` (default) the covariance is
#' Newey-West with Bartlett weights \eqn{1 - j/(L+1)} and no small-sample
#' adjustment ([sandwich::NeweyWest()] with `prewhite = FALSE, adjust = FALSE`);
#' the mean-error test uses an asymptotic normal z statistic, the joint test an
#' asymptotic chi-square with 2 degrees of freedom, and the confidence
#' intervals normal critical values. With `cov_type = "nonrobust"` the
#' covariance is the conventional OLS one, with Student t inference and the
#' classical F(2, n - 2) joint test.
#'
#' Rejection is evidence against the null; failure to reject does not prove
#' unbiasedness, since small samples can have little power. Zero mean errors
#' alone do not establish calibration, and neither test establishes forecast
#' optimality. The observations must be consecutive periods (HAC uses their
#' order); `ts` attributes are not used to detect gaps.
#'
#' @param actual Actual values: a numeric vector or univariate `ts`.
#' @param forecast Forecasts of the same length.
#' @param method `"mean_error"` or `"mincer_zarnowitz"`.
#' @param cov_type `"HAC"` or `"nonrobust"`.
#' @param maxlags HAC bandwidth L: `NULL` (default) selects
#'   \eqn{\min(n-1, \lfloor 4 (n/100)^{2/9} \rfloor)}; otherwise a non-negative
#'   integer smaller than n. Must be `NULL` with `cov_type = "nonrobust"`.
#' @param alpha Significance level; intervals have coverage `1 - alpha` and
#'   `reject_null` is `p.value < alpha`.
#' @return A list of class `forecastingtools_test` with elements `method`,
#'   `null_hypothesis`, `statistic`, `distribution` (`"normal"`, `"t"`,
#'   `"chi2"` or `"F"`), `p.value`, `alpha`, `reject_null`, `nobs`,
#'   `cov_type`, `maxlags` (`NA` for nonrobust), `coefficients` (a data frame
#'   with `term`, `estimate`, `std.error`, `conf.low`, `conf.high`; terms
#'   `intercept` and, for Mincer-Zarnowitz, `slope`), `n_restrictions`,
#'   `df_resid`, `df` (t degrees of freedom or chi-square restrictions, else
#'   `NA`), `df_num` and `df_denom` (F test only, else `NA`).
#' @examples
#' unbiasedness_test(c(11, 19, 32, 38), c(10, 20, 30, 40), cov_type = "nonrobust")
#' unbiasedness_test(c(11, 19, 32, 38), c(10, 20, 30, 40),
#'                   method = "mincer_zarnowitz", maxlags = 1)
#' @export
unbiasedness_test <- function(actual, forecast, method = "mean_error",
                              cov_type = "HAC", maxlags = NULL, alpha = 0.05) {
  method <- .check_choice(method, c("mean_error", "mincer_zarnowitz"), "method")
  cov_type <- .check_choice(cov_type, c("HAC", "nonrobust"), "cov_type")
  alpha <- .check_unit(alpha, "alpha")
  maxlags <- .check_maxlags(maxlags, cov_type)
  s <- .samples(actual = actual, forecast = forecast)
  n <- length(s$actual)
  joint <- method == "mincer_zarnowitz"
  k <- if (joint) 2L else 1L
  if (n <= k) {
    stop(method, " requires at least ", k + 1, " observations.", call. = FALSE)
  }
  if (!is.null(maxlags)) .bandwidth(maxlags, n)
  X <- if (joint) cbind(1, s$forecast) else matrix(1, n, 1)
  if (.rank(X) != k) {
    stop("The Mincer-Zarnowitz design is rank-deficient: the forecasts must vary.",
         call. = FALSE)
  }
  y <- if (joint) s$actual else .errors(s$actual, s$forecast)
  fit <- .ols(y, X)
  if (.degenerate_resid(fit, y)) {
    stop("Residual variance is zero or numerically degenerate; unbiasedness ",
         "inference is undefined.", call. = FALSE)
  }
  lags <- if (cov_type == "HAC") .bandwidth(maxlags, n) else NA_integer_
  V <- .check_vcov(.vcov(fit, cov_type, lags), k)
  b <- unname(stats::coef(fit)); se <- sqrt(diag(V))
  df_resid <- n - k
  classical <- cov_type == "nonrobust"
  if (joint) {
    w <- .wald(b, V, diag(2), c(0, 1))
    if (classical) {
      statistic <- w / 2; dist <- "F"
      p <- stats::pf(statistic, 2, df_resid, lower.tail = FALSE)
    } else {
      statistic <- w; dist <- "chi2"
      p <- stats::pchisq(statistic, 2, lower.tail = FALSE)
    }
  } else {
    statistic <- b / se
    dist <- if (classical) "t" else "normal"
    p <- if (classical) 2 * stats::pt(abs(statistic), df_resid, lower.tail = FALSE)
         else 2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
  }
  crit <- if (classical) stats::qt(1 - alpha / 2, df_resid) else stats::qnorm(1 - alpha / 2)
  tab <- .coef_table(if (joint) c("intercept", "slope") else "intercept", b, se, crit)
  if (!all(is.finite(c(statistic, p, tab$conf.low, tab$conf.high)))) {
    stop("Non-finite inference results; rescale the observations.", call. = FALSE)
  }
  .ft_test(list(
    method = method,
    null_hypothesis = if (joint) "intercept = 0, slope = 1" else "mean_error = 0",
    statistic = statistic, distribution = dist, p.value = p, alpha = alpha,
    reject_null = p < alpha, nobs = n, cov_type = cov_type, maxlags = lags,
    coefficients = tab, n_restrictions = k, df_resid = df_resid,
    df = if (dist == "t") df_resid else if (dist == "chi2") k else NA_integer_,
    df_num = if (dist == "F") k else NA_integer_,
    df_denom = if (dist == "F") df_resid else NA_integer_
  ), if (joint) "Mincer-Zarnowitz unbiasedness test" else "Mean-error unbiasedness test")
}
