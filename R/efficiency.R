# Forecast efficiency: weak efficiency and orthogonality -------------------------
#
# Matches weak_efficiency_test() and orthogonality_test() in the Python package
# Forecasting-Tools.

.check_horizon <- function(h) .check_count(h, "h", min = 1)

#' Weak forecast efficiency: serial correlation of forecast errors
#'
#' Under squared-error loss and suitable information assumptions, optimal
#' one-step forecast errors are serially uncorrelated, and optimal h-step
#' errors (from forecasts issued at consecutive origins) follow at most an
#' MA(h - 1).
#'
#' * For `h = 1`, the Ljung-Box test is applied to the forecast errors.
#' * For `h > 1`, an MA(h - 1) with a constant is estimated by exact maximum
#'   likelihood ([stats::arima()], `method = "ML"`), and the Ljung-Box test is
#'   applied to its residuals with `h - 1` degrees of freedom subtracted.
#'
#' The statistic is \eqn{Q = n(n+2) \sum_{j=1}^m r_j^2/(n-j)}, with the
#' centered residual autocorrelations \eqn{r_j}, compared with a chi-square
#' with `lags - (h - 1)` degrees of freedom. It is one joint (portmanteau)
#' test of lags 1 to m, not m separate tests.
#'
#' The test does not examine bias (the errors are centered; use
#' [unbiasedness_test()]). Failure to reject does not establish optimality, and
#' uncorrelated errors need not be independent. The reference distribution is
#' approximate, especially in small samples, with conditional
#' heteroskedasticity or MA parameters near the boundary. Horizons and dates
#' are not inferred from the data: the errors must correspond to consecutive
#' origins.
#'
#' The MA estimates can differ from those of Statsmodels in the last digits,
#' because the optimizers differ; the test statistics agree to the same
#' precision.
#'
#' @inheritParams unbiasedness_test
#' @param h Forecast horizon (positive integer).
#' @param lags Diagnostic lag m: smaller than the number of observations and
#'   greater than `h - 1`. For `h > 1` the sample must also exceed
#'   `(h - 1) + 2` observations.
#' @return A list of class `forecastingtools_test` with elements `method`
#'   (`"weak_efficiency"`), `null_hypothesis`, `h`, `ma_order` (h - 1), `lags`,
#'   `statistic`, `distribution` (`"chi2"`), `p.value`, `df` (lags - ma_order),
#'   `model_df` (ma_order), `alpha`, `reject_null`, `nobs`, `diagnostic_nobs`,
#'   `ma_coefficients` and `ma_standard_errors` (named vectors `ma1`, `ma2`,
#'   ...; empty for h = 1) and `estimation` (`NULL` for h = 1; otherwise a list
#'   with `converged`, `mean`, `sigma2`, `loglik`, `aic`, `bic`, `invertible`,
#'   `residual_burn` and the captured `warnings`).
#' @examples
#' set.seed(41)
#' e <- rnorm(100)
#' weak_efficiency_test(e, rep(0, 100), lags = 5)
#' # Two-step errors: MA(1) dependence is allowed
#' u <- rnorm(121)
#' e2 <- u[-1] + 0.4 * u[-121]
#' weak_efficiency_test(e2, rep(0, 120), h = 2, lags = 8)
#' @export
weak_efficiency_test <- function(actual, forecast, h = 1, lags = 10, alpha = 0.05) {
  h <- .check_horizon(h)
  alpha <- .check_unit(alpha, "alpha")
  lags <- .check_count(lags, "lags", min = 1)
  q <- h - 1L
  if (lags <= q) {
    stop("The diagnostic `lags` must exceed the MA order h - 1 (", q, "), so that ",
         "the Ljung-Box test has positive degrees of freedom.", call. = FALSE)
  }
  s <- .samples(actual = actual, forecast = forecast)
  e <- .errors(s$actual, s$forecast)
  n <- length(e)
  if (lags >= n) {
    stop("The diagnostic `lags` must be smaller than the number of observations (",
         n, ").", call. = FALSE)
  }
  if (q > 0 && n <= q + 2) {
    stop("Estimating the MA(", q, ") requires more than ", q + 2, " observations.",
         call. = FALSE)
  }
  .check_variation(e, "Forecast error")
  diagnostic <- e
  coefs <- stats::setNames(numeric(0), character(0)); ses <- coefs
  estimation <- NULL
  if (q > 0) {
    warns <- character(0)
    fit <- withCallingHandlers(
      tryCatch(stats::arima(e, order = c(0, 0, q), include.mean = TRUE, method = "ML"),
               error = function(err) {
                 stop("MA(", q, ") estimation failed: ", conditionMessage(err),
                      call. = FALSE)
               }),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
      })
    if (fit$code != 0L) {
      stop("MA(", q, ") estimation did not converge; revise the specification or ",
           "the sample.", call. = FALSE)
    }
    cf <- stats::coef(fit)
    V <- fit$var.coef
    se_all <- sqrt(diag(V))
    if (!all(is.finite(cf)) || !all(is.finite(se_all))) {
      stop("MA estimation returned non-finite parameters or standard errors.",
           call. = FALSE)
    }
    if (!(fit$sigma2 > 0)) {
      stop("The estimated MA innovation variance must be positive.", call. = FALSE)
    }
    ma <- paste0("ma", seq_len(q))
    coefs <- cf[ma]; ses <- se_all[ma]
    diagnostic <- as.numeric(stats::residuals(fit))
    k <- q + 2  # MA coefficients, mean and innovation variance
    ll <- fit$loglik
    roots <- polyroot(c(1, coefs))
    estimation <- list(converged = TRUE, mean = unname(cf["intercept"]),
                       sigma2 = fit$sigma2, loglik = ll,
                       aic = -2 * ll + 2 * k, bic = -2 * ll + log(n) * k,
                       invertible = all(Mod(roots) > 1), residual_burn = 0L,
                       warnings = warns)
    if (!all(is.finite(c(ll, estimation$aic, estimation$bic)))) {
      stop("MA estimation returned non-finite likelihood diagnostics.", call. = FALSE)
    }
    .check_variation(diagnostic, "MA residual")
  }
  bt <- stats::Box.test(diagnostic, lag = lags, type = "Ljung-Box", fitdf = q)
  stat <- unname(bt$statistic); p <- bt$p.value
  if (!all(is.finite(c(stat, p)))) {
    stop("The Ljung-Box inference is non-finite; revise the lag or the sample.",
         call. = FALSE)
  }
  .ft_test(list(
    method = "weak_efficiency",
    null_hypothesis = paste0("No autocorrelation of the ",
                             if (q) paste0("MA(", q, ") residuals") else "forecast errors",
                             " at lags 1 to ", lags),
    h = h, ma_order = q, lags = lags, statistic = stat, distribution = "chi2",
    p.value = p, df = lags - q, model_df = q, alpha = alpha, reject_null = p < alpha,
    nobs = n, diagnostic_nobs = length(diagnostic),
    ma_coefficients = coefs, ma_standard_errors = ses, estimation = estimation
  ), "Weak efficiency test (Ljung-Box)")
}

#' Forecast orthogonality: is information available at the origin predictive?
#'
#' Estimates \eqn{e_{t+h|t} = a + \gamma' z_t + u_{t+h}} and tests jointly
#' \eqn{H_0: a = 0, \gamma = 0}, where \eqn{z_t} are variables known when the
#' forecast was made. The intercept is part of the null, so the test includes
#' a zero-bias restriction; separate coefficient tests are not substitutes.
#'
#' **Timing is your responsibility.** The rows of `information` must already
#' correspond to the forecast-error pairs, and every value must have been
#' available at the forecast origin, not merely at the target date: use, for
#' example, `stats::lag()` or `c(y[T0], y[...])` to build \eqn{y_{t-1}} for
#' forecasts of \eqn{y_t}. Nothing is shifted or aligned; look-ahead bias
#' invalidates the test.
#'
#' **Inference.** `cov_type = "HAC"` (default): Newey-West/Bartlett
#' covariance without small-sample adjustment, an asymptotic chi-square Wald
#' test with k restrictions (k = intercept plus information coefficients) and
#' normal intervals. `cov_type = "nonrobust"`: conventional OLS covariance,
#' the classical F(k, n - k) test and t intervals. Overlapping multistep
#' errors are serially correlated: use HAC with at least h - 1 lags. The
#' automatic bandwidth is not increased with h; smaller bandwidths are allowed
#' and flagged in `bandwidth_covers_overlap`.
#'
#' Failure to reject concerns only the supplied variables and a linear form;
#' it does not prove optimality.
#'
#' @inheritParams unbiasedness_test
#' @param information Variables known at the forecast origin: a numeric vector
#'   (one variable), a matrix or a data frame (one column per variable), with
#'   one row per forecast. Column names are kept (unnamed variables become
#'   `information_1`, `information_2`, ...); they must be unique and not
#'   `"intercept"`. Do not include an intercept column. If `information` and
#'   `actual` or `forecast` are `ts` objects, they must cover the same periods.
#' @param h Forecast horizon, recorded as metadata (it does not shift rows).
#' @return A list of class `forecastingtools_test` with elements `method`
#'   (`"orthogonality"`), `null_hypothesis`, `h`, `nobs`, `statistic`,
#'   `distribution` (`"chi2"` or `"F"`), `p.value`, `alpha`, `reject_null`,
#'   `cov_type`, `maxlags`, `n_restrictions` (k), `df_resid` (n - k), `df`,
#'   `df_num`, `df_denom`, `coefficients` (terms `intercept` and the
#'   information names), `information_names` and `bandwidth_covers_overlap`
#'   (`maxlags >= h - 1` for HAC, `NA` for nonrobust).
#' @examples
#' actual <- c(11, 19, 32, 38, 53, 59)
#' forecast <- c(10, 20, 30, 40, 50, 60)
#' orthogonality_test(actual, forecast, cbind(origin_signal = c(0, 1, 0, 1, 0, 1)),
#'                    maxlags = 1)
#' orthogonality_test(actual, forecast, c(0, 1, 0, 1, 0, 1), cov_type = "nonrobust")
#' @export
orthogonality_test <- function(actual, forecast, information, h = 1,
                               cov_type = "HAC", maxlags = NULL, alpha = 0.05) {
  h <- .check_horizon(h)
  alpha <- .check_unit(alpha, "alpha")
  cov_type <- .check_choice(cov_type, c("HAC", "nonrobust"), "cov_type")
  maxlags <- .check_maxlags(maxlags, cov_type)
  s <- .samples(actual = actual, forecast = forecast)
  e <- .errors(s$actual, s$forecast)
  n <- length(e)
  if (!is.null(maxlags)) .bandwidth(maxlags, n)
  .check_same_time(actual, information, "actual", "information")
  .check_same_time(forecast, information, "forecast", "information")
  Z <- .information(information, n)
  names_z <- colnames(Z)
  X <- cbind(1, Z)
  k <- ncol(X)
  if (n <= k) {
    stop("Orthogonality inference requires more observations than regression ",
         "parameters (", k, ").", call. = FALSE)
  }
  if (.rank(X) != k) {
    stop("The orthogonality design is rank-deficient; exclude constant or ",
         "redundant information variables.", call. = FALSE)
  }
  fit <- .ols(e, X)
  if (.degenerate_resid(fit, e)) {
    stop("Residual variance is zero or numerically degenerate; orthogonality ",
         "inference is undefined.", call. = FALSE)
  }
  classical <- cov_type == "nonrobust"
  lags <- if (classical) NA_integer_ else .bandwidth(maxlags, n)
  V <- .check_vcov(.vcov(fit, cov_type, lags), k)
  b <- unname(stats::coef(fit)); se <- sqrt(diag(V))
  w <- .wald(b, V, diag(k), rep(0, k))
  if (classical) {
    statistic <- w / k; p <- stats::pf(statistic, k, n - k, lower.tail = FALSE)
    crit <- stats::qt(1 - alpha / 2, n - k)
  } else {
    statistic <- w; p <- stats::pchisq(statistic, k, lower.tail = FALSE)
    crit <- stats::qnorm(1 - alpha / 2)
  }
  tab <- .coef_table(c("intercept", names_z), b, se, crit)
  if (!all(is.finite(c(statistic, p, tab$conf.low, tab$conf.high)))) {
    stop("Non-finite inference results; rescale the observations.", call. = FALSE)
  }
  .ft_test(list(
    method = "orthogonality",
    null_hypothesis = "intercept = 0 and all information coefficients = 0",
    h = h, nobs = n, statistic = statistic,
    distribution = if (classical) "F" else "chi2", p.value = p, alpha = alpha,
    reject_null = p < alpha, cov_type = cov_type, maxlags = lags,
    n_restrictions = k, df_resid = n - k,
    df = if (classical) NA_integer_ else k,
    df_num = if (classical) k else NA_integer_,
    df_denom = if (classical) n - k else NA_integer_,
    coefficients = tab, information_names = names_z,
    bandwidth_covers_overlap = if (classical) NA else lags >= h - 1L
  ), "Forecast orthogonality test")
}

# Information variables as a named numeric matrix with n rows.
.information <- function(information, n) {
  if (is.data.frame(information)) {
    nms <- names(information)
    cols <- lapply(seq_along(information), function(j) information[[j]])
  } else if (length(dim(information)) == 2L) {
    nms <- colnames(information)
    cols <- lapply(seq_len(ncol(information)), function(j) information[, j])
  } else if (length(dim(information)) > 2L) {
    stop("`information` must be a vector, a matrix or a data frame.", call. = FALSE)
  } else {
    nms <- NULL
    cols <- list(information)
  }
  if (!length(cols)) {
    stop("`information` must contain at least one variable.", call. = FALSE)
  }
  if (is.null(nms) || any(is.na(nms) | nms == "")) {
    nms <- paste0("information_", seq_along(cols))
  }
  if (anyDuplicated(nms) || "intercept" %in% nms) {
    stop("The information variable names must be unique and not \"intercept\".",
         call. = FALSE)
  }
  vals <- lapply(seq_along(cols), function(j) {
    v <- .numeric_sample(cols[[j]], paste0("information (", nms[j], ")"))
    if (length(v) != n) {
      stop("`information` must have the same number of observations as `actual` ",
           "and `forecast` (", n, ").", call. = FALSE)
    }
    v
  })
  Z <- do.call(cbind, vals)
  colnames(Z) <- nms
  Z
}
