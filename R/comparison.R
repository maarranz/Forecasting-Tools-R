# Forecast comparison: Morgan-Granger-Newbold and Diebold-Mariano ----------------
#
# Matches mgn_test() and dm_test() in the Python package Forecasting-Tools.

.three <- function(actual, forecast1, forecast2) {
  .samples(actual = actual, forecast1 = forecast1, forecast2 = forecast2)
}

#' Morgan-Granger-Newbold test of equal forecast-error variances
#'
#' With \eqn{e_{1t} = y_t - \hat y_{1t}} and \eqn{e_{2t} = y_t - \hat y_{2t}},
#' define \eqn{u_t = e_{1t} + e_{2t}} and \eqn{v_t = e_{1t} - e_{2t}}. The test
#' estimates \eqn{u_t = a + b v_t + \varepsilon_t} and tests
#' \eqn{H_0: b = 0}, that is, \eqn{\mathrm{Cov}(u, v) = \mathrm{Var}(e_1) -
#' \mathrm{Var}(e_2) = 0}. The intercept is unrestricted.
#'
#' **Interpretation.** Only if *both* forecast errors have zero population
#' means is this the null of equal mean squared forecast errors,
#' \eqn{E(e_1^2) = E(e_2^2)}. Otherwise it tests equality of the centered error
#' variances. Small sample mean errors, or not rejecting unbiasedness, do not
#' prove that assumption. To test equal expected loss directly, use
#' [dm_test()]. The reported MSEs are descriptive only.
#'
#' **Inference.** `cov_type = "HAC"` (default): Newey-West/Bartlett covariance
#' without small-sample adjustment and an asymptotic normal test.
#' `cov_type = "nonrobust"`: the conventional OLS covariance and a Student
#' t(n - 2) test, which reproduces the classical statistic
#' \eqn{r \sqrt{(n-2)/(1-r^2)}} with \eqn{r = \mathrm{Corr}(u, v)}. Swapping the
#' models changes the sign of the statistic, not its two-sided p-value.
#' Perfectly or almost perfectly collinear u and v (\eqn{1 - r^2 \le 10^{-12}})
#' are rejected.
#'
#' @inheritParams unbiasedness_test
#' @param forecast1,forecast2 Forecasts of the two models, of the same length
#'   as `actual` (at least three observations).
#' @return A list of class `forecastingtools_test` with elements `method`
#'   (`"mgn"`), `null_hypothesis` (`"slope = 0"`), `nobs`, `statistic`,
#'   `distribution` (`"normal"` or `"t"`), `df` (`NA` or n - 2), `df_resid`,
#'   `n_restrictions`, `p.value` (two-sided), `alpha`, `reject_null`,
#'   `cov_type`, `maxlags`, `correlation` (of u and v), `coefficients` (terms
#'   `intercept` and `slope`), `mean_errors` and `mse` (named `model1`,
#'   `model2`), `mse_difference` (MSE1 - MSE2; positive means model 2 has the
#'   smaller sample MSE) and `msfe_interpretation_note`.
#' @examples
#' actual <- c(10, 20, 30, 40, 50)
#' f1 <- c(9, 22, 27, 39, 48)
#' f2 <- c(12, 19, 31, 37, 51)
#' mgn_test(actual, f1, f2, cov_type = "nonrobust")
#' mgn_test(actual, f1, f2, maxlags = 1)$p.value
#' @export
mgn_test <- function(actual, forecast1, forecast2, cov_type = "HAC", maxlags = NULL,
                     alpha = 0.05) {
  cov_type <- .check_choice(cov_type, c("HAC", "nonrobust"), "cov_type")
  alpha <- .check_unit(alpha, "alpha")
  maxlags <- .check_maxlags(maxlags, cov_type)
  s <- .three(actual, forecast1, forecast2)
  n <- length(s$actual)
  if (n < 3) {
    stop("MGN regression inference requires at least three observations.",
         call. = FALSE)
  }
  if (!is.null(maxlags)) .bandwidth(maxlags, n)
  e1 <- .errors(s$actual, s$forecast1); e2 <- .errors(s$actual, s$forecast2)
  u <- e1 + e2; v <- e1 - e2
  mse1 <- mean(e1^2); mse2 <- mean(e2^2)
  if (!all(is.finite(c(u, v, mse1, mse2, mse1 - mse2)))) {
    stop("MGN calculations overflow; rescale the observations.", call. = FALSE)
  }
  .check_variation(u, "MGN error sum u")
  .check_variation(v, "MGN error difference v")
  cu <- u - mean(u); cv <- v - mean(v)
  su <- cu / max(abs(cu)); sv <- cv / max(abs(cv))
  r <- sum(su * sv) / (.norm(su) * .norm(sv))
  r <- min(1, max(-1, r))
  if (1 - r^2 <= 1e-12) {
    stop("The MGN sum and difference are perfectly or almost perfectly collinear; ",
         "inference is undefined.", call. = FALSE)
  }
  X <- cbind(1, v)
  if (.rank(X) != 2L) {
    stop("The MGN regression design is numerically rank-deficient.", call. = FALSE)
  }
  fit <- .ols(u, X)
  if (.degenerate_resid(fit, u)) {
    stop("The MGN residual variance is zero or numerically degenerate.", call. = FALSE)
  }
  classical <- cov_type == "nonrobust"
  lags <- if (classical) NA_integer_ else .bandwidth(maxlags, n)
  V <- .check_vcov(.vcov(fit, cov_type, lags), 2L, "MGN parameter")
  b <- unname(stats::coef(fit)); se <- sqrt(diag(V))
  statistic <- b[2] / se[2]
  p <- if (classical) 2 * stats::pt(abs(statistic), n - 2, lower.tail = FALSE)
       else 2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
  crit <- if (classical) stats::qt(1 - alpha / 2, n - 2) else stats::qnorm(1 - alpha / 2)
  tab <- .coef_table(c("intercept", "slope"), b, se, crit)
  if (!all(is.finite(c(statistic, p, b, se, tab$conf.low, tab$conf.high)))) {
    stop("MGN inference results are non-finite; rescale the observations.",
         call. = FALSE)
  }
  .ft_test(list(
    method = "mgn", null_hypothesis = "slope = 0", nobs = n,
    statistic = statistic, distribution = if (classical) "t" else "normal",
    df = if (classical) n - 2L else NA_integer_, df_resid = n - 2L,
    n_restrictions = 1L, p.value = p, alpha = alpha, reject_null = p < alpha,
    cov_type = cov_type, maxlags = lags, correlation = r, coefficients = tab,
    mean_errors = c(model1 = mean(e1), model2 = mean(e2)),
    mse = c(model1 = mse1, model2 = mse2), mse_difference = mse1 - mse2,
    msfe_interpretation_note = paste(
      "Equal-MSFE interpretation requires both errors to have zero population",
      "means; otherwise this tests a covariance (centered-variance) restriction.",
      "Sample MSE summaries are descriptive only.")
  ), "Morgan-Granger-Newbold test")
}

#' Diebold-Mariano test of equal expected forecast loss
#'
#' With losses \eqn{L(e_{it})} and the loss differential
#' \eqn{d_t = L(e_{1t}) - L(e_{2t})}, tests \eqn{H_0: E(d_t) = 0}. Negative
#' mean differentials favour model 1; positive ones favour model 2.
#'
#' The statistic is computed directly as \eqn{\bar d / \sqrt{\hat\Omega/n}},
#' where \eqn{\hat\Omega = \hat\gamma_0 + 2 \sum_{k=1}^q (1 - k/(q+1))
#' \hat\gamma_k} uses the centered autocovariances of \eqn{d_t} with divisor
#' n (Bartlett weights, no small-sample adjustment). This is exactly the HAC t
#' statistic of the regression of \eqn{d_t} on a constant: see the vignette,
#' which checks the equivalence with `lm(d ~ 1)` and [sandwich::NeweyWest()].
#'
#' With `hln = TRUE`, the statistic is multiplied by the
#' Harvey-Leybourne-Newbold factor \eqn{\sqrt{(n + 1 - 2h + h(h-1)/n)/n}} and
#' compared with a Student t(n - 1); the losses, \eqn{\hat\Omega} and the
#' standard error are unchanged. Without it, the reference is the standard
#' normal.
#'
#' Unlike [mgn_test()], the test does not require zero-mean errors. HAC
#' inference needs stationarity, moments and weak dependence, and short
#' samples can make it unreliable; neither version accounts for parameter
#' estimation. Failure to reject does not prove equal accuracy.
#'
#' @inheritParams mgn_test
#' @param loss `"squared"` (\eqn{e^2}) or `"absolute"` (\eqn{|e|}).
#' @param h Forecast horizon (positive integer). Used only by the HLN
#'   correction: it does not shift observations or choose the bandwidth.
#' @param alternative `"two-sided"`, `"less"` (\eqn{E(d) < 0}: model 1 more
#'   accurate) or `"greater"` (\eqn{E(d) > 0}: model 2 more accurate).
#' @param maxlags Bartlett bandwidth q, from 0 to n - 1. `NULL` selects
#'   \eqn{\min(n-1, \lfloor 4 (n/100)^{2/9} \rfloor)}. For overlapping h-step
#'   forecasts, h - 1 is a conventional choice; explicit values are never
#'   overridden.
#' @param hln Apply the Harvey-Leybourne-Newbold correction (requires h < n).
#' @return A list of class `forecastingtools_test` with elements `method`
#'   (`"dm"`), `null_hypothesis`, `n`, `nobs`, `h`, `loss`, `alternative`,
#'   `mean_loss1`, `mean_loss2`, `mean_loss_diff`, `long_run_variance`, `se`,
#'   `dm_statistic` (uncorrected), `statistic` (used for inference), `p.value`,
#'   `maxlags`, `hln`, `hln_factor` (`NA` if not applied), `df` (`NA` or
#'   n - 1), `distribution` (`"normal"` or `"t"`), `alpha` and `reject_null`.
#' @examples
#' dm_test(c(10, 20, 30, 40, 50), c(9, 22, 27, 39, 48), c(12, 19, 31, 37, 51),
#'         maxlags = 1)
#' @export
dm_test <- function(actual, forecast1, forecast2, loss = "squared", h = 1,
                    alternative = "two-sided", maxlags = NULL, hln = FALSE,
                    alpha = 0.05) {
  h <- .check_horizon(h)
  alpha <- .check_unit(alpha, "alpha")
  loss <- .check_choice(loss, c("squared", "absolute"), "loss")
  alternative <- .check_choice(alternative, c("two-sided", "less", "greater"),
                               "alternative")
  hln <- .check_flag(hln, "hln")
  maxlags <- .check_maxlags(maxlags)
  s <- .three(actual, forecast1, forecast2)
  n <- length(s$actual)
  if (n < 2) stop("DM inference requires at least two observations.", call. = FALSE)
  q <- .bandwidth(maxlags, n)
  if (hln && h >= n) {
    stop("The HLN correction requires a horizon h smaller than the sample size.",
         call. = FALSE)
  }
  e1 <- .errors(s$actual, s$forecast1); e2 <- .errors(s$actual, s$forecast2)
  L <- if (loss == "squared") function(e) e^2 else abs
  l1 <- L(e1); l2 <- L(e2)
  d <- l1 - l2
  if (!all(is.finite(d))) {
    stop("DM calculations overflow; rescale the observations.", call. = FALSE)
  }
  dbar <- mean(d)
  centered <- d - dbar
  scale <- max(abs(d))
  if (scale == 0 || max(abs(centered)) / scale <= n * .Machine$double.eps) {
    stop("The loss differential variance is zero or numerically degenerate.",
         call. = FALSE)
  }
  gamma0 <- sum(centered^2) / n
  omega <- gamma0
  for (k in seq_len(q)) {
    omega <- omega + 2 * (1 - k / (q + 1)) *
      sum(centered[(k + 1):n] * centered[1:(n - k)]) / n
  }
  if (!is.finite(omega) || omega <= .Machine$double.eps * gamma0 || gamma0 <= 0) {
    stop("The HAC long-run variance is non-positive, non-finite or numerically ",
         "degenerate.", call. = FALSE)
  }
  se <- sqrt(omega / n)
  dm <- dbar / se
  factor <- if (hln) sqrt((n + 1 - 2 * h + h * (h - 1) / n) / n) else NA_real_
  if (hln && (!is.finite(factor) || factor <= 0)) {
    stop("The HLN correction factor must be finite and strictly positive.",
         call. = FALSE)
  }
  statistic <- if (hln) dm * factor else dm
  cdf <- if (hln) function(x, ...) stats::pt(x, n - 1, ...) else stats::pnorm
  p <- switch(alternative,
              "two-sided" = 2 * cdf(abs(statistic), lower.tail = FALSE),
              less = cdf(statistic),
              greater = cdf(statistic, lower.tail = FALSE))
  if (!all(is.finite(c(dbar, omega, se, dm, statistic, p, mean(l1), mean(l2)))) ||
      se <= 0) {
    stop("DM inference is non-finite or numerically degenerate; rescale the ",
         "observations.", call. = FALSE)
  }
  .ft_test(list(
    method = "dm", null_hypothesis = "Expected loss differential = 0",
    n = n, nobs = n, h = h, loss = loss, alternative = alternative,
    mean_loss1 = mean(l1), mean_loss2 = mean(l2), mean_loss_diff = dbar,
    long_run_variance = omega, se = se, dm_statistic = dm,
    statistic = statistic, p.value = p, maxlags = q, hln = hln,
    hln_factor = factor, df = if (hln) n - 1L else NA_integer_,
    distribution = if (hln) "t" else "normal", alpha = alpha,
    reject_null = p < alpha
  ), paste0("Diebold-Mariano test (", loss, " loss",
            if (hln) ", HLN correction" else "", ")"))
}
