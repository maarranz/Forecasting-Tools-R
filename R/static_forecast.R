# Static (sequential one-step) forecasting -----------------------------------------
#
# Matches static_forecast() in the Python package Forecasting-Tools, and adds
# models with regressors (xreg), which the Python version does not support yet.

# Regressors of a fitted model other than ARMA terms, mean and drift.
.model_regressors <- function(fit) {
  cn <- names(fit$coef)
  cn[!grepl("^s?(ar|ma)[0-9]+$", cn) & !cn %in% c("intercept", "drift")]
}

# Column names created by arimatools::outlier_regressors(): AO, LS or TC
# followed by a period label.
.intervention_pattern <- "^(AO|LS|TC)[0-9]"

.check_static_fit <- function(fit) {
  if (!inherits(fit, "Arima") || is.null(fit$x) || is.null(fit$arma)) {
    stop("`fit` must be a model fitted with forecast::Arima() or ",
         "forecast::auto.arima() (stats::arima() fits do not store the series; ",
         "refit them with forecast::Arima()).", call. = FALSE)
  }
  if (sum(fit$arma[c(3, 4, 7)]) > 0) {
    stop("Seasonal ARIMA models are not supported (as in the Python package ",
         "Forecasting-Tools).", call. = FALSE)
  }
  if (!is.null(fit$lambda)) {
    stop("Models with a Box-Cox transformation (`lambda`) are not supported: the ",
         "forecasts and their standard errors would be on different scales.",
         call. = FALSE)
  }
  invisible(fit)
}

# Evaluation regressors as a numeric matrix with the model's column names.
.evaluation_xreg <- function(fit, xreg, m, actual) {
  regs <- .model_regressors(fit)
  if (!length(regs)) {
    if (!is.null(xreg)) {
      stop("The model has no regressors (other than a mean or drift), so `xreg` is ",
           "not needed.", call. = FALSE)
    }
    return(NULL)
  }
  is_int <- grepl(.intervention_pattern, regs)
  if (is.null(xreg)) {
    stop("The model has regressors (", paste(regs, collapse = ", "), "). Supply ",
         "their values for the ", m, " evaluation periods with `xreg =`",
         if (all(is_int)) ", e.g. `arimatools::extend_interventions(X, length(actual))`",
         ".", call. = FALSE)
  }
  X <- xreg
  if (is.null(dim(X))) {
    if (length(regs) != 1L) {
      stop("`xreg` must be a matrix with columns ", paste(regs, collapse = ", "), ".",
           call. = FALSE)
    }
    X <- matrix(as.numeric(X), ncol = 1, dimnames = list(NULL, regs))
  }
  if (is.null(colnames(X))) {
    stop("`xreg` must have column names: ", paste(regs, collapse = ", "), ".",
         call. = FALSE)
  }
  miss <- setdiff(regs, colnames(X)); extra <- setdiff(colnames(X), regs)
  if (length(miss) || length(extra)) {
    stop("The columns of `xreg` must be the regressors of the model (",
         paste(regs, collapse = ", "), ")",
         if (length(miss)) paste0("; missing: ", paste(miss, collapse = ", ")),
         if (length(extra)) paste0("; not in the model: ", paste(extra, collapse = ", ")),
         ".", call. = FALSE)
  }
  if (nrow(X) != m) {
    stop("`xreg` has ", nrow(X), " rows but `actual` has ", m, " observations.",
         call. = FALSE)
  }
  M <- as.matrix(X)[, regs, drop = FALSE]
  if (!is.numeric(M) || any(!is.finite(M))) {
    stop("`xreg` must contain finite numeric values.", call. = FALSE)
  }
  if (stats::is.ts(xreg)) .check_continues(xreg, stats::as.ts(fit$x), "xreg")
  if (stats::is.ts(actual) && stats::is.ts(xreg) &&
      !isTRUE(all.equal(stats::tsp(actual), stats::tsp(xreg)))) {
    stop("`actual` and `xreg` must cover the same periods.", call. = FALSE)
  }
  M <- matrix(as.numeric(M), nrow = m, dimnames = list(NULL, regs))
  for (j in which(is_int)) {
    v <- M[, j]; type <- substr(regs[j], 1, 2)
    ok <- switch(type,
      AO = all(v == 0),
      LS = all(v == 1),
      TC = all(v > 0 & v < 1) && (length(v) < 2 || {
        g <- diff(log(v))
        all(g < 0) && isTRUE(all.equal(g, rep(g[1], length(g))))
      }))
    if (!ok) {
      rule <- switch(type, AO = "0 (an additive outlier happened once)",
                     LS = "1 (a level shift is permanent)",
                     TC = "a geometric decay continuing from the event date")
      stop("Evaluation values of ", regs[j], " must be ", rule, ". Use ",
           "arimatools::extend_interventions() to generate them. Outliers in the ",
           "evaluation period cannot be added after estimation: that would use ",
           "information from after the forecast origins.", call. = FALSE)
    }
  }
  attr(M, "conditional") <- any(!is_int)
  M
}

.check_continues <- function(x, y, arg) {
  f <- stats::frequency(y)
  expected <- stats::tsp(y)[2] + 1 / f
  if (stats::frequency(x) != f || abs(stats::tsp(x)[1] - expected) > 1e-6 / f) {
    stop("`", arg, "` must start in the period immediately after the end of the ",
         "estimation sample, with the same frequency.", call. = FALSE)
  }
  invisible(TRUE)
}

# Training regressors stored by forecast::Arima(), without the drift (which
# forecast::Arima() regenerates from the time index when `model` is given).
.training_xreg <- function(fit, regs) {
  if (!length(regs)) return(NULL)
  X <- fit$xreg
  if (is.null(X) || !all(regs %in% colnames(X))) {
    stop("The fitted model does not store its regressors; fit it with ",
         "forecast::Arima(..., xreg = ) using a matrix with column names.",
         call. = FALSE)
  }
  matrix(as.numeric(X[, regs, drop = FALSE]), nrow = nrow(X),
         dimnames = list(NULL, regs))
}

#' Static (sequential one-step) forecasts with fixed parameters
#'
#' Evaluates a fitted ARIMA model out of sample by *static* forecasting: the
#' parameters are kept fixed at their estimates, and for each evaluation period
#' the model forecasts the next observation, then that observation is added to
#' the information set before the following forecast. Each actual value is
#' used only after its own forecast. This mimics repeated one-step forecasts as
#' new data arrive, unlike a dynamic forecast issued once at the end of the
#' sample.
#'
#' Each step re-applies the fitted model to the extended series with
#' `forecast::Arima(y, model = fit)`, which re-runs the Kalman filter without
#' re-estimating anything (the coefficients and the innovation variance stay
#' as estimated), and calls [forecast::forecast()] for one step. The
#' intervals are conditional Gaussian prediction intervals: they ignore
#' parameter-estimation uncertainty. A mean, a drift and missing values in the
#' estimation sample are handled by the Kalman filter. Seasonal models and
#' models with a Box-Cox transformation are rejected.
#'
#' **Models with regressors.** The values of the regressors in each evaluation
#' period are supplied in `xreg` (the drift is handled automatically). There
#' are two cases:
#'
#' * Interventions estimated in the estimation sample (AO, LS and TC columns
#'   from `arimatools::outlier_regressors()`) have a known future path: AO is
#'   0, LS is 1 and TC keeps decaying. Generate it with
#'   `arimatools::extend_interventions(X, length(actual))`; values that
#'   contradict the definition are rejected. The forecasts are genuine ex ante
#'   forecasts.
#' * For other regressors, a static forecast of \eqn{y_t} uses the value of
#'   \eqn{x_t}. If \eqn{x_t} is not known at the forecast origin (a
#'   contemporaneous stochastic regressor), the result is a **conditional (ex
#'   post) forecast**: it overstates real-time accuracy, its intervals ignore
#'   the uncertainty about \eqn{x_t}, and efficiency and orthogonality tests no
#'   longer refer to the information at the origin. A message says so and the
#'   result has attribute `conditional = TRUE`. Lagged regressors known at the
#'   origin give ex ante forecasts.
#'
#' @param fit A nonseasonal model fitted with [forecast::Arima()] or
#'   [forecast::auto.arima()] on the estimation sample.
#' @param actual The evaluation observations, in time order: a numeric vector
#'   or a `ts` that starts in the period after the estimation sample (same
#'   frequency). Missing and non-finite values are not allowed. A zero-length
#'   vector gives an empty table.
#' @param level Coverage of the prediction intervals.
#' @param xreg Values of the model's regressors (other than drift) in the
#'   evaluation periods: a matrix or `ts` matrix with one row per element of
#'   `actual` and columns named as in estimation (a vector for a single
#'   regressor). `NULL` for models without regressors.
#' @return A data frame with one row per evaluation period and columns `time`
#'   (as in `time()` of the series), `actual`, `forecast`, `error` (actual minus
#'   forecast), `se`, `lower` and `upper`. Attributes `level` and `conditional`.
#'   The fitted model is not modified.
#' @examples
#' set.seed(20261008)
#' y <- ts(100 + arima.sim(list(ar = 0.7, ma = 0.4), n = 300, sd = 2),
#'         start = c(2000, 1), frequency = 12)
#' est <- window(y, end = c(2020, 10))     # 250 observations
#' eva <- window(y, start = c(2020, 11))   # 50 observations
#' fit <- forecast::Arima(est, order = c(1, 0, 1))
#' sf <- static_forecast(fit, eva)
#' head(sf)
#' forecast_accuracy(sf$actual, sf$forecast)
#'
#' # A model with a level shift in the estimation sample
#' set.seed(2)
#' z <- ts(cumsum(rnorm(120)), start = c(1990, 1), frequency = 4)
#' z[60:120] <- z[60:120] + 5                     # level shift in 2004Q4
#' X <- cbind(LS2004Q4 = as.numeric(seq_len(100) >= 60))
#' fit_ls <- forecast::Arima(window(z, end = c(2014, 4)), order = c(1, 1, 0), xreg = X)
#' head(static_forecast(fit_ls, window(z, start = c(2015, 1)),
#'                      xreg = cbind(LS2004Q4 = rep(1, 20))))
#' @export
static_forecast <- function(fit, actual, level = 0.95, xreg = NULL) {
  .check_static_fit(fit)
  level <- .check_unit(level, "level")
  y <- stats::as.ts(fit$x)   # a plain vector when fitted to a vector
  f <- stats::frequency(y)
  if (stats::is.ts(actual)) .check_continues(actual, y, "actual")
  m <- length(actual)
  empty <- data.frame(time = numeric(0), actual = numeric(0), forecast = numeric(0),
                      error = numeric(0), se = numeric(0), lower = numeric(0),
                      upper = numeric(0))
  if (m == 0L) {
    if (is.data.frame(actual) || length(dim(actual)) > 1L || !is.numeric(actual)) {
      stop("`actual` must be a numeric vector or a univariate ts.", call. = FALSE)
    }
    return(structure(empty, level = level, conditional = FALSE))
  }
  a <- .numeric_sample(actual, "actual")
  Xe <- .evaluation_xreg(fit, xreg, m, actual)
  regs <- colnames(Xe)
  Xt <- .training_xreg(fit, regs)
  n0 <- length(y)
  start <- stats::tsp(y)[1]
  history <- as.numeric(y)
  crit <- stats::qnorm(1 - (1 - level) / 2)
  fc <- se <- numeric(m)
  for (i in seq_len(m)) {
    yi <- stats::ts(c(history, a[seq_len(i - 1)]), start = start, frequency = f)
    if (is.null(Xe)) {
      ri <- forecast::Arima(yi, model = fit)
      out <- forecast::forecast(ri, h = 1, level = 95)
    } else {
      Xi <- rbind(Xt, Xe[seq_len(i - 1), , drop = FALSE])
      ri <- forecast::Arima(yi, model = fit, xreg = Xi)
      out <- forecast::forecast(ri, h = 1, level = 95,
                                xreg = Xe[i, , drop = FALSE])
    }
    fc[i] <- as.numeric(out$mean)
    se[i] <- (as.numeric(out$upper) - fc[i]) / stats::qnorm(0.975)
  }
  conditional <- isTRUE(attr(Xe, "conditional"))
  if (conditional) {
    message("The forecasts use the evaluation-period values of regressors that ",
            "are not interventions (", paste(regs[!grepl(.intervention_pattern, regs)],
                                             collapse = ", "),
            "). If these were not known at the forecast origins, these are ",
            "conditional (ex post) forecasts; see ?static_forecast.")
  }
  tt <- start + (n0 + seq_len(m) - 1) / f
  out <- data.frame(time = tt, actual = a, forecast = fc, error = a - fc, se = se,
                    lower = fc - crit * se, upper = fc + crit * se)
  structure(out, level = level, conditional = conditional)
}
