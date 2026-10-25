# Forecast accuracy measures ----------------------------------------------------
#
# Matches me(), mae(), mse(), rmse(), mape() and forecast_accuracy() in the
# Python package Forecasting-Tools. Errors are actual minus forecast; the
# denominator is always n.

.mape <- function(actual, forecast) {
  if (any(actual == 0)) {
    stop("MAPE is undefined when any actual observation is zero.", call. = FALSE)
  }
  100 * mean(abs((actual - forecast) / actual))
}

.acc <- function(actual, forecast) {
  s <- .samples(actual = actual, forecast = forecast)
  list(actual = s$actual, forecast = s$forecast, e = .errors(s$actual, s$forecast))
}

#' Forecast accuracy measures
#'
#' With forecast errors \eqn{e_i = y_i - \hat y_i}, \eqn{i = 1, \dots, n}:
#' \deqn{\mathrm{ME} = \frac1n \sum e_i, \quad
#'       \mathrm{MAE} = \frac1n \sum |e_i|, \quad
#'       \mathrm{MSE} = \frac1n \sum e_i^2,}
#' \deqn{\mathrm{RMSE} = \sqrt{\mathrm{MSE}}, \quad
#'       \mathrm{MAPE} = \frac{100}{n} \sum \left| \frac{e_i}{y_i} \right|.}
#'
#' **ME** measures signed bias: a positive value means underforecasting on
#' average, a negative value overforecasting. Positive and negative errors
#' cancel, so a small ME does not mean accurate forecasts. **MAE**, **MSE** and
#' **RMSE** measure the size of the errors; MSE and RMSE give more weight to
#' large errors (and always rank models in the same way). All four are
#' scale-dependent. **MAPE** is a relative measure, in percent; it is undefined
#' when an actual value is zero (an error is raised), unstable near zero, and
#' most meaningful for positive data with a natural zero.
#'
#' The denominator is n, with no degrees-of-freedom adjustment. Inputs are
#' paired by position; nothing is removed or aligned, so missing or non-finite
#' values are errors. If both inputs are `ts` objects, they must cover the same
#' periods.
#'
#' @param actual Actual values: a numeric vector or univariate `ts`.
#' @param forecast Forecasts of the same length.
#' @return A number (in percent for `mape()`). `forecast_accuracy()` returns a
#'   named numeric vector with elements `ME`, `MAE`, `MSE`, `RMSE` and `MAPE`,
#'   in that order.
#' @examples
#' actual <- c(10, 20, 30)
#' forecast <- c(8, 23, 26)          # errors 2, -3, 4
#' me(actual, forecast)              # 1
#' mae(actual, forecast)             # 3
#' mse(actual, forecast)             # 29/3
#' forecast_accuracy(actual, forecast)
#' @name accuracy
NULL

#' @rdname accuracy
#' @export
me <- function(actual, forecast) mean(.acc(actual, forecast)$e)

#' @rdname accuracy
#' @export
mae <- function(actual, forecast) mean(abs(.acc(actual, forecast)$e))

#' @rdname accuracy
#' @export
mse <- function(actual, forecast) mean(.acc(actual, forecast)$e^2)

#' @rdname accuracy
#' @export
rmse <- function(actual, forecast) sqrt(mse(actual, forecast))

#' @rdname accuracy
#' @export
mape <- function(actual, forecast) {
  a <- .acc(actual, forecast)
  .mape(a$actual, a$forecast)
}

#' @rdname accuracy
#' @export
forecast_accuracy <- function(actual, forecast) {
  a <- .acc(actual, forecast)
  pct <- .mape(a$actual, a$forecast)
  m2 <- mean(a$e^2)
  c(ME = mean(a$e), MAE = mean(abs(a$e)), MSE = m2, RMSE = sqrt(m2), MAPE = pct)
}
