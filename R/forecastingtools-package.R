#' forecastingtools: teaching-oriented tools for forecast evaluation
#'
#' R port of the Python package Forecasting-Tools (ATSE, UC3M).
#' [static_forecast()] produces sequential one-step forecasts from a model
#' fitted with [forecast::Arima()] or [forecast::auto.arima()]. Every other
#' function takes plain numeric vectors (or `ts` objects) of actual values and
#' forecasts, so it can evaluate forecasts from any source.
#'
#' Forecast errors are always actual minus forecast. Function names and
#' arguments match the Python API. Tests return a list of class
#' `forecastingtools_test` with the same fields as the Python dictionaries
#' (`pvalue` is called `p.value`, and the coefficient estimates, standard
#' errors and confidence intervals are collected in one `coefficients` data
#' frame). HAC inference uses [sandwich::NeweyWest()] with Bartlett weights,
#' no prewhitening and no small-sample adjustment, which reproduces the
#' Statsmodels results of the Python package.
#'
#' Part of the ATSE software suite, with arimatools (ARIMA-Tools-R),
#' ARIMA-Tools and Forecasting-Tools (Python).
#'
#' @keywords internal
"_PACKAGE"
