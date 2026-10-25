# forecastingtools (R)

Teaching-oriented R tools for forecast evaluation, developed for **Applied Time Series Econometrics (ATSE)** at **Universidad Carlos III de Madrid**. This is the R port of the Python package [Forecasting-Tools](https://github.com/maarranz/Forecasting-Tools): static (sequential one-step) ARIMA forecasting, accuracy measures, and tests of forecast unbiasedness, efficiency and comparison. It complements [arimatools](https://github.com/maarranz/ARIMA-Tools-R), which covers estimation, diagnostics and dynamic forecasting.

## Installation

```r
# install.packages("remotes")
remotes::install_github("maarranz/Forecasting-Tools-R")
library(forecastingtools)
```

With the vignette (requires [pandoc](https://pandoc.org), bundled with RStudio) and the suggested `arimatools`:

```r
remotes::install_github("maarranz/Forecasting-Tools-R", dependencies = TRUE, build_vignettes = TRUE)
vignette("forecasting-tools", package = "forecastingtools")
```

Requires R ≥ 4.1. Imports `forecast` and `sandwich`; `arimatools` is optional (intervention regressors in the vignette and tests).

## Quick start

```r
library(forecast)
library(forecastingtools)

set.seed(20261008)
y <- ts(100 + arima.sim(list(ar = 0.7, ma = 0.4), n = 300, sd = 2))
est <- window(y, end = 250)
eva <- window(y, start = 251)

arma <- static_forecast(Arima(est, order = c(1, 0, 1)), eva)
ar   <- static_forecast(Arima(est, order = c(1, 0, 0)), eva)

forecast_accuracy(eva, arma$forecast)
unbiasedness_test(eva, arma$forecast, method = "mincer_zarnowitz")
weak_efficiency_test(eva, arma$forecast, lags = 10)
dm_test(eva, arma$forecast, ar$forecast, hln = TRUE)
encompassing_test(eva, arma$forecast, ar$forecast)
```

## Public API

The 13 functions of the Python package, with the same names and arguments (except that `static_forecast()` takes the fitted model as `fit` and also accepts `xreg`). Forecast errors are always actual minus forecast.

| Function | Purpose | Returns |
|:--|:--|:--|
| `static_forecast()` | Sequential one-step ARIMA forecasts with fixed parameters and conditional Gaussian intervals | data frame |
| `me()`, `mae()`, `mse()`, `rmse()`, `mape()` | Mean error (bias), mean absolute, mean squared, root mean squared and mean absolute percentage error | number |
| `forecast_accuracy()` | ME, MAE, MSE, RMSE and MAPE together | named vector |
| `unbiasedness_test()` | Mean-error or joint Mincer–Zarnowitz test; HAC or conventional inference | test |
| `weak_efficiency_test()` | Ljung–Box on the errors, or on MA(h−1) residuals for multistep forecasts | test |
| `orthogonality_test()` | Joint test that information available at the origin does not predict the errors | test |
| `mgn_test()` | Morgan–Granger–Newbold covariance restriction (equal MSFE only with zero-mean errors) | test |
| `dm_test()` | Diebold–Mariano test of equal expected squared or absolute loss, optional HLN correction | test |
| `encompassing_test()` | Difference, cross-forecast or unrestricted encompassing, in either direction | test |

Tests return a list of class `forecastingtools_test`, printed as a short report; every element can be extracted (`r$p.value`, `r$coefficients`, ...).

[API documentation — rendered HTML](docs/index.html) · [Quarto source](docs/forecastingtools_Documentation.qmd) · in R: `?dm_test`, etc.

**Vignette:** `vignette("forecasting-tools", package = "forecastingtools")`, the R version of the Python demonstration notebook. It includes the check that the Diebold–Mariano statistic equals the HAC t statistic of the loss differential on a constant (`lm(d ~ 1)` with `sandwich::NeweyWest()`), and static forecasting of a model with a level shift.

## Inference conventions

HAC inference uses `sandwich::NeweyWest(fit, lag = L, prewhite = FALSE, adjust = FALSE)`: Bartlett weights 1 − j/(L+1), no prewhitening and no small-sample factor (Statsmodels' `cov_type="HAC"`, `use_correction=False`). The default bandwidth is min(n − 1, ⌊4(n/100)^(2/9)⌋). Note that `NeweyWest()` prewhitens and adjusts by default, so calling it with its defaults gives different standard errors.

## Agreement with the Python version

`tests/testthat/test-python-parity.R` compares every function with values computed by the Python package on the same data (`data-raw/python_reference.py`). Statistics, p-values, coefficients, standard errors and intervals agree to about 1e-12 or better, and static forecasts and their standard errors to about 1e-14. The exception is `weak_efficiency_test()` with h > 1: the MA(h−1) is estimated by different optimizers (`stats::arima` and Statsmodels), so estimates and statistics agree to about three or four digits.

## Differences from the Python version

| Topic | Python (Forecasting-Tools) | R (forecastingtools) |
|:--|:--|:--|
| Fitted model | nonseasonal Statsmodels `ARIMAResults` | nonseasonal `forecast::Arima()` or `forecast::auto.arima()` |
| Static forecasting mechanism | `results.append(..., refit=False)` | `forecast::Arima(y, model = fit)` (fixed coefficients and variance) |
| Regressors in `static_forecast()` | rejected | supported via `xreg`: interventions (from `arimatools::extend_interventions()`) give ex ante forecasts; other regressors give conditional forecasts, with a message |
| Static forecast output | `Actual`, `Forecast`, `Error`, `SE`, `Lower`, `Upper`, indexed by date | `time`, `actual`, `forecast`, `error`, `se`, `lower`, `upper` |
| Test results | dictionaries | `forecastingtools_test` lists with the same fields; `pvalue` is `p.value`; `parameters`, `standard_errors` and `confidence_intervals` form one `coefficients` data frame; `None` is `NA` |
| Index checks | pandas indexes must match | two `ts` inputs must cover the same periods; other inputs pair by position |
| MA estimation (weak efficiency, h > 1) | Statsmodels state-space ML | `stats::arima(method = "ML")`; reports `loglik`, `aic`, `bic` (no iteration count) |

## Licence

MIT © Miguel A. Arranz
