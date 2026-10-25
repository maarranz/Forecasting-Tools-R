# forecastingtools 0.1.0

First version: R port of the Python package Forecasting-Tools 0.5.0.

* `static_forecast()`: sequential one-step forecasts with fixed parameters for
  nonseasonal models fitted with `forecast::Arima()` or
  `forecast::auto.arima()`. Unlike the Python version, it accepts models with
  regressors (`xreg`): interventions from `arimatools::outlier_regressors()`, continued
  with `arimatools::extend_interventions()`, and other regressors, whose
  forecasts are flagged as conditional.
* Accuracy measures: `me()`, `mae()`, `mse()`, `rmse()`, `mape()`,
  `forecast_accuracy()`.
* Tests: `unbiasedness_test()`, `weak_efficiency_test()`,
  `orthogonality_test()`, `mgn_test()`, `dm_test()`, `encompassing_test()`,
  with a print method for their results.
* testthat suite, including a comparison with the Python package on identical
  data, and a vignette that mirrors the Python demonstration notebook.
