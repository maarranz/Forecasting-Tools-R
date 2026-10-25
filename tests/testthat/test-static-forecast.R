set.seed(738)
y <- ts(100 + arima.sim(list(ar = 0.55, ma = 0.3), n = 100), start = c(2000, 1),
        frequency = 4)
est <- window(y, end = c(2019, 4))      # 80 observations
eva <- window(y, start = c(2020, 1))    # 20 observations

# Independent oracle: filter the whole series once with the estimated
# parameters; its in-sample one-step predictions are the static forecasts.
oracle <- function(fit, full) {
  refit <- forecast::Arima(full, model = fit)
  as.numeric(stats::fitted(refit))
}

test_that("static forecasts equal one-step predictions of the full-sample filter", {
  for (order in list(c(1, 0, 1), c(1, 1, 1))) {
    fit <- forecast::Arima(est, order = order)
    coefs <- coef(fit)
    sf <- static_forecast(fit, eva, level = 0.9)
    expect_equal(sf$forecast, oracle(fit, y)[81:100], tolerance = 1e-8)
    expect_equal(sf$error, as.numeric(eva) - sf$forecast)
    expect_equal(sf$lower, sf$forecast - qnorm(0.95) * sf$se)
    expect_equal(sf$upper, sf$forecast + qnorm(0.95) * sf$se)
    expect_equal(sf$time, as.numeric(time(eva)))
    expect_equal(attr(sf, "level"), 0.9)
    expect_identical(coef(fit), coefs)                 # fit unchanged
    expect_true(all(sf$se >= sqrt(fit$sigma2) - 1e-12))
  }
})

test_that("first static forecast is the one-step dynamic forecast", {
  fit <- forecast::Arima(est, order = c(1, 0, 1))
  fc <- forecast::forecast(fit, h = 1, level = 95)
  sf <- static_forecast(fit, eva)
  expect_equal(sf$forecast[1], as.numeric(fc$mean))
  expect_equal(sf$upper[1], as.numeric(fc$upper))
})

test_that("later observations do not affect earlier forecasts", {
  fit <- forecast::Arima(est, order = c(1, 0, 1))
  sf <- static_forecast(fit, eva)
  altered <- eva; altered[6:20] <- altered[6:20] + 50
  expect_equal(static_forecast(fit, altered)$forecast[1:6], sf$forecast[1:6])
})

test_that("drift, vectors and auto.arima fits", {
  z <- ts(cumsum(0.3 + arima.sim(list(ar = 0.4), n = 100)), frequency = 4)
  fit <- forecast::Arima(window(z, end = c(20, 4)), order = c(1, 1, 0),
                         include.drift = TRUE)
  sf <- static_forecast(fit, as.numeric(window(z, start = c(21, 1))))
  expect_equal(sf$forecast, oracle(fit, z)[81:100], tolerance = 1e-8)
  fv <- forecast::Arima(as.numeric(est), order = c(1, 0, 1))
  sv <- static_forecast(fv, as.numeric(eva))
  expect_equal(sv$time, 81:100)
  fa <- forecast::auto.arima(est, seasonal = FALSE)
  expect_equal(nrow(static_forecast(fa, eva)), 20L)
})

test_that("models with interventions and other regressors", {
  set.seed(2)
  z <- ts(cumsum(rnorm(120)), start = c(1990, 1), frequency = 4)
  z[60:120] <- z[60:120] + 5
  X <- cbind(LS2004Q4 = as.numeric(seq_len(120) >= 60))
  fit <- forecast::Arima(window(z, end = c(2014, 4)), order = c(1, 1, 0),
                         xreg = X[1:100, , drop = FALSE])
  ev <- window(z, start = c(2015, 1))
  sf <- static_forecast(fit, ev, xreg = X[101:120, , drop = FALSE])
  refit <- forecast::Arima(z, model = fit, xreg = X)
  expect_equal(sf$forecast, as.numeric(fitted(refit))[101:120], tolerance = 1e-8)
  expect_false(attr(sf, "conditional"))
  expect_equal(static_forecast(fit, ev, xreg = rep(1, 20))$forecast, sf$forecast)
  # Rules for the future values of interventions
  expect_error(static_forecast(fit, ev, xreg = cbind(LS2004Q4 = rep(0, 20))),
               "must be 1")
  expect_error(static_forecast(fit, ev), "Supply their values")
  expect_error(static_forecast(fit, ev, xreg = cbind(LS2004Q4 = rep(1, 19))), "19 rows")
  expect_error(static_forecast(fit, ev, xreg = cbind(other = rep(1, 20))), "missing")
  # A stochastic regressor gives conditional forecasts, with a message
  x <- rnorm(120)
  fx <- forecast::Arima(window(z, end = c(2014, 4)), order = c(1, 1, 0),
                        xreg = cbind(x = x[1:100]))
  expect_message(sx <- static_forecast(fx, ev, xreg = cbind(x = x[101:120])),
                 "conditional")
  expect_true(attr(sx, "conditional"))
  refit <- forecast::Arima(z, model = fx, xreg = cbind(x = x))
  expect_equal(sx$forecast, as.numeric(fitted(refit))[101:120], tolerance = 1e-8)
})

test_that("arimatools interventions work directly", {
  skip_if_not_installed("arimatools")
  set.seed(2)
  z <- ts(cumsum(rnorm(120)), start = c(1990, 1), frequency = 4)
  z[60:120] <- z[60:120] + 5
  z[90] <- z[90] + 4
  zest <- window(z, end = c(2014, 4))
  X <- arimatools::outlier_regressors(zest, dates = c("2004Q4", "2012Q2"),
                                      types = c("LS", "TC"), delta = 0.7)
  fit <- forecast::Arima(zest, order = c(1, 1, 0), xreg = X)
  Xf <- arimatools::extend_interventions(X, 20)
  sf <- static_forecast(fit, window(z, start = c(2015, 1)), xreg = Xf)
  expect_false(attr(sf, "conditional"))
  expect_equal(sf$forecast[1], as.numeric(forecast::forecast(fit, h = 1, xreg = Xf[1, , drop = FALSE])$mean))
})

test_that("validation and the empty case", {
  fit <- forecast::Arima(est, order = c(1, 0, 1))
  expect_error(static_forecast(fit, eva, level = 1), "level")
  expect_error(static_forecast(fit, c(1, NA)), "missing")
  expect_error(static_forecast(fit, c(1, Inf)), "non-finite")
  expect_error(static_forecast(fit, "a"), "real numeric")
  expect_error(static_forecast(fit, window(y, start = c(2020, 2))), "immediately after")
  expect_error(static_forecast(fit, ts(1:3, start = c(2020, 1), frequency = 12)),
               "immediately after")
  expect_error(static_forecast(fit, eva, xreg = rep(1, 20)), "not needed")
  expect_error(static_forecast(arima(est, order = c(1, 0, 0)), eva), "forecast::Arima")
  expect_error(static_forecast(lm(est ~ 1), eva), "forecast::Arima")
  for (s in list(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))) {
    fs <- forecast::Arima(est, order = c(1, 0, 0), seasonal = s)
    expect_error(static_forecast(fs, eva), "Seasonal ARIMA models are not supported")
  }
  fl <- forecast::Arima(est, order = c(1, 0, 0), lambda = 0)
  expect_error(static_forecast(fl, eva), "Box-Cox")
  empty <- static_forecast(fit, numeric(0))
  expect_identical(dim(empty), c(0L, 7L))
  expect_named(empty, c("time", "actual", "forecast", "error", "se", "lower", "upper"))
})
