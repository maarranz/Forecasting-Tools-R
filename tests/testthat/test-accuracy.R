test_that("hand-calculated example, for vectors and ts", {
  # Errors 2, -3, 4: sum 3, absolute sum 9, squared sum 29; MAPE terms
  # 20, 15, 40/3 with mean 145/9.
  expected <- c(ME = 1, MAE = 3, MSE = 29 / 3, RMSE = sqrt(29 / 3), MAPE = 145 / 9)
  a <- c(10, 20, 30); f <- c(8, 23, 26)
  expect_equal(forecast_accuracy(a, f), expected)
  expect_equal(c(me(a, f), mae(a, f), mse(a, f), rmse(a, f), mape(a, f)),
               unname(expected))
  expect_equal(forecast_accuracy(ts(a, start = 2000), ts(f, start = 2000)), expected)
  expect_equal(forecast_accuracy(ts(a, start = 2000), f), expected)   # positional
  expect_equal(forecast_accuracy(1:3 * 10L, f), expected)             # integers
})

test_that("signs, perfect forecasts and one observation", {
  expect_equal(forecast_accuracy(c(1, -2, 3), c(1, -2, 3)),
               c(ME = 0, MAE = 0, MSE = 0, RMSE = 0, MAPE = 0))
  expect_equal(me(c(10, 20), c(8, 16)), 3)       # underforecasting
  expect_equal(me(c(10, 20), c(12, 24)), -3)     # overforecasting
  expect_equal(me(c(10, 20), c(8, 22)), 0)       # cancellation
  expect_equal(mae(c(10, 20), c(8, 22)), 2)
  expect_equal(forecast_accuracy(4, 3), c(ME = 1, MAE = 1, MSE = 1, RMSE = 1, MAPE = 25))
})

test_that("MAPE rejects zero actuals; other measures accept them", {
  expect_error(mape(c(0, 2), c(1, 1)), "MAPE is undefined")
  expect_error(forecast_accuracy(c(0, 2), c(1, 1)), "MAPE is undefined")
  expect_equal(mse(c(0, 2), c(0, 2)), 0)
  expect_equal(mape(c(1, 2), c(0, 0)), 100)
})

test_that("invalid samples are rejected, not cleaned", {
  for (fun in list(me, mae, mse, rmse, mape, forecast_accuracy)) {
    expect_error(fun(c(1, 2), 1), "equal lengths")
    expect_error(fun(numeric(0), numeric(0)), "nonempty")
    expect_error(fun(c(1, NA), c(1, 2)), "missing")
    expect_error(fun(c(1, 2), c(1, Inf)), "non-finite")
    expect_error(fun(c("1", "2"), c(1, 2)), "real numeric")
    expect_error(fun(c(TRUE, FALSE), c(1, 2)), "real numeric")
    expect_error(fun(c(1 + 2i, 2), c(1, 2)), "real numeric")
    expect_error(fun(matrix(1:4, 2), 1:4), "one-dimensional")
    expect_error(fun(data.frame(a = 1:2), 1:2), "one-dimensional")
    expect_error(fun(ts(1:3, start = 2000), ts(1:3, start = 2001)), "different periods")
  }
})
