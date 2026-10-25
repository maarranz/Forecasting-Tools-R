# Generates the data used to compare forecastingtools with the Python package
# Forecasting-Tools. Run from the package root:
#   Rscript data-raw/reference_data.R
#   python data-raw/python_reference.py /path/to/Forecasting-Tools
# The second step writes tests/testthat/helper-python-reference.R.
# Full double precision, so that both languages read identical numbers.
write_exact <- function(df, file) {
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], sprintf, fmt = "%.17g")
  write.csv(df, file, row.names = FALSE, quote = FALSE)
}

set.seed(20261025)
n <- 60
u <- rnorm(n + 1)
actual <- 50 + cumsum(rnorm(n, sd = 0.3)) + u[-1]
forecast1 <- actual - (0.3 + u[-1] + 0.5 * u[-(n + 1)])   # biased, MA(1) errors
forecast2 <- actual - rnorm(n, sd = 1.3)
information <- cbind(info_a = c(0, head(actual, -1) - 50), info_b = rnorm(n))
write_exact(data.frame(actual, forecast1, forecast2, information),
            "data-raw/reference_samples.csv")

# ARIMA series for static forecasts: estimation 120, evaluation 30
y <- 100 + arima.sim(list(ar = 0.6, ma = 0.3), n = 150, sd = 1.5)
fit1 <- forecast::Arima(y[1:120], order = c(1, 0, 1))
z <- cumsum(0.2 + arima.sim(list(ar = 0.5), n = 150))
fit2 <- forecast::Arima(z[1:120], order = c(1, 1, 1), include.drift = TRUE)
write_exact(data.frame(y = as.numeric(y), z = as.numeric(z)),
            "data-raw/reference_series.csv")
params <- data.frame(
  model = c("arma11", "arima111_drift"),
  ar1 = c(coef(fit1)[["ar1"]], coef(fit2)[["ar1"]]),
  ma1 = c(coef(fit1)[["ma1"]], coef(fit2)[["ma1"]]),
  trend = c(coef(fit1)[["intercept"]], coef(fit2)[["drift"]]),
  sigma2 = c(fit1$sigma2, fit2$sigma2))
write_exact(params, "data-raw/reference_params.csv")
