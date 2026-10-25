set.seed(41)
e <- rnorm(100)

test_that("one-step weak efficiency is Box.test on the errors", {
  r <- weak_efficiency_test(e, rep(0, 100), lags = 5)
  bt <- Box.test(e, lag = 5, type = "Ljung-Box")
  expect_equal(r$statistic, unname(bt$statistic))
  expect_equal(r$p.value, bt$p.value)
  expect_identical(c(r$df, r$ma_order, r$model_df), c(5L, 0L, 0L))
  expect_null(r$estimation)
  expect_length(r$ma_coefficients, 0)
  # Manual Ljung-Box with centered errors
  x <- e - mean(e); n <- 100
  rho <- sapply(1:5, function(j) sum(x[(j + 1):n] * x[1:(n - j)]) / sum(x^2))
  expect_equal(r$statistic, n * (n + 2) * sum(rho^2 / (n - 1:5)))
})

test_that("serial correlation raises the statistic", {
  ar <- as.numeric(stats::filter(e, 0.8, method = "recursive"))
  expect_gt(weak_efficiency_test(ar, rep(0, 100))$statistic,
            weak_efficiency_test(e, rep(0, 100))$statistic)
})

test_that("multistep: MA(h-1) residual diagnostic with h-1 degrees of freedom removed", {
  set.seed(7)
  u <- rnorm(301)
  e2 <- u[-1] + 0.4 * u[-301]
  r <- weak_efficiency_test(e2, rep(0, 300), h = 2, lags = 8)
  fit <- arima(e2, order = c(0, 0, 1), method = "ML")
  bt <- Box.test(residuals(fit), lag = 8, type = "Ljung-Box", fitdf = 1)
  expect_equal(r$statistic, unname(bt$statistic))
  expect_equal(r$p.value, bt$p.value)
  expect_identical(r$df, 7L)
  expect_named(r$ma_coefficients, "ma1")
  expect_equal(unname(r$ma_coefficients), 0.4, tolerance = 0.3)
  expect_true(r$estimation$converged)
  expect_true(r$estimation$invertible)
  expect_equal(r$estimation$aic, -2 * fit$loglik + 2 * 3)
  expect_identical(r$diagnostic_nobs, 300L)
})

test_that("weak efficiency validation", {
  expect_error(weak_efficiency_test(e, rep(0, 100), h = 0), "`h`")
  expect_error(weak_efficiency_test(e, rep(0, 100), lags = 0), "`lags`")
  expect_error(weak_efficiency_test(e, rep(0, 100), h = 3, lags = 2), "exceed the MA order")
  expect_error(weak_efficiency_test(e[1:5], rep(0, 5), lags = 5), "smaller than")
  expect_error(weak_efficiency_test(e[1:4], rep(0, 4), h = 3, lags = 3), "more than 4")
  expect_error(weak_efficiency_test(rep(1, 20), rep(0, 20)), "degenerate")
})

set.seed(11)
n <- 80
a <- cumsum(rnorm(n)); f <- c(0, head(a, -1)) + rnorm(n, sd = 0.3)
z <- cbind(prev = c(0, head(a, -1)), other = rnorm(n))

test_that("orthogonality: HAC Wald test of all coefficients, and classical F", {
  r <- orthogonality_test(a, f, z, maxlags = 2)
  fit <- lm((a - f) ~ z)
  V <- sandwich::NeweyWest(fit, lag = 2, prewhite = FALSE, adjust = FALSE)
  b <- coef(fit)
  expect_equal(r$coefficients$estimate, unname(b))
  expect_equal(r$coefficients$std.error, unname(sqrt(diag(V))))
  expect_equal(r$statistic, as.numeric(t(b) %*% solve(V, b)))
  expect_identical(r$coefficients$term, c("intercept", "prev", "other"))
  expect_identical(r$information_names, c("prev", "other"))
  expect_identical(r$df, 3L)
  rc <- orthogonality_test(a, f, z, cov_type = "nonrobust")
  Vc <- vcov(fit)
  expect_equal(rc$statistic, as.numeric(t(b) %*% solve(Vc, b)) / 3)
  expect_identical(c(rc$df_num, rc$df_denom), c(3L, 77L))
  expect_true(is.na(rc$bandwidth_covers_overlap))
})

test_that("orthogonality accepts vectors, matrices and data frames", {
  v <- orthogonality_test(a, f, z[, 1])
  m <- orthogonality_test(a, f, z[, 1, drop = FALSE])
  df <- orthogonality_test(a, f, data.frame(prev = z[, 1]))
  expect_identical(v$information_names, "information_1")
  expect_identical(m$information_names, "prev")
  expect_equal(v$statistic, m$statistic)
  expect_equal(m$statistic, df$statistic)
  expect_false(orthogonality_test(a, f, z, h = 4, maxlags = 2)$bandwidth_covers_overlap)
  expect_true(orthogonality_test(a, f, z, h = 3, maxlags = 2)$bandwidth_covers_overlap)
})

test_that("orthogonality validation", {
  expect_error(orthogonality_test(a, f, z[-1, ]), "same number")
  expect_error(orthogonality_test(a, f, cbind(z, z)), "unique")
  expect_error(orthogonality_test(a, f, cbind(intercept = 1:n)), "unique")
  expect_error(orthogonality_test(a, f, cbind(k = rep(1, n))), "rank-deficient")
  expect_error(orthogonality_test(a, f, cbind(z, z2 = 2 * z[, 1])), "rank-deficient")
  expect_error(orthogonality_test(a, f, c(z[-1, 1], NA)), "missing")
  expect_error(orthogonality_test(a[1:2], f[1:2], 1:2), "more observations")
  expect_error(orthogonality_test(ts(a, start = 1), f, ts(z[, 1], start = 2)),
               "different periods")
})
