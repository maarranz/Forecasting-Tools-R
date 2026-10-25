set.seed(9)
n <- 60
a <- 20 + cumsum(rnorm(n))
f1 <- a + rnorm(n, sd = 0.8)
f2 <- 0.5 * a + 0.5 * f1 + rnorm(n)

nw <- function(fit, L) sandwich::NeweyWest(fit, lag = L, prewhite = FALSE, adjust = FALSE)

test_that("difference regression, both directions, with and without intercept", {
  d <- f1 - f2
  for (dir in c("1_encompasses_2", "2_encompasses_1")) for (ic in c(FALSE, TRUE)) {
    r <- encompassing_test(a, f1, f2, direction = dir, include_intercept = ic, maxlags = 2)
    e <- if (dir == "1_encompasses_2") a - f1 else a - f2
    fit <- if (ic) lm(e ~ d) else lm(e ~ 0 + d)
    k <- length(coef(fit))
    expect_equal(r$coefficients$estimate, unname(coef(fit)))
    expect_equal(r$statistic, unname(coef(fit)[k] / sqrt(nw(fit, 2)[k, k])))
    expect_identical(r$include_intercept, ic)
    expect_identical(tail(r$coefficients$term, 1), "forecast_difference")
  }
})

test_that("cross-forecast and unrestricted regressions", {
  r <- encompassing_test(a, f1, f2, "cross_forecast", "2_encompasses_1")
  fit <- lm((a - f2) ~ 0 + f1)
  expect_equal(r$coefficients$estimate, unname(coef(fit)))
  expect_identical(r$coefficients$term, "competing_forecast")
  ru <- encompassing_test(a, f1, f2, "unrestricted", maxlags = 3)
  fit <- lm(a ~ f1 + f2)
  V <- nw(fit, 3)
  R <- rbind(c(0, 1, 0), c(0, 0, 1)); dd <- R %*% coef(fit) - c(1, 0)
  expect_equal(ru$statistic, as.numeric(t(dd) %*% solve(R %*% V %*% t(R), dd)))
  expect_identical(ru$df, 2L)
  expect_identical(ru$null_hypothesis, "forecast1 = 1 and forecast2 = 0")
})

test_that("unrestricted direction 2 equals direction 1 with the forecasts swapped", {
  r2 <- encompassing_test(a, f1, f2, "unrestricted", "2_encompasses_1")
  r1 <- encompassing_test(a, f2, f1, "unrestricted", "1_encompasses_2")
  expect_equal(r2$statistic, r1$statistic)
  rd2 <- encompassing_test(a, f1, f2, "difference", "2_encompasses_1")
  rd1 <- encompassing_test(a, f2, f1, "difference", "1_encompasses_2")
  expect_equal(rd2$statistic, -rd1$statistic)   # d changes sign
})

test_that("encompassing validation", {
  expect_error(encompassing_test(a, f1, f2, method = "combined"), "method")
  expect_error(encompassing_test(a, f1, f2, direction = "both"), "direction")
  expect_error(encompassing_test(a, f1, f2, "unrestricted", include_intercept = TRUE),
               "always includes")
  expect_error(encompassing_test(a, f1, f1), "rank-deficient")
  expect_error(encompassing_test(a, f1, f1, "unrestricted"), "rank-deficient")
  expect_error(encompassing_test(1:3, 1:3 + 0.5, 3:1, "unrestricted"), "more observations")
  expect_error(encompassing_test(a, f1, f2, maxlags = -1), "maxlags")
})
