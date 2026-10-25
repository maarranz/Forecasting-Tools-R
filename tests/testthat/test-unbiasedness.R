# Independent oracle: OLS and the Bartlett HAC sandwich written out by hand.
hac_oracle <- function(y, X, L) {
  b <- solve(crossprod(X), crossprod(X, y))
  u <- as.numeric(y - X %*% b)
  g <- X * u
  S <- crossprod(g)
  for (j in seq_len(L)) {
    G <- crossprod(g[(j + 1):nrow(X), , drop = FALSE], g[1:(nrow(X) - j), , drop = FALSE])
    S <- S + (1 - j / (L + 1)) * (G + t(G))
  }
  B <- solve(crossprod(X))
  list(b = as.numeric(b), V = B %*% S %*% B)
}

set.seed(5)
a <- 10 + cumsum(rnorm(40)); f <- a + rnorm(40, 0.2)

test_that("mean-error test equals the hand-computed HAC z test", {
  r <- unbiasedness_test(a, f)
  L <- floor(4 * (40 / 100)^(2 / 9))
  o <- hac_oracle(a - f, matrix(1, 40), L)
  expect_equal(r$maxlags, L)
  expect_equal(r$coefficients$estimate, o$b)
  expect_equal(r$coefficients$estimate, me(a, f))
  expect_equal(r$coefficients$std.error, sqrt(o$V[1, 1]))
  expect_equal(r$statistic, o$b / sqrt(o$V[1, 1]))
  expect_equal(r$p.value, 2 * pnorm(-abs(r$statistic)))
  expect_equal(r$coefficients$conf.low, o$b - qnorm(0.975) * sqrt(o$V[1, 1]))
  expect_identical(r$distribution, "normal")
  expect_identical(r$df, NA_integer_)
})

test_that("Mincer-Zarnowitz: HAC chi-square and classical F", {
  X <- cbind(1, f)
  r <- unbiasedness_test(a, f, "mincer_zarnowitz", maxlags = 2)
  o <- hac_oracle(a, X, 2)
  d <- o$b - c(0, 1)
  expect_equal(r$statistic, as.numeric(t(d) %*% solve(o$V, d)))
  expect_equal(r$p.value, pchisq(r$statistic, 2, lower.tail = FALSE))
  expect_identical(r$coefficients$term, c("intercept", "slope"))
  rc <- unbiasedness_test(a, f, "mincer_zarnowitz", cov_type = "nonrobust")
  fit <- lm(a ~ f)
  expect_equal(rc$coefficients$estimate, unname(coef(fit)))
  expect_equal(rc$coefficients$std.error, unname(sqrt(diag(vcov(fit)))))
  restricted <- sum((a - f)^2); unrestricted <- sum(resid(fit)^2)
  expect_equal(rc$statistic, ((restricted - unrestricted) / 2) / (unrestricted / 38))
  expect_identical(c(rc$df_num, rc$df_denom), c(2L, 38L))
  expect_equal(rc$coefficients$conf.high, unname(confint(fit)[, 2]))
})

test_that("conventional mean-error test is the one-sample t test", {
  r <- unbiasedness_test(a, f, cov_type = "nonrobust")
  tt <- t.test(a - f)
  expect_equal(r$statistic, unname(tt$statistic))
  expect_equal(r$p.value, tt$p.value)
  expect_identical(r$df, 39L)
  expect_true(is.na(r$maxlags))
})

test_that("HAC with zero lags is White's covariance", {
  r <- unbiasedness_test(a, f, "mincer_zarnowitz", maxlags = 0)
  o <- hac_oracle(a, cbind(1, f), 0)
  expect_equal(r$coefficients$std.error, unname(sqrt(diag(o$V))))
})

test_that("hand-calculated small example", {
  r <- unbiasedness_test(c(11, 19, 32, 38), c(10, 20, 30, 40), cov_type = "nonrobust")
  expect_equal(r$coefficients$estimate, 0)
  expect_false(r$reject_null)
})

test_that("invalid options and samples", {
  expect_error(unbiasedness_test(a, f, method = "mz"), "method")
  expect_error(unbiasedness_test(a, f, cov_type = "HC0"), "cov_type")
  expect_error(unbiasedness_test(a, f, alpha = 1), "alpha")
  expect_error(unbiasedness_test(a, f, maxlags = -1), "maxlags")
  expect_error(unbiasedness_test(a, f, maxlags = 1.5), "maxlags")
  expect_error(unbiasedness_test(a, f, maxlags = 40), "smaller than")
  expect_error(unbiasedness_test(a, f, cov_type = "nonrobust", maxlags = 2), "only used")
  expect_error(unbiasedness_test(1, 2), "at least 2")
  expect_error(unbiasedness_test(1:2, 2:3, "mincer_zarnowitz"), "at least 3")
  expect_error(unbiasedness_test(1:5, rep(3, 5), "mincer_zarnowitz"), "rank-deficient")
  expect_error(unbiasedness_test(1:5, 0:4), "degenerate")
  expect_error(unbiasedness_test(1:5, 1:5 + 0, "mincer_zarnowitz"), "degenerate")
})

test_that("result structure and printing", {
  r <- unbiasedness_test(a, f)
  expect_s3_class(r, "forecastingtools_test")
  expect_named(r, c("method", "null_hypothesis", "statistic", "distribution",
                    "p.value", "alpha", "reject_null", "nobs", "cov_type", "maxlags",
                    "coefficients", "n_restrictions", "df_resid", "df", "df_num",
                    "df_denom"))
  expect_output(print(r), "Mean-error unbiasedness test")
  expect_output(print(r), "HAC \\(Bartlett")
})
