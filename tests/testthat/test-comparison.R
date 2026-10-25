set.seed(3)
n <- 70
a <- 50 + cumsum(rnorm(n))
f1 <- a + rnorm(n, 0.1)
f2 <- a + rnorm(n, sd = 1.4)

test_that("MGN slope test, HAC and classical", {
  e1 <- a - f1; e2 <- a - f2; u <- e1 + e2; v <- e1 - e2
  fit <- lm(u ~ v)
  r <- mgn_test(a, f1, f2, maxlags = 3)
  V <- sandwich::NeweyWest(fit, lag = 3, prewhite = FALSE, adjust = FALSE)
  expect_equal(r$statistic, unname(coef(fit)[2] / sqrt(V[2, 2])))
  expect_equal(r$coefficients$estimate, unname(coef(fit)))
  rc <- mgn_test(a, f1, f2, cov_type = "nonrobust")
  rr <- cor(u, v)
  expect_equal(rc$statistic, rr * sqrt((n - 2) / (1 - rr^2)))
  expect_equal(rc$p.value, cor.test(u, v)$p.value)
  expect_equal(rc$correlation, rr)
  expect_equal(rc$mse_difference, mse(a, f1) - mse(a, f2))
  expect_equal(rc$mean_errors, c(model1 = me(a, f1), model2 = me(a, f2)))
})

test_that("MGN: swapping the models flips the sign only", {
  r12 <- mgn_test(a, f1, f2); r21 <- mgn_test(a, f2, f1)
  expect_equal(r12$statistic, -r21$statistic)
  expect_equal(r12$p.value, r21$p.value)
})

test_that("MGN validation", {
  expect_error(mgn_test(1:2, 1:2, 2:3), "at least three")
  expect_error(mgn_test(a, f1, f1), "degenerate")
  expect_error(mgn_test(a, f1, f2, cov_type = "nonrobust", maxlags = 1), "only used")
  expect_error(mgn_test(a, f1, f2[-1]), "equal lengths")
})

test_that("DM equals the HAC t statistic of d on a constant", {
  for (loss in c("squared", "absolute")) for (q in c(0, 2, 5)) {
    r <- dm_test(a, f1, f2, loss = loss, maxlags = q)
    L <- if (loss == "squared") function(e) e^2 else abs
    d <- L(a - f1) - L(a - f2)
    fit <- lm(d ~ 1)
    V <- sandwich::NeweyWest(fit, lag = q, prewhite = FALSE, adjust = FALSE)
    expect_equal(r$mean_loss_diff, unname(coef(fit)))
    expect_equal(r$se, sqrt(V[1, 1]))
    expect_equal(r$dm_statistic, unname(coef(fit) / sqrt(V[1, 1])))
    expect_equal(r$p.value, 2 * pnorm(-abs(r$statistic)))
  }
})

test_that("HLN correction rescales the statistic and uses t(n-1)", {
  r0 <- dm_test(a, f1, f2, h = 3, maxlags = 2)
  r <- dm_test(a, f1, f2, h = 3, maxlags = 2, hln = TRUE)
  k <- sqrt((n + 1 - 2 * 3 + 3 * 2 / n) / n)
  expect_equal(r$hln_factor, k)
  expect_equal(r$statistic, k * r0$dm_statistic)
  expect_equal(r$p.value, 2 * pt(-abs(r$statistic), n - 1))
  expect_equal(r$df, n - 1)
  expect_equal(r$se, r0$se)
  expect_true(is.na(r0$hln_factor))
})

test_that("one-sided alternatives", {
  r <- dm_test(a, f1, f2)
  expect_equal(dm_test(a, f1, f2, alternative = "less")$p.value, pnorm(r$statistic))
  expect_equal(dm_test(a, f1, f2, alternative = "greater")$p.value,
               pnorm(r$statistic, lower.tail = FALSE))
})

test_that("hand-calculated short example", {
  r <- dm_test(c(10, 20, 30, 40, 50), c(9, 22, 27, 39, 48), c(12, 19, 31, 37, 51),
               maxlags = 1)
  expect_equal(r$mean_loss_diff, 0.6)
  expect_identical(r$distribution, "normal")
})

test_that("h does not change the bandwidth or the data", {
  expect_equal(dm_test(a, f1, f2, h = 4)$statistic, dm_test(a, f1, f2)$statistic)
})

test_that("DM validation", {
  expect_error(dm_test(a, f1, f2, loss = "linex"), "loss")
  expect_error(dm_test(a, f1, f2, alternative = "two.sided"), "alternative")
  expect_error(dm_test(a, f1, f2, hln = NA), "hln")
  expect_error(dm_test(a, f1, f2, maxlags = n), "smaller than")
  expect_error(dm_test(1, 1, 2), "at least two")
  expect_error(dm_test(1:3, 1:3, 1:3), "degenerate")
  expect_error(dm_test(a[1:3], f1[1:3], f2[1:3], h = 3, hln = TRUE), "smaller than the sample")
})
