# Agreement with the Python package Forecasting-Tools on identical data.
# Reference values: helper-python-reference.R, written by
# data-raw/python_reference.py. Tolerances are relative.

d <- py_ref$data
y <- d$actual; f1 <- d$forecast1; f2 <- d$forecast2
tol <- 1e-9

expect_matches_python <- function(r, p, tolerance = tol) {
  expect_equal(r$statistic, p$statistic, tolerance = tolerance)
  expect_equal(r$p.value, p$pvalue, tolerance = tolerance)
  expect_identical(r$reject_null, p$reject_null)
  expect_identical(r$distribution, p$distribution)
  expect_identical(r$nobs, p$nobs)
  for (f in c("df", "df_num", "df_denom", "n_restrictions", "df_resid", "maxlags")) {
    if (f %in% names(p)) {
      if (is.null(p[[f]])) expect_true(is.na(r[[f]]), label = f)
      else expect_equal(r[[f]], p[[f]], label = f)
    }
  }
  if (!is.null(p$parameters)) {
    tab <- r$coefficients
    expect_identical(tab$term, names(p$parameters))
    expect_equal(tab$estimate, unname(unlist(p$parameters)), tolerance = tolerance)
    expect_equal(tab$std.error, unname(unlist(p$standard_errors)), tolerance = tolerance)
    ci <- do.call(rbind, p$confidence_intervals)
    expect_equal(cbind(tab$conf.low, tab$conf.high), unname(ci), tolerance = tolerance)
  }
}

test_that("accuracy measures match Python", {
  expect_equal(forecast_accuracy(y, f1), unlist(py_ref$accuracy), tolerance = tol)
})

test_that("unbiasedness tests match Python", {
  for (m in c("mean_error", "mincer_zarnowitz")) {
    expect_matches_python(unbiasedness_test(y, f1, m, "HAC"),
                          py_ref$unbiasedness[[paste0(m, "_HAC")]])
    expect_matches_python(unbiasedness_test(y, f1, m, "nonrobust"),
                          py_ref$unbiasedness[[paste0(m, "_nonrobust")]])
    expect_matches_python(unbiasedness_test(y, f1, m, maxlags = 3),
                          py_ref$unbiasedness[[paste0(m, "_HAC3")]])
  }
})

test_that("weak efficiency matches Python (MA estimates up to optimizer precision)", {
  expect_matches_python(weak_efficiency_test(y, f1, h = 1, lags = 10), py_ref$weak$h1)
  for (case in list(list("h2", f1, 2L), list("h3", f2, 3L))) {
    p <- py_ref$weak[[case[[1]]]]
    r <- weak_efficiency_test(y, case[[2]], h = case[[3]], lags = 8)
    expect_equal(r$statistic, p$statistic, tolerance = 1e-3)
    expect_equal(r$p.value, p$pvalue, tolerance = 1e-3)
    expect_equal(unname(r$ma_coefficients), unname(unlist(p$ma_coefficients)),
                 tolerance = 1e-3)
    expect_identical(r$df, p$df)
    expect_identical(r$diagnostic_nobs, p$diagnostic_nobs)
  }
})

test_that("orthogonality tests match Python", {
  expect_matches_python(orthogonality_test(y, f1, d$info_a), py_ref$orthogonality$a_HAC)
  expect_matches_python(orthogonality_test(y, f1, d$info_a, cov_type = "nonrobust"),
                        py_ref$orthogonality$a_nonrobust)
  r <- orthogonality_test(y, f1, data.frame(info_a = d$info_a, info_b = d$info_b),
                          h = 3, maxlags = 2)
  expect_matches_python(r, py_ref$orthogonality$ab_HAC2)
  expect_identical(r$information_names, py_ref$orthogonality$ab_HAC2$information_names)
  expect_identical(r$bandwidth_covers_overlap,
                   py_ref$orthogonality$ab_HAC2$bandwidth_covers_overlap)
})

test_that("MGN tests match Python", {
  for (k in c("HAC", "nonrobust")) {
    r <- mgn_test(y, f1, f2, cov_type = k); p <- py_ref$mgn[[k]]
    expect_matches_python(r, p)
    expect_equal(r$correlation, p$correlation, tolerance = tol)
    expect_equal(r$mse_difference, p$mse_difference, tolerance = tol)
  }
  expect_matches_python(mgn_test(y, f1, f2, maxlags = 0), py_ref$mgn$HAC0)
})

test_that("DM tests match Python", {
  for (k in names(py_ref$dm)) {
    p <- py_ref$dm[[k]]
    r <- dm_test(y, f1, f2, loss = p$loss, h = p$h, alternative = p$alternative,
                 maxlags = p$maxlags, hln = p$hln)
    expect_matches_python(r, p)
    for (f in c("mean_loss1", "mean_loss2", "mean_loss_diff", "long_run_variance",
                "se", "dm_statistic")) {
      expect_equal(r[[f]], p[[f]], tolerance = tol, label = paste(k, f))
    }
    if (p$hln) expect_equal(r$hln_factor, p$hln_factor, tolerance = tol)
  }
})

test_that("encompassing tests match Python", {
  for (k in names(py_ref$encompassing)) {
    p <- py_ref$encompassing[[k]]
    r <- encompassing_test(y, f1, f2, method = p$method, direction = p$direction,
                           include_intercept = grepl("_True_", k) ,
                           maxlags = if (grepl("_2$", k)) 2 else NULL)
    expect_matches_python(r, p, tolerance = 1e-8)
    expect_identical(r$null_hypothesis, p$null_hypothesis)
    expect_identical(r$include_intercept, p$include_intercept)
  }
})

test_that("static forecasts match Python's append(refit=False) loop", {
  s <- py_ref$series
  # The models are refitted with the coefficients and innovation variance
  # given to Python, so that the comparison does not depend on the optimizer.
  fixed_fit <- function(x, order, p, drift) {
    fit <- suppressWarnings(forecast::Arima(
      x, order = order, include.drift = drift, include.mean = !drift,
      fixed = c(p$ar1, p$ma1, p$trend)))
    fit$sigma2 <- p$sigma2
    fit
  }
  fits <- list(
    arma11 = fixed_fit(s$y[1:120], c(1, 0, 1), py_ref$params$arma11, FALSE),
    arima111_drift = fixed_fit(s$z[1:120], c(1, 1, 1), py_ref$params$arima111_drift, TRUE))
  for (k in names(fits)) {
    a <- (if (k == "arma11") s$y else s$z)[121:150]
    sf <- static_forecast(fits[[k]], a, level = 0.9)
    p <- py_ref$static[[k]]
    for (col in c("actual", "forecast", "error", "se", "lower", "upper")) {
      expect_equal(sf[[col]], p[[col]], tolerance = 1e-9, label = paste(k, col))
    }
  }
})
