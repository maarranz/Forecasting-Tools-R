# Internal helpers -------------------------------------------------------------
#
# None of these functions is exported. They centralize the conventions shared
# by the public API: validation of options and samples, the automatic HAC
# bandwidth, OLS with conventional or Newey-West covariance, Wald tests and the
# result objects. They match the private helpers of the Python package
# Forecasting-Tools (_accuracy_inputs, _efficiency_options, ...).

# ---- Options -----------------------------------------------------------------

.is_flag <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)

.check_flag <- function(x, arg) {
  if (!.is_flag(x)) stop("`", arg, "` must be TRUE or FALSE.", call. = FALSE)
  x
}

.check_choice <- function(x, allowed, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !x %in% allowed) {
    stop("`", arg, "` must be one of: ",
         paste0("\"", allowed, "\"", collapse = ", "), ".", call. = FALSE)
  }
  x
}

.check_unit <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x >= 1) {
    stop("`", arg, "` must be a number strictly between 0 and 1.", call. = FALSE)
  }
  as.numeric(x)
}

.is_count <- function(x, min = 0) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && x == round(x) && x >= min
}

.check_count <- function(x, arg, min = 0) {
  if (!.is_count(x, min)) {
    what <- if (min >= 1) "a positive integer" else "a non-negative integer"
    stop("`", arg, "` must be ", what, ".", call. = FALSE)
  }
  as.integer(x)
}

# Validates `maxlags` before the sample size is known; `.bandwidth()` then
# checks it against n and applies the automatic rule.
.check_maxlags <- function(maxlags, cov_type = "HAC") {
  if (is.null(maxlags)) return(NULL)
  if (!.is_count(maxlags, 0)) {
    stop("`maxlags` must be NULL or a non-negative integer.", call. = FALSE)
  }
  if (cov_type == "nonrobust") {
    stop("`maxlags` is only used with HAC inference; leave it NULL with ",
         "cov_type = \"nonrobust\".", call. = FALSE)
  }
  as.integer(maxlags)
}

# Newey-West automatic bandwidth used by Statsmodels: floor(4 (n/100)^(2/9)),
# capped at n - 1.
.auto_lags <- function(n) as.integer(min(n - 1, floor(4 * (n / 100)^(2 / 9))))

.bandwidth <- function(maxlags, n) {
  if (is.null(maxlags)) return(.auto_lags(n))
  if (maxlags >= n) {
    stop("The HAC bandwidth `maxlags` must be smaller than the number of ",
         "observations (", n, ").", call. = FALSE)
  }
  maxlags
}

# ---- Samples -----------------------------------------------------------------

# One real numeric sample, as a plain vector. Nothing is removed, reordered or
# aligned: missing and non-finite values are errors.
.numeric_sample <- function(x, arg) {
  if (is.data.frame(x) || length(dim(x)) > 1L) {
    stop("`", arg, "` must be a one-dimensional numeric vector or a univariate ts.",
         call. = FALSE)
  }
  if (is.logical(x) || is.character(x) || is.complex(x) || is.factor(x) ||
      inherits(x, c("Date", "POSIXt", "difftime")) || !is.numeric(x)) {
    stop("`", arg, "` must contain real numeric observations (not logical, ",
         "character, complex or dates).", call. = FALSE)
  }
  if (length(x) == 0L) stop("`", arg, "` must be nonempty.", call. = FALSE)
  v <- as.numeric(x)
  if (anyNA(v)) {
    stop("`", arg, "` contains missing observations; all values must be finite.",
         call. = FALSE)
  }
  if (!all(is.finite(v))) {
    stop("`", arg, "` contains non-finite observations; all values must be finite.",
         call. = FALSE)
  }
  v
}

# Two time series must cover the same periods; anything else pairs by position.
.check_same_time <- function(x, y, ax, ay) {
  if (stats::is.ts(x) && stats::is.ts(y) &&
      !isTRUE(all.equal(stats::tsp(x), stats::tsp(y)))) {
    stop("`", ax, "` and `", ay, "` are time series covering different periods; ",
         "their time indexes must match.", call. = FALSE)
  }
  invisible(TRUE)
}

# Validates several samples (named list) and returns them as plain vectors.
.samples <- function(...) {
  args <- list(...)
  nms <- names(args)
  out <- Map(.numeric_sample, args, nms)
  for (i in seq_along(args)[-1]) for (j in seq_len(i - 1)) {
    .check_same_time(args[[j]], args[[i]], nms[j], nms[i])
  }
  n <- lengths(out)
  if (length(unique(n)) != 1L) {
    stop("`", paste(nms, collapse = "`, `"), "` must have equal lengths (",
         paste(n, collapse = ", "), ").", call. = FALSE)
  }
  out
}

.errors <- function(actual, forecast) {
  e <- actual - forecast
  if (!all(is.finite(e))) {
    stop("Forecast errors overflow; rescale the observations.", call. = FALSE)
  }
  e
}

# Absent or numerically negligible variation (as in the Python package).
.check_variation <- function(x, label) {
  if (!all(is.finite(x))) stop(label, " contains non-finite values.", call. = FALSE)
  centered <- x - mean(x)
  if (!all(is.finite(centered)) ||
      .norm(centered) <= length(x) * .Machine$double.eps * max(1, .norm(x))) {
    stop(label, " variance is zero or numerically degenerate.", call. = FALSE)
  }
  invisible(TRUE)
}

# ---- Linear algebra and regression ------------------------------------------

.norm <- function(x) sqrt(sum(x^2))

# Numerical rank with NumPy's default tolerance (matrix_rank).
.rank <- function(X) {
  X <- as.matrix(X)
  if (!all(is.finite(X))) return(NA_integer_)
  s <- svd(X, nu = 0, nv = 0)$d
  if (!length(s)) return(0L)
  as.integer(sum(s > max(s) * max(dim(X)) * .Machine$double.eps))
}

# OLS of y on the columns of X (X contains the intercept column if any).
.ols <- function(y, X) {
  X <- as.matrix(X)
  colnames(X) <- paste0("x", seq_len(ncol(X)))
  stats::lm(y ~ 0 + X)
}

.degenerate_resid <- function(fit, y) {
  r <- stats::residuals(fit)
  .norm(r) <= length(y) * .Machine$double.eps * max(1, .norm(y))
}

# Conventional OLS covariance, or Newey-West HAC with Bartlett weights
# 1 - j/(L+1), no prewhitening and no small-sample adjustment (Statsmodels'
# cov_type="HAC", kernel="bartlett", use_correction=False).
.vcov <- function(fit, cov_type, lags) {
  V <- if (cov_type == "HAC") {
    sandwich::NeweyWest(fit, lag = lags, prewhite = FALSE, adjust = FALSE)
  } else {
    stats::vcov(fit)
  }
  unname(as.matrix(V))
}

.check_vcov <- function(V, k, what = "Parameter") {
  if (!all(is.finite(V)) || .rank(V) != k) {
    stop(what, " covariance is non-finite or singular; inference is undefined.",
         call. = FALSE)
  }
  ok <- tryCatch({ chol(V); TRUE }, error = function(e) FALSE)
  if (!ok) {
    stop(what, " covariance is not positive definite; inference is undefined.",
         call. = FALSE)
  }
  invisible(V)
}

# Wald quadratic form for H0: R b = r.
.wald <- function(b, V, R, r) {
  d <- R %*% b - r
  as.numeric(t(d) %*% solve(R %*% V %*% t(R), d))
}

.coef_table <- function(terms, b, se, crit) {
  data.frame(term = terms, estimate = as.numeric(b), std.error = as.numeric(se),
             conf.low = as.numeric(b - crit * se), conf.high = as.numeric(b + crit * se),
             row.names = NULL, stringsAsFactors = FALSE)
}

# ---- Result objects ----------------------------------------------------------

.ft_test <- function(x, title) structure(x, class = "forecastingtools_test", title = title)

.fmt <- function(x, digits) {
  if (is.null(x) || (length(x) == 1L && is.na(x))) return("-")
  if (is.numeric(x) && !is.integer(x)) format(signif(x, digits)) else format(x)
}

#' Print a forecast test
#'
#' Prints the hypothesis, statistic, p-value and decision of a test from
#' forecastingtools, followed by its coefficient table when there is one.
#' Every element of the list can also be extracted directly, for example
#' `result$p.value` or `result$coefficients`.
#'
#' @param x A result of one of the test functions.
#' @param digits Significant digits.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.forecastingtools_test <- function(x, digits = 4, ...) {
  title <- attr(x, "title")
  cat(title, "\n", strrep("-", nchar(title)), "\n", sep = "")
  dist <- x$distribution
  df <- if (identical(dist, "F")) paste0("(", x$df_num, ", ", x$df_denom, ")")
        else if (!is.null(x$df) && !is.na(x$df)) paste0("(", x$df, ")") else ""
  label <- switch(dist, normal = "z", t = "t", chi2 = "Chi-square", F = "F", dist)
  cat("H0:          ", x$null_hypothesis, "\n", sep = "")
  if (!is.null(x$alternative)) cat("Alternative: ", x$alternative, "\n", sep = "")
  cat("Statistic:   ", .fmt(x$statistic, digits), "  [", label, df, "]\n", sep = "")
  cat("p-value:     ", .fmt(x$p.value, digits), "\n", sep = "")
  cat("Decision:    ", if (isTRUE(x$reject_null)) "reject H0" else "do not reject H0",
      " at ", format(100 * x$alpha), "%\n", sep = "")
  if (!is.null(x$cov_type)) {
    cov <- if (x$cov_type == "HAC") paste0("HAC (Bartlett, ", x$maxlags, " lags)")
           else "conventional OLS"
    cat("Covariance:  ", cov, "\n", sep = "")
  } else if (!is.null(x$maxlags) && !is.na(x$maxlags)) {
    cat("HAC lags:    ", x$maxlags, "\n", sep = "")
  }
  cat("N:           ", x$nobs, "\n", sep = "")
  if (is.data.frame(x$coefficients) && nrow(x$coefficients)) {
    cat("\n")
    tab <- x$coefficients
    num <- vapply(tab, is.numeric, logical(1))
    tab[num] <- lapply(tab[num], signif, digits = digits)
    print(tab, row.names = FALSE)
  }
  invisible(x)
}
