library(testthat)

# Helper: run expr and re-emit ONLY the warning matching `pattern`;
# all other warnings are muffled so test output stays clean.
emit_only_warning <- function(expr, pattern) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl(pattern, msg)) {
        # Re-emit your target warning for expect_warning to catch:
        warning(msg, call. = FALSE)
      }
      # Muffle all warnings so nothing else escapes:
      invokeRestart("muffleWarning")
    }
  )
}

test_that("log1mexp computes correct values for basic inputs", {
  # Known identity: log1mexp(log(3)) == log(2/3)
  expect_equal(log1mexp(log(3)), log(2 / 3), tolerance = 1e-15)

  # At x = 0: log(1 - exp(0)) = log(0) = -Inf
  expect_identical(log1mexp(0), -Inf)

  # At large positive x: value is close to 0 from below
  v <- log1mexp(50)
  expect_lte(v, 0)
  expect_lt(abs(v), 1e-20)
})

test_that("log1mexp uses stable branches around the cutoff", {
  cutoff <- log(2)

  # Just below and just above cutoff to hit each branch
  x1 <- cutoff * (1 - 1e-12) # left branch: log(-expm1(-x))
  x2 <- cutoff * (1 + 1e-12) # right branch: log1p(-exp(-x))

  f1 <- log1mexp(x1)
  f2 <- log1mexp(x2)

  ref1 <- log(-expm1(-x1))
  ref2 <- log1p(-exp(-x2))

  expect_equal(f1, ref1, tolerance = 1e-14)
  expect_equal(f2, ref2, tolerance = 1e-14)
})

test_that("log1mexp is vectorized and positions align", {
  x <- c(log(3), 0, 10, 1e-8)
  res <- log1mexp(x)

  # Element-wise checks
  expect_equal(res[1], log(2 / 3), tolerance = 1e-15)
  expect_identical(res[2], -Inf)

  # For x = 10, compare to a stable reference and assert small magnitude
  expect_equal(res[3], log1p(-exp(-10)), tolerance = 1e-15)
  expect_lte(res[3], 0)
  expect_lt(abs(res[3]), 1e-4) # realistic cross-platform tolerance

  # For very small positive x, prefer the expm1-stable expression
  expect_equal(res[4], log(-expm1(-1e-8)), tolerance = 1e-15)
})

test_that("log1mexp preserves NAs and names", {
  x <- c(a = log(3), b = 0, c = NA_real_)
  res <- log1mexp(x)

  # NA preserved at the same position/name
  expect_true(is.na(res["c"]))
  expect_identical(names(res), names(x))

  # Single-element extraction with '[' keeps the name; compare numerically
  expect_identical(unname(res["b"]), -Inf)
})

test_that("log1mexp preserves names even when there are no NAs", {
  x <- c(a = log(3), b = 0)
  res <- log1mexp(x)

  # Values correct
  expect_equal(unname(res), c(log(2 / 3), -Inf), tolerance = 1e-15)

  # Names are preserved by current implementation (r <- a)
  expect_identical(names(res), names(x))
})

test_that("log1mexp warns for negative inputs (current behavior)", {
  # Only emit the "'a' >= 0 needed" warning; suppress math NaNs warnings
  expect_warning(
    res <- emit_only_warning(log1mexp(c(-1, 0)), "'a' >= 0 needed"),
    "'a' >= 0 needed"
  )
  # Spot check: first is NaN (negative), second is -Inf (zero)
  expect_true(is.nan(res[1]))
  expect_identical(res[2], -Inf)
})

test_that("log1mexp handles Inf and NAs consistently with current code path", {
  x <- c(-Inf, Inf, NA_real_)
  # Negative inputs include -Inf. We still expect the custom warning;
  # muffle spurious "NaNs produced" warnings.
  expect_warning(
    res <- emit_only_warning(log1mexp(x), "'a' >= 0 needed"),
    "'a' >= 0 needed"
  )

  # -Inf -> NaN (log(1 - exp(-(-Inf))) = log(1 - Inf) = NaN)
  expect_true(is.nan(res[1]))
  # +Inf -> 0 (log(1 - 0) = 0)
  expect_identical(res[2], 0)
  # NA remains NA
  expect_true(is.na(res[3]))
})

test_that("log1mexp respects a custom cutoff but same math", {
  cf <- 0.1
  # Branch selection
  # should use left branch
  x_small <- cf * 0.5
  # should use right branch
  x_large <- cf * 2

  r_small <- log1mexp(x_small, cutoff = cf)
  r_large <- log1mexp(x_large, cutoff = cf)

  expect_equal(r_small, log(-expm1(-x_small)), tolerance = 1e-14)
  expect_equal(r_large, log1p(-exp(-x_large)), tolerance = 1e-14)
})

test_that("log1mexp(a) equals log(1 - exp(-a)) on the stable range", {
  a <- seq(1e-3, 5, length.out = 51)
  expect_equal(log1mexp(a), log1p(-exp(-a)), tolerance = 1e-12)
})

test_that("log1mexp uses the accurate branch for tiny a", {
  # direct log(1 - exp(-a)) loses precision for small a; log1mexp must not.
  a <- 1e-8
  expect_equal(log1mexp(a), log(-expm1(-a)), tolerance = 1e-16)
})

test_that("log1mexp warns on negative a and preserves NA", {
  expect_warning(log1mexp(-1), "'a' >= 0")
  expect_equal(log1mexp(c(1, NA)), c(log1p(-exp(-1)), NA))
})