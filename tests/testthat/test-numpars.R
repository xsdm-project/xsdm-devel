library(testthat)

test_that("num_par() computes the correct sequence and is strictly increasing", {
  # Known small values
  expect_equal(num_par(1), 5L)
  expect_equal(num_par(2), 9L)
  expect_equal(num_par(3), 14L)

  # Monotonicity over a reasonable range
  vals <- num_par(1:50)
  expect_true(all(diff(vals) > 0))
})

test_that("num_env_var() inverts num_par() exactly (round-trip)", {
  for (p in 1:50) {
    n <- num_par(p)
    expect_identical(num_env_var(n), as.integer(p))
  }
})

test_that("num_env_var() rejects invalid n (message fragments from checkmate)", {
  # NA -> checkmate says "Contains missing values"
  expect_error(num_env_var(NA_integer_), "Contains missing values")

  # length > 1 -> checkmate says "Must have length 1"
  expect_error(num_env_var(c(5L, 9L)), "Must have length 1")
})

test_that("num_env_var() rejects values not on the num_par sequence", {
  expect_error(num_env_var(6), "Invalid 'n'")
  expect_error(num_env_var(7), "Invalid 'n'")
  expect_error(num_env_var(8), "Invalid 'n'")
})

test_that("boundary case: p = 1 maps to n = 5 and back", {
  expect_equal(num_par(1), 5L)
  expect_identical(num_env_var(5L), 1L)
})
