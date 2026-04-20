library(testthat)

test_that("expit returns correct values for numeric input", {
  expect_equal(expit(0), 0.5)
  expect_equal(round(expit(1), 6), round(exp(1) / (1 + exp(1)), 6))
})

test_that("expit handles vector input correctly", {
  x <- c(-1, 0, 1)
  result <- expit(x)
  expect_equal(length(result), 3)
  expect_true(all(result > 0 & result < 1))
})

test_that("expit caps values at 1 for large positive input", {
  expect_equal(expit(1000), 1)
  expect_equal(expit(200), 1)
})

test_that("expit approaches 0 for large negative input", {
  expect_true(expit(-1000) < 1e-6)
})

test_that("expit fails on non-numeric input", {
  expect_error(expit("a"), "numeric")
})


expect_na <- function(object) {
  expect_true(is.na(object), info = "Expected value to be NA")
}

test_that("expit returns NA values", {
  expect_na(expit(NA))
})

