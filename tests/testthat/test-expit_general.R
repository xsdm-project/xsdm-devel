library(testthat)

test_that("expit_general returns value between bounds", {
  result <- expit_general(0, -1, 1, 0)
  expect_true(result >= -1 && result <= 1)
})

test_that("expit_general returns midpoint at x0", {
  result <- expit_general(0, -1, 1, 0)
  expect_equal(round(result, 6), 0) # midpoint should be 0
})

test_that("expit_general handles vector input", {
  x <- c(-10, 0, 10)
  result <- expit_general(x, -1, 1, 0)
  expect_equal(length(result), 3)
  expect_true(all(result >= -1 & result <= 1))
})

test_that("expit_general caps at upper bound for large x", {
  expect_equal(expit_general(200, -1, 1, 0), 1)
})

test_that("expit_general caps at lower bound for very negative x", {
  expect_equal(expit_general(-200, -1, 1, 0), -1)
})

test_that("expit_general fails on non-numeric input", {
  expect_error(expit_general("a", -1, 1, 0), "numeric")
})

test_that("expit_general fails when bounds are not numeric", {
  expect_error(expit_general(0, "low", 1, 0), "numeric")
})
