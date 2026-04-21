# tests/testthat/test-permutations.R
library(testthat)

test_that("permutations() - repeats.allowed = TRUE hits all inner branches", {
  # Case 1: r == 1 (first branch inside sub() when repeats are allowed)
  res_r1 <- xsdm:::permutations(
    n = 3, r = 1, v = letters[1:3],
    set = TRUE, repeats.allowed = TRUE
  )
  expect_true(is.matrix(res_r1))
  expect_equal(dim(res_r1), c(3, 1))
  expect_equal(as.vector(res_r1), letters[1:3])
  
  # Case 2: n == 1 (second branch inside sub() when repeats are allowed)
  # Here n = 1 => v has length 1; with repeats allowed we can request r > 1
  res_n1 <- xsdm:::permutations(
    n = 1, r = 3, v = "x",
    set = TRUE, repeats.allowed = TRUE
  )
  expect_true(is.matrix(res_n1))
  expect_equal(dim(res_n1), c(1, 3))
  expect_equal(as.vector(res_n1), rep("x", 3))
  
  # Case 3: general recursive branch (neither r == 1 nor n == 1)
  # n = 2, r = 2, repeats allowed: all 2^2 = 4 ordered pairs with replacement
  res_gen <- xsdm:::permutations(
    n = 2, r = 2, v = c(10, 20),
    set = TRUE, repeats.allowed = TRUE
  )
  expect_true(is.matrix(res_gen))
  expect_equal(dim(res_gen), c(4, 2))
  
  # Robust check of the set of rows, ignoring order:
  # (10,10), (10,20), (20,10), (20,20)
  exp_df  <- expand.grid(c(10, 20), c(10, 20))  # columns Var1, Var2
  got_key <- apply(res_gen, 1, paste, collapse = ",")
  exp_key <- apply(as.matrix(exp_df[, 2:1]), 1, paste, collapse = ",") # same col order as res_gen
  expect_setequal(got_key, exp_key)
})

test_that("permutations() - regular path without repeats", {
  # n = 3, r = 2, no repeats => 3 * 2 = 6 rows
  res <- xsdm:::permutations(
    n = 3, r = 2, v = c("a", "b", "c"),
    set = TRUE, repeats.allowed = FALSE
  )
  expect_true(is.matrix(res))
  expect_equal(dim(res), c(6, 2))
  
  # All permutations of size 2 without repetition from {a, b, c}
  expected <- rbind(
    c("a","b"), c("a","c"),
    c("b","a"), c("b","c"),
    c("c","a"), c("c","b")
  )
  expect_true(all(
    apply(expected, 1, paste, collapse = ",") %in%
      apply(res, 1, paste, collapse = ",")
  ))
})

test_that("permutations() - errors: bad value of n", {
  # n must be a single positive integer (numeric mode)
  expect_error(xsdm:::permutations(n = 0,   r = 1), "bad value of n")
  expect_error(xsdm:::permutations(n = -2,  r = 1), "bad value of n")
  expect_error(xsdm:::permutations(n = 1.5, r = 1), "bad value of n")
  expect_error(xsdm:::permutations(n = c(2,3), r = 1), "bad value of n")
  expect_error(xsdm:::permutations(n = "2", r = 1), "bad value of n")
})

test_that("permutations() - errors: bad value of r", {
  # r must be a single positive integer (numeric mode)
  expect_error(xsdm:::permutations(n = 2, r = 0),        "bad value of r")
  expect_error(xsdm:::permutations(n = 2, r = -1),       "bad value of r")
  expect_error(xsdm:::permutations(n = 2, r = 1.5),      "bad value of r")
  expect_error(xsdm:::permutations(n = 2, r = c(1,2)),   "bad value of r")
  expect_error(xsdm:::permutations(n = 2, r = "1"),      "bad value of r")
})

test_that("permutations() - errors: v non-atomic or too short", {
  # v must be atomic and length >= n
  expect_error(
    xsdm:::permutations(n = 2, r = 1, v = list(1, 2)),
    "v is either non-atomic or too short"
  )
  expect_error(
    xsdm:::permutations(n = 3, r = 1, v = c(1, 2)),
    "v is either non-atomic or too short"
  )
})

test_that("permutations() - errors: r > n and repeats.allowed = FALSE", {
  expect_error(
    xsdm:::permutations(
      n = 2, r = 3, v = 1:3, set = TRUE, repeats.allowed = FALSE
    ),
    "r > n and repeats.allowed=FALSE"
  )
})

test_that("permutations() - errors: set=TRUE with too few distinct elements", {
  # unique(sort(v)) becomes shorter than n, which should error
  expect_error(
    xsdm:::permutations(
      n = 3, r = 2, v = c(5, 5, 5), set = TRUE, repeats.allowed = TRUE
    ),
    "too few different elements"
  )
})

test_that("permutations() - r > n allowed when repeats are allowed", {
  # With repeats.allowed = TRUE, this should work and produce n^r rows
  res <- xsdm:::permutations(
    n = 2, r = 3, v = c("x","y"),
    set = TRUE, repeats.allowed = TRUE
  )
  expect_true(is.matrix(res))
  expect_equal(dim(res), c(8, 3)) # 2^3 = 8
  # Minimal check: all values are from the provided alphabet
  expect_true(all(apply(res, 2, function(col) all(col %in% c("x","y")))))
})

