# tests/testthat/test-loglik_math_cpp.R
#
# Parity: xsdm:::loglik_math_cpp() must agree with loglik_math() — including the
# `negative` flag and mask handling.

test_that("loglik_math_cpp: matches loglik_math on examples fixture", {
  pv <- examples$par_vec
  env_dat <- examples$env_array
  occ <- as.integer(examples$occ_vec)

  for (neg in c(TRUE, FALSE)) {
    ref <- loglik_math(pv, env_dat, occ, negative = neg)
    got <- xsdm:::loglik_math_cpp(pv, env_dat, occ, negative = neg)
    expect_equal(got, ref, tolerance = 1e-10,
                 info = paste0("negative = ", neg))
  }
})

test_that("loglik_math_cpp: matches loglik_math on a small p = 1 fixture", {
  set.seed(5)
  n <- 8
  Tt <- 4
  env_dat <- array(runif(n * Tt * 1, -2, 2), dim = c(n, Tt, 1))
  occ <- rep(c(1L, 0L), length.out = n)

  pv <- make_mask_names(1)
  pv[] <- c(mu1 = 0.1, sigltil1 = log(1.1), sigrtil1 = log(0.9),
            ctil = -0.5, pd = 0.3)

  ref <- loglik_math(pv, env_dat, occ, negative = TRUE)
  got <- xsdm:::loglik_math_cpp(pv, env_dat, occ, negative = TRUE)
  expect_equal(got, ref, tolerance = 1e-10)
})

test_that("loglik_math_cpp: matches loglik_math with a mask", {
  set.seed(9)
  n <- 10
  Tt <- 5
  p <- 2
  env_dat <- array(runif(n * Tt * p, -1, 1), dim = c(n, Tt, p))
  occ <- rep(c(1L, 0L), length.out = n)

  # free params; mask fixes the rest
  free <- c(mu1 = 0.2, mu2 = -0.3, sigltil1 = 0.1, sigltil2 = -0.1,
            sigrtil1 = 0.2, sigrtil2 = -0.2, ctil = -0.4, o_par1 = 0.5)
  mask <- c(pd = -1.0)

  ref <- loglik_math(free, env_dat, occ, mask = mask, negative = TRUE)
  got <- xsdm:::loglik_math_cpp(free, env_dat, occ, mask = mask, negative = TRUE)
  expect_equal(got, ref, tolerance = 1e-10)
})

test_that("loglik_math_cpp: random p = 2 parity", {
  set.seed(2028)
  for (trial in 1:5) {
    n <- 10
    Tt <- 6
    p <- 2
    env_dat <- array(runif(n * Tt * p, -1, 1), dim = c(n, Tt, p))
    occ <- rbinom(n, 1, 0.5)

    pv <- make_mask_names(p)
    pv[] <- rnorm(length(pv), sd = 0.5)

    ref <- loglik_math(pv, env_dat, occ, negative = TRUE)
    got <- xsdm:::loglik_math_cpp(pv, env_dat, occ, negative = TRUE)
    expect_equal(got, ref, tolerance = 1e-10)
  }
})

test_that("loglik_math_cpp: num_threads determinism", {
  pv <- examples$par_vec
  env_dat <- examples$env_array
  occ <- as.integer(examples$occ_vec)
  r1 <- xsdm:::loglik_math_cpp(pv, env_dat, occ, num_threads = 1L)
  r2 <- xsdm:::loglik_math_cpp(pv, env_dat, occ, num_threads = 2L)
  expect_equal(r1, r2, tolerance = 1e-12)
})
