# tests/testthat/test-loglik_math_r_vs_cpp.R
#
# Parity: the canonical exported `loglik_math()` (thin C++ wrapper) must
# agree numerically with the pure-R reference `loglik_math_r()` kept in
# R/internals.R. Tolerance is 1e-6, which is well inside ecological
# precision while still catching any algorithmic divergence between the
# C++ kernel and the historical R implementation.

TOL <- 1e-6

test_that("loglik_math vs loglik_math_r: examples fixture, both `negative` flags", {
  pv      <- examples$par_vec
  env_dat <- examples$env_array
  occ     <- as.integer(examples$occ_vec)

  for (neg in c(TRUE, FALSE)) {
    cpp_val <- loglik_math(pv, env_dat, occ, negative = neg)
    r_val   <- xsdm:::loglik_math_r(pv, env_dat, occ, negative = neg)
    expect_equal(cpp_val, r_val, tolerance = TOL,
                 info = paste0("negative = ", neg))
  }
})

test_that("loglik_math vs loglik_math_r: p = 1 random sweep", {
  set.seed(101)
  for (trial in 1:5) {
    n  <- sample(5:12, 1)
    Tt <- sample(3:8, 1)
    env_dat <- array(runif(n * Tt * 1, -2, 2), dim = c(n, Tt, 1))
    occ <- rbinom(n, 1, 0.5)
    pv <- c(mu1 = runif(1, -1, 1),
            sigltil1 = runif(1, -0.5, 0.5),
            sigrtil1 = runif(1, -0.5, 0.5),
            ctil = runif(1, -1, 1),
            pd = runif(1, -1, 1))
    cpp_val <- loglik_math(pv, env_dat, occ, negative = TRUE)
    r_val   <- xsdm:::loglik_math_r(pv, env_dat, occ, negative = TRUE)
    expect_equal(cpp_val, r_val, tolerance = TOL,
                 info = paste0("p = 1, trial ", trial))
  }
})

test_that("loglik_math vs loglik_math_r: p = 2 random sweep", {
  set.seed(202)
  for (trial in 1:5) {
    n  <- sample(8:15, 1)
    Tt <- sample(3:8, 1)
    p  <- 2
    env_dat <- array(runif(n * Tt * p, -1, 1), dim = c(n, Tt, p))
    occ <- rbinom(n, 1, 0.5)
    pv <- make_mask_names(p)
    pv[] <- c(runif(p, -1, 1),               # mu
              runif(p, -0.5, 0.5),           # sigltil
              runif(p, -0.5, 0.5),           # sigrtil
              runif(1, -1, 1),               # ctil
              runif(1, -1, 1),               # pd
              runif(p * (p - 1) / 2, -1, 1)) # o_par
    cpp_val <- loglik_math(pv, env_dat, occ, negative = FALSE)
    r_val   <- xsdm:::loglik_math_r(pv, env_dat, occ, negative = FALSE)
    expect_equal(cpp_val, r_val, tolerance = TOL,
                 info = paste0("p = 2, trial ", trial))
  }
})

test_that("loglik_math vs loglik_math_r: p = 3 with mask", {
  set.seed(303)
  n  <- 12
  Tt <- 6
  p  <- 3
  env_dat <- array(runif(n * Tt * p, -1, 1), dim = c(n, Tt, p))
  occ <- rbinom(n, 1, 0.4)

  # Free params: leave out one mu and pd -> mask fixes them.
  full <- make_mask_names(p)
  full[] <- c(runif(p, -1, 1),
              runif(p, -0.5, 0.5),
              runif(p, -0.5, 0.5),
              runif(1, -1, 1),
              runif(1, -1, 1),
              runif(p * (p - 1) / 2, -1, 1))
  mask <- full[c("mu2", "pd")]
  free <- full[setdiff(names(full), names(mask))]

  cpp_val <- loglik_math(free, env_dat, occ, mask = mask, negative = TRUE)
  r_val   <- xsdm:::loglik_math_r(free, env_dat, occ, mask = mask, negative = TRUE)
  expect_equal(cpp_val, r_val, tolerance = TOL)
})

test_that("loglik_math vs loglik_math_r: all-presence and all-absence", {
  set.seed(404)
  n  <- 10
  Tt <- 5
  p  <- 2
  env_dat <- array(runif(n * Tt * p, -1, 1), dim = c(n, Tt, p))
  pv <- make_mask_names(p)
  pv[] <- c(runif(p, -1, 1), runif(p, -0.5, 0.5), runif(p, -0.5, 0.5),
            runif(1, -1, 1), runif(1, -1, 1), runif(p * (p - 1) / 2, -1, 1))

  for (occ in list(rep(1L, n), rep(0L, n))) {
    cpp_val <- loglik_math(pv, env_dat, occ, negative = TRUE)
    r_val   <- xsdm:::loglik_math_r(pv, env_dat, occ, negative = TRUE)
    expect_equal(cpp_val, r_val, tolerance = TOL,
                 info = paste0("occ = ", paste(occ, collapse = "")))
  }
})
