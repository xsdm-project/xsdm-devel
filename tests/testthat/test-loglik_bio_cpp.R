# tests/testthat/test-loglik_bio_cpp.R
#
# Parity: xsdm:::loglik_bio_cpp() must match loglik_bio(..., sum_log_p = TRUE,
# return_prob = FALSE) to near machine precision.

env_flat <- function(env_dat) {
  list(
    vec  = as.numeric(env_dat),
    dims = as.integer(dim(env_dat))
  )
}

test_that("loglik_bio_cpp: matches loglik_bio on example_1 fixture", {
  env <- env_flat(example_1$env_array)
  ref <- loglik_bio(
    env_dat = example_1$env_array,
    occ     = example_1$occ_vec,
    mu      = example_1$par_list$mu,
    sigltil = example_1$par_list$sigltil,
    sigrtil = example_1$par_list$sigrtil,
    o_mat   = example_1$par_list$o_mat,
    ctil    = example_1$par_list$ctil,
    pd      = example_1$par_list$pd
  )

  got <- xsdm:::loglik_bio_cpp(
    env_dat_vec  = env$vec,
    env_dat_dims = env$dims,
    occ          = as.integer(example_1$occ_vec),
    mu           = example_1$par_list$mu,
    sigltil      = example_1$par_list$sigltil,
    sigrtil      = example_1$par_list$sigrtil,
    o_mat        = example_1$par_list$o_mat,
    ctil         = example_1$par_list$ctil,
    pd           = example_1$par_list$pd
  )

  expect_equal(got, ref, tolerance = 1e-10)
})

test_that("loglik_bio_cpp: matches loglik_bio on a small p = 1 fixture", {
  set.seed(7)
  n <- 6
  Tt <- 4
  env_dat <- array(runif(n * Tt * 1, -2, 2), dim = c(n, Tt, 1))
  occ <- rep(c(1L, 0L), length.out = n)

  mu <- 0.3
  sigltil <- 1.2
  sigrtil <- 0.9
  o_mat <- matrix(1, 1, 1)
  ctil <- -0.4
  pd <- 0.7

  ref <- loglik_bio(env_dat, occ, mu, sigltil, sigrtil, o_mat, ctil, pd)
  env <- env_flat(env_dat)
  got <- xsdm:::loglik_bio_cpp(env$vec, env$dims, occ,
                       mu, sigltil, sigrtil, o_mat, ctil, pd)
  expect_equal(got, ref, tolerance = 1e-10)
})

test_that("loglik_bio_cpp: matches loglik_bio on random p = 2 draws", {
  set.seed(2027)
  for (trial in 1:5) {
    n <- 10
    Tt <- 5
    p <- 2
    env_dat <- array(runif(n * Tt * p, -1, 1), dim = c(n, Tt, p))
    occ <- rbinom(n, 1, 0.4)

    pv <- make_mask_names(p)
    pv[] <- rnorm(length(pv), sd = 0.6)
    bp <- math_to_bio(pv)

    ref <- loglik_bio(env_dat, occ,
                      bp$mu, bp$sigltil, bp$sigrtil,
                      bp$o_mat, bp$ctil, bp$pd)

    env <- env_flat(env_dat)
    got <- xsdm:::loglik_bio_cpp(env$vec, env$dims, as.integer(occ),
                          bp$mu, bp$sigltil, bp$sigrtil,
                          bp$o_mat, bp$ctil, bp$pd)
    expect_equal(got, ref, tolerance = 1e-10)
  }
})

test_that("loglik_bio_cpp: honours num_threads and matches 1-thread result", {
  set.seed(11)
  n <- 12
  Tt <- 5
  p <- 2
  env_dat <- array(runif(n * Tt * p, -1, 1), dim = c(n, Tt, p))
  occ <- rbinom(n, 1, 0.5)
  pv  <- make_mask_names(p)
  pv[] <- rnorm(length(pv))
  bp <- math_to_bio(pv)
  env <- env_flat(env_dat)

  r1 <- xsdm:::loglik_bio_cpp(env$vec, env$dims, as.integer(occ),
                       bp$mu, bp$sigltil, bp$sigrtil, bp$o_mat,
                       bp$ctil, bp$pd, num_threads = 1L)
  r2 <- xsdm:::loglik_bio_cpp(env$vec, env$dims, as.integer(occ),
                       bp$mu, bp$sigltil, bp$sigrtil, bp$o_mat,
                       bp$ctil, bp$pd, num_threads = 2L)
  expect_equal(r1, r2, tolerance = 1e-12)
})
