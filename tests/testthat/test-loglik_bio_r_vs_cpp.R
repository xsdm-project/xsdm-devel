# tests/testthat/test-loglik_bio_r_vs_cpp.R
#
# Parity: the canonical exported `loglik_bio()` (thin C++ wrapper) must
# agree numerically with the pure-R reference `loglik_bio_r()` kept in
# R/internals.R. Tolerance is 1e-6 (ecological precision).

TOL <- 1e-6

# Helper to draw a random biological-scale parameter set.
draw_bio <- function(p, seed) {
  set.seed(seed)
  list(
    mu      = runif(p, -1, 1),
    sigltil = exp(runif(p, -0.5, 0.5)),       # positive
    sigrtil = exp(runif(p, -0.5, 0.5)),       # positive
    ctil    = runif(1, -1, 1),
    pd      = runif(1, 0.1, 0.9),             # in (0, 1)
    o_mat   = build_orthogonal_matrix(if (p == 1) NULL else runif(p * (p - 1) / 2, -1, 1))
  )
}

test_that("loglik_bio vs loglik_bio_r: example_1 fixture, default flags", {
  ll_cpp <- loglik_bio(
    env_dat = example_1$env_array,
    occ     = example_1$occ_df$presence,
    mu      = example_1$par_list$mu,
    sigltil = example_1$par_list$sigltil,
    sigrtil = example_1$par_list$sigrtil,
    o_mat   = example_1$par_list$o_mat,
    ctil    = example_1$par_list$ctil,
    pd      = example_1$par_list$pd
  )
  ll_r <- xsdm:::loglik_bio_r(
    env_dat = example_1$env_array,
    occ     = example_1$occ_df$presence,
    mu      = example_1$par_list$mu,
    sigltil = example_1$par_list$sigltil,
    sigrtil = example_1$par_list$sigrtil,
    o_mat   = example_1$par_list$o_mat,
    ctil    = example_1$par_list$ctil,
    pd      = example_1$par_list$pd
  )
  expect_equal(ll_cpp, ll_r, tolerance = TOL)
})

test_that("loglik_bio vs loglik_bio_r: random p in {1, 2, 3}, default flags", {
  seed_base <- 1100
  for (p in 1:3) {
    n  <- 10
    Tt <- 5
    env_dat <- array(runif(n * Tt * p, -1, 1), dim = c(n, Tt, p))
    occ <- rbinom(n, 1, 0.5)
    bp  <- draw_bio(p, seed = seed_base + p)

    ll_cpp <- loglik_bio(env_dat, occ, bp$mu, bp$sigltil, bp$sigrtil,
                         bp$o_mat, bp$ctil, bp$pd)
    ll_r   <- xsdm:::loglik_bio_r(env_dat, occ, bp$mu, bp$sigltil, bp$sigrtil,
                                  bp$o_mat, bp$ctil, bp$pd)
    expect_equal(ll_cpp, ll_r, tolerance = TOL,
                 info = paste0("p = ", p))
  }
})

test_that("loglik_bio vs loglik_bio_r: sum_log_p = FALSE returns matching per-location vector", {
  set.seed(7)
  p  <- 2
  n  <- 12
  Tt <- 6
  env_dat <- array(runif(n * Tt * p, -1, 1), dim = c(n, Tt, p))
  occ <- rbinom(n, 1, 0.5)
  bp  <- draw_bio(p, seed = 1207)

  v_cpp <- loglik_bio(env_dat, occ, bp$mu, bp$sigltil, bp$sigrtil,
                      bp$o_mat, bp$ctil, bp$pd, sum_log_p = FALSE)
  v_r   <- xsdm:::loglik_bio_r(env_dat, occ, bp$mu, bp$sigltil, bp$sigrtil,
                               bp$o_mat, bp$ctil, bp$pd, sum_log_p = FALSE)
  expect_length(v_cpp, n)
  expect_equal(v_cpp, v_r, tolerance = TOL)
})

test_that("loglik_bio vs loglik_bio_r: return_prob = TRUE gives matching linear-scale value", {
  set.seed(11)
  p  <- 2
  n  <- 8
  Tt <- 4
  env_dat <- array(runif(n * Tt * p, -1, 1), dim = c(n, Tt, p))
  occ <- rbinom(n, 1, 0.5)
  bp  <- draw_bio(p, seed = 1108)

  v_cpp <- loglik_bio(env_dat, occ, bp$mu, bp$sigltil, bp$sigrtil,
                      bp$o_mat, bp$ctil, bp$pd, return_prob = TRUE)
  v_r   <- xsdm:::loglik_bio_r(env_dat, occ, bp$mu, bp$sigltil, bp$sigrtil,
                               bp$o_mat, bp$ctil, bp$pd, return_prob = TRUE)
  expect_equal(v_cpp, v_r, tolerance = TOL)
})
