# tests/testthat/test-log_prob_detect_r_vs_cpp.R
#
# Parity: the canonical exported `log_prob_detect()` (thin C++ wrapper) must
# agree numerically with the pure-R reference `log_prob_detect_r()` kept in
# R/internals.R. Tolerance is 1e-6 (ecological precision).

TOL <- 1e-6

draw_bio <- function(p, seed) {
  set.seed(seed)
  list(
    mu      = runif(p, -1, 1),
    sigltil = exp(runif(p, -0.5, 0.5)),
    sigrtil = exp(runif(p, -0.5, 0.5)),
    ctil    = runif(1, -1, 1),
    pd      = runif(1, 0.1, 0.9),
    o_mat   = build_orthogonal_matrix(if (p == 1) NULL else runif(p * (p - 1) / 2, -1, 1))
  )
}

test_that("log_prob_detect vs log_prob_detect_r: examples fixture", {
  v_cpp <- log_prob_detect(
    env_dat = examples$env_array,
    mu      = examples$par_list$mu,
    sigltil = examples$par_list$sigltil,
    sigrtil = examples$par_list$sigrtil,
    o_mat   = examples$par_list$o_mat,
    ctil    = examples$par_list$ctil,
    pd      = examples$par_list$pd
  )
  v_r <- xsdm:::log_prob_detect_r(
    env_dat = examples$env_array,
    mu      = examples$par_list$mu,
    sigltil = examples$par_list$sigltil,
    sigrtil = examples$par_list$sigrtil,
    o_mat   = examples$par_list$o_mat,
    ctil    = examples$par_list$ctil,
    pd      = examples$par_list$pd
  )
  expect_equal(v_cpp, v_r, tolerance = TOL)
})

test_that("log_prob_detect vs log_prob_detect_r: random p in {1, 2, 3}", {
  for (p in 1:3) {
    set.seed(2200 + p)
    n  <- 10
    Tt <- 5
    env_dat <- array(runif(n * Tt * p, -1, 1), dim = c(n, Tt, p))
    bp <- draw_bio(p, seed = 2200 + p)

    v_cpp <- log_prob_detect(env_dat, bp$mu, bp$sigltil, bp$sigrtil,
                             bp$o_mat, bp$ctil, bp$pd)
    v_r   <- xsdm:::log_prob_detect_r(env_dat, bp$mu, bp$sigltil, bp$sigrtil,
                                      bp$o_mat, bp$ctil, bp$pd)
    expect_length(v_cpp, n)
    expect_equal(v_cpp, v_r, tolerance = TOL,
                 info = paste0("p = ", p))
  }
})

test_that("log_prob_detect vs log_prob_detect_r: return_prob = TRUE", {
  set.seed(99)
  p  <- 2
  n  <- 12
  Tt <- 6
  env_dat <- array(runif(n * Tt * p, -1, 1), dim = c(n, Tt, p))
  bp <- draw_bio(p, seed = 9921)

  v_cpp <- log_prob_detect(env_dat, bp$mu, bp$sigltil, bp$sigrtil,
                           bp$o_mat, bp$ctil, bp$pd, return_prob = TRUE)
  v_r   <- xsdm:::log_prob_detect_r(env_dat, bp$mu, bp$sigltil, bp$sigrtil,
                                    bp$o_mat, bp$ctil, bp$pd, return_prob = TRUE)
  expect_equal(v_cpp, v_r, tolerance = TOL)
})
