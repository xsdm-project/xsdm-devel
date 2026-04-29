library(testthat)
library(xsdm)

# Parity tests between the exported `dist_between_params` (C++ backend) and
# the pure-R reference `xsdm:::dist_between_params_r` (clue::solve_LSAP
# backend). Tolerance 1e-6 matches the other *_r_vs_cpp tests in the suite.

skip_if_not_clue <- function() {
  if (!requireNamespace("clue", quietly = TRUE)) {
    skip("clue not installed")
  }
}

make_random_param_list <- function(p, seed) {
  set.seed(seed)
  raw_o <- matrix(rnorm(p * p), p, p)
  qrd   <- qr(raw_o)
  o_mat <- qr.Q(qrd)
  list(
    mu      = rnorm(p),
    sigltil = exp(rnorm(p, mean = 0.3, sd = 0.3)),
    sigrtil = exp(rnorm(p, mean = 0.3, sd = 0.3)),
    ctil    = rnorm(1),
    pd      = runif(1, 0.1, 0.9),
    o_mat   = o_mat
  )
}

test_that("dist_between_params: distance matches _r on random bio-scale inputs", {
  skip_if_not_clue()
  for (p in 1:3) {
    for (seed in c(11L, 23L, 42L, 101L)) {
      p1 <- make_random_param_list(p, seed = seed)
      p2 <- make_random_param_list(p, seed = seed + 1000L)

      d_cpp <- dist_between_params(p1, p2)
      d_r   <- xsdm:::dist_between_params_r(p1, p2)

      expect_equal(d_cpp, d_r, tolerance = 1e-6,
                   label = sprintf("dist_between_params (p=%d, seed=%d)", p, seed))
    }
  }
})

test_that("dist_between_params: give_closest_rep matches _r (distance + rep)", {
  skip_if_not_clue()
  for (p in 2:3) {
    for (seed in c(7L, 19L, 37L)) {
      p1 <- make_random_param_list(p, seed = seed)
      p2 <- make_random_param_list(p, seed = seed + 500L)

      out_cpp <- dist_between_params(p1, p2, give_closest_rep = TRUE)
      out_r   <- xsdm:::dist_between_params_r(p1, p2, give_closest_rep = TRUE)

      expect_named(out_cpp, c("distance", "representative"))
      expect_equal(out_cpp$distance, out_r$distance, tolerance = 1e-6)

      # Representative fields that are invariant across the equivalence class
      expect_equal(out_cpp$representative$mu,   p1$mu)
      expect_equal(out_cpp$representative$ctil, p1$ctil)
      expect_equal(out_cpp$representative$pd,   p1$pd)

      # The representative produced by C++ must reproduce the distance when
      # scored against p2 through the brute-force reference.
      d_via_rep <- xsdm:::distance_between_params_r(out_cpp$representative, p2)
      expect_equal(d_via_rep, out_cpp$distance, tolerance = 1e-6)
    }
  }
})

test_that("dist_between_params: math-scale inputs match _r", {
  skip_if_not_clue()
  v1 <- examples$optim_par_vec
  v2 <- examples$optim_par_vec_equivalent

  expect_equal(dist_between_params(v1, v2),
               xsdm:::dist_between_params_r(v1, v2),
               tolerance = 1e-6)

  # Should be ~0 since examples expose an equivalent representative.
  expect_lt(dist_between_params(v1, v2), 1e-6)
})

test_that(".solve_lsap_cpp matches clue::solve_LSAP on random cost matrices", {
  skip_if_not_clue()
  for (n in c(1L, 2L, 5L, 10L)) {
    for (seed in c(1L, 2L, 3L, 4L)) {
      set.seed(seed + n * 1000L)
      cost <- matrix(runif(n * n, min = 0, max = 10), n, n)
      # Keep integer to reduce tie ambiguity for larger matrices.
      cost <- round(cost, 6)

      p_cpp  <- as.integer(xsdm:::.solve_lsap_cpp(cost))
      p_clue <- as.integer(clue::solve_LSAP(cost))

      # Both must produce optimal assignments; the *perm* can differ when the
      # optimum is attained by multiple pairings, so compare total cost.
      cost_cpp  <- sum(cost[cbind(seq_len(n), p_cpp)])
      cost_clue <- sum(cost[cbind(seq_len(n), p_clue)])
      expect_equal(cost_cpp, cost_clue, tolerance = 1e-10,
                   label = sprintf("solve_lsap_cpp (n=%d, seed=%d)", n, seed))
    }
  }
})

test_that("dist_between_params: canonical == unexported _r for examples fixture", {
  skip_if_not_clue()
  par_list            <- math_to_bio(examples$optim_par_vec)
  par_list_equivalent <- math_to_bio(examples$optim_par_vec_equivalent)

  expect_equal(
    dist_between_params(par_list, par_list_equivalent),
    xsdm:::dist_between_params_r(par_list, par_list_equivalent),
    tolerance = 1e-6
  )
})
