library(testthat)

# ---------------------------------------------------------------------------
# Tests for the canonical wrapper `log_prob_detect()`.
#
# These tests focus on contract-level invariants that hold regardless of the
# specific value returned by the like_neg_ltsgr / log_prob_detect_cpp kernel:
#
#   * shape contract  (length of returned vector == n_loc)
#   * pd-boundary behaviour    (pd = 0 -> log = -Inf; prob = 0)
#   * return_prob duality      (prob == exp(log_prob), elementwise)
#   * pd-scaling identity      (changing only pd shifts log_prob by log ratio)
#
# Numerical equality of the wrapper to a pure-R reference implementation is
# covered separately by tests/testthat/test-log_prob_detect_r_vs_cpp.R, which
# sweeps p in {1,2,3} at tolerance 1e-6.
#
# Note: prior versions of this file used `with_mocked_bindings()` to override
# `like_neg_ltsgr_cpp` to a constant return value, then asserted the
# implementation-aware formula log(pd) - log1pexp(ctil + 2.0). After
# `log_prob_detect()` was refactored to delegate directly to
# `log_prob_detect_cpp`, that internal call is no longer reachable through
# the public wrapper and the mock has no effect; the tests have therefore
# been rewritten to assert invariants instead of implementation details.
# ---------------------------------------------------------------------------

# A canonical, well-shaped fixture used by all three tests below: 5 locations,
# 5 time steps, p = 2 environmental variables. Matches mu = c(11, 5),
# 2x2 o_mat, and length-2 sigltil/sigrtil.
make_fixture <- function() {
  n_loc <- 5L
  n_time <- 5L
  p <- 2L
  set.seed(1L)
  list(
    env_dat = array(stats::runif(n_loc * n_time * p, 5, 20),
                    dim = c(n_loc, n_time, p)),
    mu      = c(11, 5),
    sigltil = c(1, 2),
    sigrtil = c(2, 1),
    ctil    = -2,
    o_mat   = matrix(c(-0.4, 0.9, -0.9, -0.4), ncol = 2)
  )
}

test_that("log_prob_detect: returns a length-n_loc numeric vector", {
  f <- make_fixture()
  out <- log_prob_detect(f$env_dat, f$mu, f$sigltil, f$sigrtil,
                         f$o_mat, f$ctil, pd = 0.9)
  expect_type(out, "double")
  expect_length(out, dim(f$env_dat)[1])
  expect_false(anyNA(out))
})

test_that("log_prob_detect: return_prob = TRUE matches exp(log_prob)", {
  f <- make_fixture()
  log_p <- log_prob_detect(f$env_dat, f$mu, f$sigltil, f$sigrtil,
                           f$o_mat, f$ctil, pd = 0.9, return_prob = FALSE)
  prob  <- log_prob_detect(f$env_dat, f$mu, f$sigltil, f$sigrtil,
                           f$o_mat, f$ctil, pd = 0.9, return_prob = TRUE)
  expect_equal(prob, exp(log_p))
})

test_that("log_prob_detect: pd = 0 yields -Inf log-prob and 0 probability", {
  f <- make_fixture()
  log_p <- log_prob_detect(f$env_dat, f$mu, f$sigltil, f$sigrtil,
                           f$o_mat, f$ctil, pd = 0, return_prob = FALSE)
  prob  <- log_prob_detect(f$env_dat, f$mu, f$sigltil, f$sigrtil,
                           f$o_mat, f$ctil, pd = 0, return_prob = TRUE)
  expect_true(all(is.infinite(log_p) & log_p < 0))
  expect_equal(prob, rep(0, length(prob)))
})

test_that("log_prob_detect: changing only pd shifts log-prob by log(ratio)", {
  # The wrapper computes log(pd) - log1pexp(ctil + kernel(...)). The kernel
  # term does not depend on `pd`, so multiplying pd by a constant must shift
  # every entry of log_prob by exactly log(constant), independent of the
  # underlying kernel value. This invariant survives any future change to
  # the kernel implementation.
  f <- make_fixture()
  a <- log_prob_detect(f$env_dat, f$mu, f$sigltil, f$sigrtil,
                       f$o_mat, f$ctil, pd = 1.0)
  b <- log_prob_detect(f$env_dat, f$mu, f$sigltil, f$sigrtil,
                       f$o_mat, f$ctil, pd = 0.5)
  expect_equal(b - a, rep(log(0.5), length(a)))
})

test_that("log_prob_detect: rejects mismatched mu length (C++ contract)", {
  f <- make_fixture()
  bad_mu <- c(f$mu, 0)  # length p+1 instead of p
  expect_error(
    log_prob_detect(f$env_dat, bad_mu, f$sigltil, f$sigrtil,
                    f$o_mat, f$ctil, pd = 0.9),
    "mu"
  )
})
