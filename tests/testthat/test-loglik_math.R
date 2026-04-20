library(testthat)

test_that("loglik_math: positive path mirrors loglik_bio (sign + threads)", {
  # --- minimal fixtures (as in test-loglik_bio) ---
  n_loc <- 2L
  t_len <- 2L
  n_env <- 2L
  env_dat <- array(0, dim = c(n_loc, t_len, n_env))
  occ <- c(1, 0)

  mu <- c(1, 1)
  sigltil <- c(1, 1)
  sigrtil <- c(1, 1)
  ctil <- 0.1
  pd <- 0.5
  o_mat <- diag(n_env)

  # --- we need the math <-> bio converter(s) present in the package ---
  # math_to_bio is required to build a valid param_vector"

  # To do: bio_to_math Build a math-space parameter vector from a known-good
  # biological list
  # listofpars To bio_to_math To param_vector

  param_vector <- c(
    mu1 = 1,
    mu2 = 1,
    sigltil1 = 1,
    sigltil2 = 1,
    sigrtil1 = 1,
    sigrtil2 = 1,
    ctil = 1,
    pd = 1,
    o_par1 = 1
  )


  # By default negative = TRUE, so loglik_math should return a positive value
  # of log-lik

  expect_silent(
    ll_math_neg <- loglik_math(param_vector,
      env_dat,
      occ,
      negative = TRUE
    )
  )
  expect_type(ll_math_neg, "double")
  expect_length(ll_math_neg, 1L)
  expect_true(ll_math_neg > 0)
  #
  # # Exercise num_threads branch (without asserting on global state)
  expect_silent(loglik_math(param_vector, env_dat, occ, num_threads = 1L))
})

test_that("loglik_math: opt vector fixes parameters (no-ops still run)", {
  n_loc <- 2L
  t_len <- 2L
  n_env <- 2L
  env_dat <- array(0, dim = c(n_loc, t_len, n_env))
  occ <- c(1, 0)

  mu <- c(1, 1)
  sigltil <- c(1, 1)
  sigrtil <- c(1, 1)
  ctil <- 0.1
  pd <- 0.5
  o_mat <- diag(n_env)


  param_vector <- c(
    mu1 = 13.1,
    mu2 = 5.4,
    sigltil1 = 0.9,
    sigltil2 = -0.4,
    sigrtil1 = 0.3,
    sigrtil2 = -0.5,
    ctil = -4.4,
    o_par1 = -9.5
  )


  # Provide an opt vector that "fixes" everything to the same values (no-op),
  # just to cover the replacement branch safely.
  mask <- c(pd = -1.7)

  expect_silent(
    res <- loglik_math(param_vector,
      env_dat,
      occ,
      mask = mask,
      negative = FALSE
    )
  )
  # Should remain a valid numeric scalar
  expect_true(is.numeric(res) && length(res) == 1L && !is.na(res))
})

test_that("loglik_math: input validation errors on bad env_dat and occ", {
  n_loc <- 2L
  t_len <- 2L
  n_env <- 2L
  env_dat <- array(0, dim = c(n_loc, t_len, n_env))
  occ <- c(1, 0)

  mu <- c(1, 1)
  sigltil <- c(1, 1)
  sigrtil <- c(1, 1)
  ctil <- 0.1
  pd <- 0.5
  o_mat <- diag(n_env)


  param_vector <- c(
    mu1 = 13.1,
    mu2 = 5.4,
    sigltil1 = 0.9,
    sigltil2 = -0.4,
    sigrtil1 = 0.3,
    sigrtil2 = -0.5,
    ctil = -4.4,
    pd = -1.7,
    o_par1 = -9.5
  )

  # Bad env_dat
  expect_error(loglik_math(param_vector, 1:5, occ), regexp = "env_dat")

  # Bad occ
  expect_error(loglik_math(param_vector, env_dat, c(-1, 0)), regexp = "occ")
})

test_that(
  "loglik_math fails before conversion if bio params are accessed",
  {
    n_loc <- 2L
    t_len <- 2L
    n_env <- 2L
    env_dat <- array(0, dim = c(n_loc, t_len, n_env))
    occ <- c(1, 0)

    param_vector <- 0 # intentionally malformed (length 1)

    # Current behavior: early input validation on param_vector
    expect_error(
      loglik_math(param_vector, env_dat, occ),
      regexp = "Assertion on 'param_vector' failed|param_vector.*length",
      fixed = FALSE
    )
  }
)


test_that("function uses non-NA opt values to override internal param_vector", {
  set.seed(1)
  n_loc <- 2L
  t_len <- 2L
  n_env <- 2L
  env_dat <- array(0, dim = c(n_loc, t_len, n_env))
  occ <- c(1, 0)


  param_vector_a <- c(
    mu1 = 13.1,
    mu2 = 5,
    sigltil1 = 0.9,
    sigltil2 = -0.4,
    sigrtil1 = 0.3,
    sigrtil2 = -0.5,
    ctil = -4.4,
    pd = -1.7,
    o_par1 = -9.5
  )

  param_vector_b <- c(
    mu1 = 13.1,
    sigltil1 = 0.9,
    sigltil2 = -0.4,
    sigrtil1 = 0.3,
    sigrtil2 = -0.5,
    ctil = -4.4,
    pd = -1.7,
    o_par1 = -9.5
  )

  # Two opts differing only by mu2 non-NA value
  mask_b <- c(mu2 = 5)


  ll_a <- loglik_math(
    occ = occ,
    env_dat = env_dat,
    param_vector = param_vector_a
  )
  ll_b <- loglik_math(
    occ = occ,
    env_dat = env_dat,
    param_vector = param_vector_b,
    mask = mask_b
  )

  # We expect the same output because in case A
  # There is no mask, and in case b
  # we extract manually mu2
  expect_true(all.equal(ll_a, ll_b))
})
