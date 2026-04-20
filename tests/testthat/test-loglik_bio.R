library(testthat)
test_that("loglik_bio: full coverage including uncovered lines", {
  # Minimal valid data with n_env = 2 (matches internal checks)
  n_loc <- 2
  t_len <- 2
  n_env <- 2
  env_dat <- array(0, dim = c(n_loc, t_len, n_env))
  occ <- c(1, 0)
  # length matches n_env
  mu <- c(1, 1)
  sigl <- c(1, 1)
  sigr <- c(1, 1)
  ctil <- 0.1
  pd <- 0.5
  # 2 x 2 orthogonal matrix
  o_mat <- diag(n_env)

  # ---- Positive case: covers input validation and default thread options ----
  expect_silent(
    loglik_bio(env_dat, occ, mu, sigl, sigr, o_mat, ctil, pd)
  )

  # ---- Covers num_threads argument (thread options lines) ----
  expect_silent(
    loglik_bio(env_dat, occ, mu, sigl, sigr, o_mat, ctil, pd, num_threads = 1L)
  )

  # ---- Covers sum_log_p = FALSE branch ----
  res_vec <- loglik_bio(env_dat, occ, mu, sigl, sigr, o_mat, ctil, pd,
    sum_log_p = FALSE, return_prob = FALSE
  )
  expect_length(res_vec, length(occ))


  # ---- Covers return_prob = TRUE branch ----
  res_prob <- loglik_bio(env_dat,
    occ,
    mu,
    sigl,
    sigr,
    o_mat,
    ctil,
    pd,
    sum_log_p = TRUE,
    return_prob = TRUE
  )

  expect_true(length(res_prob) == 1 && is.numeric(res_prob) && !is.na(res_prob))


  # ---- Negative cases: trigger each assertion error ----
  expect_error(
    loglik_bio(1:5, occ, mu, sigl, sigr, o_mat, ctil, pd),
    regexp = "env_dat"
  )
  expect_error(
    loglik_bio(env_dat, c(-1, 0), mu, sigl, sigr, o_mat, ctil, pd),
    regexp = "occ"
  )
  expect_error(
    loglik_bio(env_dat, occ, mu = NA_real_, sigl, sigr, o_mat, ctil, pd),
    regexp = "mu"
  )
  expect_error(
    loglik_bio(env_dat, occ, mu, sigl = NA_real_, sigr, o_mat, ctil, pd),
    regexp = "sigl"
  )
  expect_error(
    loglik_bio(env_dat, occ, mu, sigl, sigr = NA_real_, o_mat, ctil, pd),
    regexp = "sigr"
  )
  expect_error(
    loglik_bio(env_dat, occ, mu, sigl, sigr, o_mat, ctil = c(0.1, 0.2), pd),
    regexp = "ctil"
  )
  expect_error(
    loglik_bio(env_dat, occ, mu, sigl, sigr, o_mat, ctil, pd = c(0.9, 0.8)),
    regexp = "pd"
  )
  expect_error(
    loglik_bio(env_dat, occ, mu, sigl, sigr, o_mat = 1:4, ctil, pd),
    regexp = "o_mat"
  )
})

# ---------------------------------------------------------------------------
# Tests for sigltil / sigrtil Inf cases
# ---------------------------------------------------------------------------

test_that("loglik_bio returns finite numeric when sigltil = Inf", {
  n_loc <- 2
  t_len <- 2
  n_env <- 2
  env_dat <- array(runif(n_loc * t_len * n_env), dim = c(n_loc, t_len, n_env))
  occ <- c(1, 0)
  mu <- c(1, 1)
  sigltil <- c(Inf, Inf)
  sigrtil <- c(1, 1)
  ctil <- 0.1
  pd <- 0.5
  o_mat <- diag(n_env)

  result <- loglik_bio(env_dat, occ, mu, sigltil, sigrtil, o_mat, ctil, pd,
                       num_threads = 1L)
  expect_true(is.numeric(result))
  expect_length(result, 1L)
  expect_true(is.finite(result))

  # Also check vector output (sum_log_p = FALSE)
  res_vec <- loglik_bio(env_dat, occ, mu, sigltil, sigrtil, o_mat, ctil, pd,
                        sum_log_p = FALSE, num_threads = 1L)
  expect_length(res_vec, n_loc)
  expect_true(all(is.finite(res_vec)))
})

test_that("loglik_bio returns finite numeric when sigrtil = Inf", {
  n_loc <- 2
  t_len <- 2
  n_env <- 2
  env_dat <- array(runif(n_loc * t_len * n_env), dim = c(n_loc, t_len, n_env))
  occ <- c(1, 0)
  mu <- c(1, 1)
  sigltil <- c(1, 1)
  sigrtil <- c(Inf, Inf)
  ctil <- 0.1
  pd <- 0.5
  o_mat <- diag(n_env)

  result <- loglik_bio(env_dat, occ, mu, sigltil, sigrtil, o_mat, ctil, pd,
                       num_threads = 1L)
  expect_true(is.numeric(result))
  expect_length(result, 1L)
  expect_true(is.finite(result))

  # Also check vector output (sum_log_p = FALSE)
  res_vec <- loglik_bio(env_dat, occ, mu, sigltil, sigrtil, o_mat, ctil, pd,
                        sum_log_p = FALSE, num_threads = 1L)
  expect_length(res_vec, n_loc)
  expect_true(all(is.finite(res_vec)))
})

test_that("loglik_bio returns finite numeric when both sigltil and sigrtil = Inf", {
  n_loc <- 2
  t_len <- 2
  n_env <- 2
  env_dat <- array(runif(n_loc * t_len * n_env), dim = c(n_loc, t_len, n_env))
  occ <- c(1, 0)
  mu <- c(1, 1)
  sigltil <- c(Inf, Inf)
  sigrtil <- c(Inf, Inf)
  ctil <- 0.1
  pd <- 0.5
  o_mat <- diag(n_env)

  result <- loglik_bio(env_dat, occ, mu, sigltil, sigrtil, o_mat, ctil, pd,
                       num_threads = 1L)
  expect_true(is.numeric(result))
  expect_length(result, 1L)
  expect_true(is.finite(result))

  # Also check vector output (sum_log_p = FALSE)
  res_vec <- loglik_bio(env_dat, occ, mu, sigltil, sigrtil, o_mat, ctil, pd,
                        sum_log_p = FALSE, num_threads = 1L)
  expect_length(res_vec, n_loc)
  expect_true(all(is.finite(res_vec)))
})

test_that("loglik_bio handles mixed finite/Inf sigltil and sigrtil (p=2)", {
  n_loc <- 2
  t_len <- 2
  n_env <- 2
  env_dat <- array(runif(n_loc * t_len * n_env), dim = c(n_loc, t_len, n_env))
  occ <- c(1, 0)
  mu <- c(1, 1)
  # var1: left=Inf (no left constraint); var2: right=Inf (no right constraint)
  sigltil <- c(Inf, 1)
  sigrtil <- c(1, Inf)
  ctil <- 0.1
  pd <- 0.5
  o_mat <- diag(n_env)

  result <- loglik_bio(env_dat, occ, mu, sigltil, sigrtil, o_mat, ctil, pd,
                       num_threads = 1L)
  expect_true(is.numeric(result))
  expect_length(result, 1L)
  expect_true(is.finite(result))

  # Also check vector output (sum_log_p = FALSE)
  res_vec <- loglik_bio(env_dat, occ, mu, sigltil, sigrtil, o_mat, ctil, pd,
                        sum_log_p = FALSE, num_threads = 1L)
  expect_length(res_vec, n_loc)
  expect_true(all(is.finite(res_vec)))
})

test_that("loglik_bio handles Inf sigltil/sigrtil for p=1", {
  n_loc <- 3
  t_len <- 2
  n_env <- 1
  env_dat <- array(runif(n_loc * t_len * n_env), dim = c(n_loc, t_len, n_env))
  occ <- c(1, 0, 1)
  mu <- 0.5
  o_mat <- matrix(1, 1, 1)
  ctil <- 0.1
  pd <- 0.5

  # sigltil = Inf only
  result_l <- loglik_bio(env_dat, occ, mu, sigltil = Inf, sigrtil = 1,
                         o_mat, ctil, pd, num_threads = 1L)
  expect_true(is.numeric(result_l) && is.finite(result_l))

  # sigrtil = Inf only
  result_r <- loglik_bio(env_dat, occ, mu, sigltil = 1, sigrtil = Inf,
                         o_mat, ctil, pd, num_threads = 1L)
  expect_true(is.numeric(result_r) && is.finite(result_r))

  # Both Inf
  result_both <- loglik_bio(env_dat, occ, mu, sigltil = Inf, sigrtil = Inf,
                            o_mat, ctil, pd, num_threads = 1L)
  expect_true(is.numeric(result_both) && is.finite(result_both))
})
