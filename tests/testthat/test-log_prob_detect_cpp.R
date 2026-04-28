library(testthat)

# Test 1: log_prob_detect_cpp matches R log_prob_detect() on example env_array
test_that("log_prob_detect_cpp matches R log_prob_detect() on example data", {
  env <- examples$env_array
  pl  <- examples$par_list

  r_result <- log_prob_detect(
    env_dat = env,
    mu = pl$mu, sigltil = pl$sigltil, sigrtil = pl$sigrtil,
    o_mat = pl$o_mat, ctil = pl$ctil, pd = pl$pd
  )

  cpp_result <- xsdm:::log_prob_detect_cpp(
    env_dat_vec  = as.vector(env),
    env_dat_dims = as.integer(dim(env)),
    mu      = pl$mu,
    sigltil = pl$sigltil,
    sigrtil = pl$sigrtil,
    o_mat   = pl$o_mat,
    ctil    = pl$ctil,
    pd      = pl$pd,
    return_prob = FALSE
  )

  expect_equal(cpp_result, r_result, tolerance = 1e-10)
})

# Test 2: return_prob = TRUE matches exp() of log result
test_that("log_prob_detect_cpp return_prob=TRUE matches exp(log result)", {
  env <- examples$env_array
  pl  <- examples$par_list

  log_result  <- xsdm:::log_prob_detect_cpp(as.vector(env), as.integer(dim(env)),
                                      pl$mu, pl$sigltil, pl$sigrtil,
                                      pl$o_mat, pl$ctil, pl$pd,
                                      return_prob = FALSE)
  prob_result <- xsdm:::log_prob_detect_cpp(as.vector(env), as.integer(dim(env)),
                                      pl$mu, pl$sigltil, pl$sigrtil,
                                      pl$o_mat, pl$ctil, pl$pd,
                                      return_prob = TRUE)

  expect_equal(prob_result, exp(log_result), tolerance = 1e-12)
})

# Test 3: output length equals number of locations
test_that("log_prob_detect_cpp returns vector of length n_loc", {
  env <- examples$env_array
  pl  <- examples$par_list
  n_loc <- dim(env)[1]

  result <- xsdm:::log_prob_detect_cpp(as.vector(env), as.integer(dim(env)),
                                 pl$mu, pl$sigltil, pl$sigrtil,
                                 pl$o_mat, pl$ctil, pl$pd)
  expect_length(result, n_loc)
})

# Test 4: output values are in valid range for log-probabilities
test_that("log_prob_detect_cpp output is <= 0 (log probabilities)", {
  env <- examples$env_array
  pl  <- examples$par_list
  result <- xsdm:::log_prob_detect_cpp(as.vector(env), as.integer(dim(env)),
                                 pl$mu, pl$sigltil, pl$sigrtil,
                                 pl$o_mat, pl$ctil, pl$pd,
                                 return_prob = FALSE)
  expect_true(all(result <= 0))
})

# Test 5: probabilities are in [0, 1]
test_that("log_prob_detect_cpp probabilities are in [0,1]", {
  env <- examples$env_array
  pl  <- examples$par_list
  result <- xsdm:::log_prob_detect_cpp(as.vector(env), as.integer(dim(env)),
                                 pl$mu, pl$sigltil, pl$sigrtil,
                                 pl$o_mat, pl$ctil, pl$pd,
                                 return_prob = TRUE)
  expect_true(all(result >= 0 & result <= 1))
})

# Test 6: single location works (edge case n_loc = 1)
test_that("log_prob_detect_cpp works for single location", {
  env <- examples$env_array[1, , , drop = FALSE]
  pl  <- examples$par_list
  result <- xsdm:::log_prob_detect_cpp(as.vector(env), as.integer(dim(env)),
                                 pl$mu, pl$sigltil, pl$sigrtil,
                                 pl$o_mat, pl$ctil, pl$pd)
  expect_length(result, 1)
  r_ref <- log_prob_detect(env, pl$mu, pl$sigltil, pl$sigrtil,
                            pl$o_mat, pl$ctil, pl$pd)
  expect_equal(result, r_ref, tolerance = 1e-10)
})
