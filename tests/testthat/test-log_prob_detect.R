library(testthat)

# This test checks that log_prob_detect() works correctly when using mocked data
test_that("log_prob_detect works with realistic env_dat array", {
  # Create a 3D array (2 rows, 5 columns, 5 slices) with rounded integer values
  # This simulates environmental data passed to the function
  env_dat <- array(c(
    # slice 1
    13, 13, 12, 12, 11,
    2, 3, 3, 3, 4,
    # slice 2
    15, 16, 15, 15, 14,
    3, 5, 3, 3, 6,
    # slice 3
    14, 14, 13, 13, 12,
    2, 3, 3, 3, 4,
    # slice 4
    13, 13, 12, 12, 11,
    2, 3, 4, 4, 4,
    # slice 5
    10, 10, 9, 9, 9,
    3, 3, 4, 4, 4
  ), dim = c(2, 5, 5))

  # Define simple numeric parameters for the function
  mu <- c(11, 5)
  sigltil <- c(1, 2)
  sigrtil <- c(2, 1)
  ctil <- -2
  pd <- 0.9
  o_mat <- matrix(c(-0.4, 0.9, -0.9, -0.4), ncol = 2)

  # with_mocked_bindings() temporarily replaces a function with a mock version
  # Here, we mock like_neg_ltsgr_cpp() to always return 2.0
  # This isolates log_prob_detect() from its C++ dependency for predictable
  # testing
  with_mocked_bindings(
    `like_neg_ltsgr_cpp` = function(env_dat, mu, sigltil, sigrtil, o_mat) 2.0,
    {
      # Call the function with return_prob = FALSE
      # Should return log probability
      result_log <- log_prob_detect(env_dat,
        mu,
        sigltil,
        sigrtil,
        o_mat,
        ctil,
        pd,
        return_prob = FALSE
      )
      # Compute expected value manually using the formula in the function
      expected_log <- log(pd) - log1pexp(ctil + 2.0)

      # Check that the result matches the expected value
      expect_equal(result_log, expected_log)

      # Call the function with return_prob = TRUE (should return probability)
      result_prob <- log_prob_detect(env_dat,
        mu,
        sigltil,
        sigrtil,
        o_mat,
        ctil,
        pd,
        return_prob = TRUE
      )

      # Check that the result matches exp(expected_log)
      expect_equal(result_prob, exp(expected_log))
    }
  )
})

# Second test: checks edge cases for pd (probability of detection)
test_that("log_prob_detect handles edge cases for pd", {
  # Use a simpler env_dat array for this test
  env_dat <- array(1:50, dim = c(2, 5, 5))

  # Same parameters as before
  mu <- c(11, 5)
  sigltil <- c(1, 2)
  sigrtil <- c(2, 1)
  ctil <- -2
  o_mat <- matrix(c(-0.4, 0.9, -0.9, -0.4), ncol = 2)

  # Mock like_neg_ltsgr_cpp()
  with_mocked_bindings(
    `like_neg_ltsgr_cpp` = function(env_dat, mu, sigltil, sigrtil, o_mat) 2.0,
    {
      # pd = 1 (maximum detection probability)
      expect_equal(
        log_prob_detect(env_dat, mu, sigltil, sigrtil, o_mat, ctil, 1),
        log(1) - log1pexp(ctil + 2.0)
      )

      # pd = 0.5 (50% detection probability)
      expect_equal(
        log_prob_detect(env_dat, mu, sigltil, sigrtil, o_mat, ctil, 0.5),
        log(0.5) - log1pexp(ctil + 2.0)
      )
    }
  )
})

test_that("log_prob_detect returns -Inf (or 0) when pd = 0", {
  # Minimal env_dat array
  env_dat <- array(1:50, dim = c(2, 5, 5))
  mu <- c(11, 5)
  sigltil <- c(1, 2)
  sigrtil <- c(2, 1)
  ctil <- -2
  o_mat <- matrix(c(-0.4, 0.9, -0.9, -0.4), ncol = 2)
  
  # Mock like_neg_ltsgr_cpp() to return a finite value
  with_mocked_bindings(
    `like_neg_ltsgr_cpp` = function(env_dat, mu, sigltil, sigrtil, o_mat) 2.0,
    {
      # Case: return_prob = FALSE -> should give -Inf
      result_log <- log_prob_detect(
        env_dat, mu, sigltil, sigrtil, o_mat, ctil, pd = 0,
        return_prob = FALSE
      )
      expect_true(all(is.infinite(result_log) & result_log < 0))
      
      # Case: return_prob = TRUE -> should give 0
      result_prob <- log_prob_detect(
        env_dat, mu, sigltil, sigrtil, o_mat, ctil, pd = 0,
        return_prob = TRUE
      )
      expect_equal(result_prob, rep(0, length(result_prob)))
    }
  )
})
