library(testthat)

test_that("like_neg_ltsgr_cpp covers p == 1 branch with real like_ltsg", {
  # Prepare inputs for p = 1
  env_dat <- array(runif(3 * 4 * 1), dim = c(3, 4, 1)) # n=1, ts_length=4, p=1
  mu <- 0.7
  sigl <- 1.1
  sigr <- 1.3
  o_mat <- matrix(1, nrow = 1, ncol = 1)

  # Call the actual function
  result <- xsdm:::like_neg_ltsgr_cpp(env_dat, mu, sigl, sigr, o_mat, num_threads = 1)

  # Assertions
  expect_type(result, "double")
  # matches n
  expect_equal(length(result), dim(env_dat)[1])
  # no NA or Inf
  expect_true(all(is.finite(result)))
  # Likelihood should be non-negative
  expect_true(all(result >= 0))
})

test_that("like_neg_ltsgr_cpp handles Inf entries in sigltil and sigrtil", {
  # ---- 1D case with explicit array construction ----
  # 2 locations, 2 time steps, 1 variable
  # Location 1: values (1, 2) ; Location 2: values (3, 4)
  env_dat <- array(NA, dim = c(2, 2, 1))
  env_dat[1, , 1] <- c(1, 2)
  env_dat[2, , 1] <- c(3, 4)
  
  mu <- 0
  o_mat <- matrix(1, 1, 1)
  
  # Case 1: sigltil = Inf, sigrtil finite
  sigl <- Inf
  sigr <- 1
  res <- xsdm:::like_neg_ltsgr_cpp(env_dat, mu, sigl, sigr, o_mat, num_threads = 1)
  # Computation:
  # - dl_inv = 0, dr_inv = 1
  # - usym = max(0, dot_product) * 1
  # - loc1 dot_product = (1,2) -> squares = 1+4=5 -> / (2*q=4) = 1.25
  # - loc2 dot_product = (3,4) -> squares = 9+16=25 -> /4 = 6.25
  expect_equal(res, c(1.25, 6.25))
  
  
  # Case 2: sigrtil = Inf, sigltil finite
  sigl <- 1
  sigr <- Inf
  res <- xsdm:::like_neg_ltsgr_cpp(env_dat, mu, sigl, sigr, o_mat, num_threads = 1)
  # Only negative deviations matter, all positive => zero
  expect_equal(res, c(0, 0))
  
  # Case 3: both Inf -> zero
  sigl <- Inf
  sigr <- Inf
  res <- xsdm:::like_neg_ltsgr_cpp(env_dat, mu, sigl, sigr, o_mat, num_threads = 1)
  expect_equal(res, c(0, 0))
  
  # ---- 2D case with mixed Inf ----
  # 2 locations, 2 time steps, 2 variables
  # For clarity, fill each location's time series explicitly
  env_dat_2d <- array(NA, dim = c(2, 2, 2))
  # Location 1: var1 = (1,2), var2 = (1,2)
  env_dat_2d[1, , 1] <- c(1, 2)
  env_dat_2d[1, , 2] <- c(1, 2)
  # Location 2: var1 = (3,4), var2 = (3,4)
  env_dat_2d[2, , 1] <- c(3, 4)
  env_dat_2d[2, , 2] <- c(3, 4)
  
  mu <- c(0, 0)
  o_mat <- diag(2)
  sigl <- c(1, Inf)   # var1: left=1; var2: left=Inf -> effectively only right matters
  sigr <- c(Inf, 1)   # var1: right=Inf -> left only; var2: right=1 -> right only
  
  res <- xsdm:::like_neg_ltsgr_cpp(env_dat_2d, mu, sigl, sigr, o_mat, num_threads = 1)
  # var1: only left branch (negative deviations) -> all positive => 0 contribution
  # var2: only right branch (positive deviations) -> values = (1,2) for loc1, (3,4) for loc2
  # => squares sum = 5 for loc1, 25 for loc2, divided by (2*2) = 1.25 and 6.25
  expect_equal(res, c(1.25, 6.25))
})
