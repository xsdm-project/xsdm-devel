# Test to ensure the C++ implementation (like_neg_ltsgr_cpp) and
# the  internal R implementation (xsdm:::like_neg_ltsgr_r) return
# identical results when fed parameters from `param_table_example`,
# using an orthogonal matrix built via `build_orthogonal_matrix`.

library(testthat)

test_that("R and C++ likelihood workers are numerically identical on param_table_example", {
  # Load packaged parameter table
  data("param_table_example", package = "xsdm")
  expect_true(exists("param_table_example"))

  set.seed(123)

  # Dimensions for the synthetic, deterministic environmental data
  n_locations <- 3
  ts_length <- 4
  p <- 2

  # Choose several rows from the parameter table (only those that exist)
  rows_to_test <- c(1, 2, 3, 10, 20, 30, 40, 50)
  rows_to_test <- rows_to_test[rows_to_test <= nrow(param_table_example)]

  for (rr in rows_to_test) {
    row <- param_table_example[rr, ]
    param_list <- math_to_bio(unlist(row))
    # Parameters from the table
    mu <- c(param_list$mu)
    sigltil <- c(param_list$sigltil)
    sigrtil <- c(param_list$sigrtil)
    ctil <- c(param_list$ctil)
    pd <- c(param_list$pd)

    # Build 2x2 orthogonal matrix from the single entry (skew-symmetric exp)
    # Use the package's function directly
    o_mat <- param_list$o_mat


    # Compute outputs
    r_out <- xsdm:::like_neg_ltsgr_r(example_1_env_array, mu, sigltil, sigrtil, o_mat)
    cpp_out <- like_neg_ltsgr_cpp(example_1_env_array, mu, sigltil, sigrtil, o_mat)

    # Output should be a numeric vector of length equal to number of locations
    expect_type(r_out, "double")
    expect_type(cpp_out, "double")
    expect_equal(length(r_out), dim(example_1_env_array)[1])
    expect_equal(length(cpp_out), dim(example_1_env_array)[1])

    # Parity check with tight tolerance
    expect_equal(cpp_out, r_out,
      tolerance = 1e-14,
      info = paste("Mismatch at param_table_example row", rr)
    )
  }
})

test_that("R and C++ implementations agree with Inf in sigltil/sigrtil", {
  # ---- 1D case (p = 1) ----
  # Build array: 2 locations, 2 time steps, 1 variable
  env_dat_1d <- array(NA, dim = c(2, 2, 1))
  env_dat_1d[1, , 1] <- c(1, 2)   # location 1: values (1,2)
  env_dat_1d[2, , 1] <- c(3, 4)   # location 2: values (3,4)
  mu <- 0
  o_mat <- matrix(1, 1, 1)
  
  # Case 1: sigltil = Inf, sigrtil finite
  sigl <- Inf
  sigr <- 1
  r1 <- xsdm:::like_neg_ltsgr_r(env_dat_1d, mu, sigl, sigr, o_mat)
  cpp1 <- like_neg_ltsgr_cpp(env_dat_1d, mu, sigl, sigr, o_mat, num_threads = 1)
  expect_equal(cpp1, r1, tolerance = 1e-14)
  
  # Case 2: sigrtil = Inf, sigltil finite
  sigl <- 1
  sigr <- Inf
  r2 <- xsdm:::like_neg_ltsgr_r(env_dat_1d, mu, sigl, sigr, o_mat)
  cpp2 <- like_neg_ltsgr_cpp(env_dat_1d, mu, sigl, sigr, o_mat, num_threads = 1)
  expect_equal(cpp2, r2, tolerance = 1e-14)
  
  # Case 3: both Inf
  sigl <- Inf
  sigr <- Inf
  r3 <- xsdm:::like_neg_ltsgr_r(env_dat_1d, mu, sigl, sigr, o_mat)
  cpp3 <- like_neg_ltsgr_cpp(env_dat_1d, mu, sigl, sigr, o_mat, num_threads = 1)
  expect_equal(cpp3, r3, tolerance = 1e-14)
  
  # ---- 2D case (p = 2) ----
  # Build array: 2 locations, 2 time steps, 2 variables
  env_dat_2d <- array(NA, dim = c(2, 2, 2))
  env_dat_2d[1, , 1] <- c(1, 2)   # loc1 var1
  env_dat_2d[1, , 2] <- c(1, 2)   # loc1 var2
  env_dat_2d[2, , 1] <- c(3, 4)   # loc2 var1
  env_dat_2d[2, , 2] <- c(3, 4)   # loc2 var2
  mu <- c(0, 0)
  o_mat <- diag(2)
  
  # Mixed Inf: var1 only left, var2 only right
  sigl <- c(1, Inf)
  sigr <- c(Inf, 1)
  r4 <- xsdm:::like_neg_ltsgr_r(env_dat_2d, mu, sigl, sigr, o_mat)
  cpp4 <- like_neg_ltsgr_cpp(env_dat_2d, mu, sigl, sigr, o_mat, num_threads = 1)
  expect_equal(cpp4, r4, tolerance = 1e-14)
})
