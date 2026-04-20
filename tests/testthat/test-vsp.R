library(testthat)
library(terra)
library(tibble)
library(checkmate)

test_that("vsp returns tibble when return_raster = FALSE", {
  # Create mock environmental data
  r1 <- rast(nrows = 2, ncols = 2, vals = c(10, 12, 14, 16))
  r2 <- rast(nrows = 2, ncols = 2, vals = c(100, 120, 140, 160))
  env_data <- list(bio1 = r1, bio12 = r2)

  # Mock parameter list
  param_list <- list(
    mu = c(10, 100),
    sigltil = c(2, 10),
    sigrtil = c(2, 10),
    ctil = c(0.5, 0.5),
    pd = 0.8,
    o_mat = matrix(0, nrow = 2, ncol = 2)
  )

  # Run function
  result <- vsp(env_data, param_list, return_raster = FALSE)

  # Assertions
  expect_s3_class(result, "tbl_df")
  expect_true(all(c("x", "y", "probs") %in% colnames(result)))
  expect_equal(nrow(result), ncell(r1))
})

test_that("vsp returns SpatRaster when return_raster = TRUE", {
  r1 <- rast(nrows = 2, ncols = 2, vals = c(10, 12, 14, 16))
  r2 <- rast(nrows = 2, ncols = 2, vals = c(100, 120, 140, 160))
  env_data <- list(bio1 = r1, bio12 = r2)

  param_list <- list(
    mu = c(10, 100),
    sigltil = c(2, 10),
    sigrtil = c(2, 10),
    ctil = c(0.5, 0.5),
    pd = 0.8,
    o_mat = matrix(0, nrow = 2, ncol = 2)
  )

  result <- vsp(env_data, param_list, return_raster = TRUE)

  expect_s4_class(result, "SpatRaster")
  expect_equal(ncell(result), ncell(r1))
})

# ---- Invalid input tests ----
test_that("vsp fails with invalid env_data", {
  bad_env <- list(matrix(1:4, 2, 2)) # Not SpatRaster
  param_list <- list(
    mu = 1,
    sigltil = 1,
    sigrtil = 1,
    ctil = 1,
    pd = 1,
    o_mat = matrix(0, 1, 1)
  )
  expect_error(vsp(bad_env, param_list, return_raster = FALSE), "SpatRaster")
})


test_that("vsp fails with missing parameters", {
  r1 <- rast(nrows = 2, ncols = 2, vals = c(10, 12, 14, 16))
  env_data <- list(bio1 = r1)

  # Missing required names
  bad_params <- list(mu = 1)
  expect_error(
    vsp(env_data, bad_params, return_raster = FALSE),
    "param_list must contain"
  )
})

test_that("vsp fails with invalid return_raster type", {
  r1 <- rast(nrows = 2, ncols = 2, vals = c(10, 12, 14, 16))
  env_data <- list(bio1 = r1)
  param_list <- list(
    mu = 1,
    sigltil = 1,
    sigrtil = 1,
    ctil = 1,
    pd = 1,
    o_mat = matrix(0, 1, 1)
  )

  expect_error(
    vsp(env_data, param_list, return_raster = "yes"),
    regexp = "logical flag"
  )
})
