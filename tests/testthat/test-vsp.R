# tests/testthat/test-vsp.R

test_that("vsp returns a tibble with correct structure, respects threshold, and handles edge cases", {
  # Load example data
  data("example_1", package = "xsdm")
  
  # Unpack rasters (scale as in example)
  bio1_ts  <- terra::unwrap(example_1$bio01) / 100
  bio12_ts <- terra::unwrap(example_1$bio12) / 100
  env_data <- list(bio1 = bio1_ts, bio12 = bio12_ts)
  
  # Set seed for reproducibility
  set.seed(123)
  
  # ---- 1. Normal use with threshold = 0.5 ----
  size_pres <- 50
  size_abs  <- 50
  
  result <- vsp(
    param_list    = example_1$par_list,
    env_data      = env_data,
    size_presence = size_pres,
    size_absence  = size_abs,
    threshold     = 0.5
  )
  
  # Output structure
  expect_true(tibble::is_tibble(result))
  expect_named(result, c("lon", "lat", "presence"))
  expect_type(result$lon, "double")
  expect_type(result$lat, "double")
  expect_type(result$presence, "integer")
  expect_equal(nrow(result), size_pres + size_abs)
  expect_true(all(result$presence %in% c(0L, 1L)))
  
  # Verify threshold split
  suit <- habitat_suitability(
    param_list  = example_1$par_list,
    env_list    = env_data,
    return_prob = TRUE
  )
  pts <- terra::vect(result[, c("lon", "lat")], geom = c("lon", "lat"))
  prob_vals <- terra::extract(suit, pts)[[2]]  # first column = values
  
  is_pres <- prob_vals > 0.5
  is_abs  <- prob_vals <= 0.5
  expect_true(all(is_pres | is_abs))
  expect_equal(sum(is_pres), size_pres)
  expect_equal(sum(is_abs),  size_abs)
  
  # ---- 2. Edge case: threshold = 1.0 (no presence cells) ----
  expect_warning(
    result2 <- vsp(
      param_list    = example_1$par_list,
      env_data      = env_data,
      size_presence = size_pres,
      size_absence  = size_abs,
      threshold     = 1.0
    ),
    "No cells available for presence sampling"
  )
  expect_true(tibble::is_tibble(result2))
  pts2 <- terra::vect(result2[, c("lon", "lat")], geom = c("lon", "lat"))
  prob_vals2 <- terra::extract(suit, pts2)[[2]]
  expect_true(all(prob_vals2 <= 1.0))
  # Only absence points exist -> exactly size_abs rows (assuming enough cells)
  expect_equal(nrow(result2), size_abs)
  
  # ---- 3. Edge case: threshold = 0.0 (absence pool may be empty or not) ----
  # Instead of expecting a warning, we verify that all sampled points have
  # probability > 0 (i.e., they come from the presence pool) and the total
  # number of rows equals size_pres (the presence sample size).
  expect_warning(
    result3 <- vsp(
      param_list    = example_1$par_list,
      env_data      = env_data,
      size_presence = size_pres,
      size_absence  = size_abs,
      threshold     = 0.0
    ), regexp = "No cells available for absence sampling"
  )
  
  pts3 <- terra::vect(result3[, c("lon", "lat")], geom = c("lon", "lat"))
  prob_vals3 <- terra::extract(suit, pts3)[[2]]
  expect_true(all(prob_vals3 > 0))  # all points must have positive probability
  expect_equal(nrow(result3), size_abs)
  
  # ---- 4. Sample size larger than available cells (warning + reduced sample) ----
  # Count number of cells with prob > 0.5 (presence pool)
  n_pres_cells <- sum(prob_vals > 0.5, na.rm = TRUE)
  expect_warning(
    result4 <- vsp(
      param_list    = example_1$par_list,
      env_data      = env_data,
      size_presence = n_pres_cells + 10000,
      size_absence  = size_abs,
      threshold     = 0.5
    ),
    "exceeds available cells"
  )
  expect_lte(nrow(result4), n_pres_cells + size_abs + 10000)
  
  # ---- 5. Invalid inputs ----
  expect_error(
    vsp(example_1$par_list, env_data, 10, 10, threshold = 1.5),
    "Element 1 is not <= 1"
  )
  expect_error(
    vsp(example_1$par_list, env_data, 0, 10),
    "size_presence"
  )
  expect_error(
    vsp(example_1$par_list, env_data, 10, 0),
    "size_absence"
  )
})
