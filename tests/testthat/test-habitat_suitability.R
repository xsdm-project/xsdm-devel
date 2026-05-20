library(testthat)

# ---------------------------------------------------------------------------
# habitat_suitability: contract & parity with vsp() on example data
# ---------------------------------------------------------------------------


# Test 2: returns SpatRaster with correct dimensions
test_that("habitat_suitability returns SpatRaster with 1 layer", {
  skip_if_not_installed("terra")
  bio1  <- terra::unwrap(example_1$bio01) / 100
  bio12 <- terra::unwrap(example_1$bio12) / 100

  out <- habitat_suitability(
    example_1$true_par_list,
    list(bio1, bio12),
    return_prob = TRUE
  )

  expect_equal(terra::nlyr(out), 1L)
  expect_equal(terra::nrow(out), terra::nrow(bio1))
  expect_equal(terra::ncol(out), terra::ncol(bio1))
})

# Test 2: log-probability output is <= 0
test_that("habitat_suitability log-prob output is <= 0", {
  skip_if_not_installed("terra")
  bio1  <- terra::unwrap(example_1$bio01) / 100
  bio12 <- terra::unwrap(example_1$bio12) / 100

  out <- habitat_suitability(
    example_1$true_par_list,
    list(bio1, bio12),
    return_prob = FALSE
  )

  vals <- terra::values(out)
  # allow NA for cells that are NA in input
  expect_true(all(vals[!is.na(vals)] <= 0))
})

# Test 3: probability output is in [0, 1]
test_that("habitat_suitability probabilities are in [0, 1]", {
  skip_if_not_installed("terra")
  bio1  <- terra::unwrap(example_1$bio01) / 100
  bio12 <- terra::unwrap(example_1$bio12) / 100

  out <- habitat_suitability(
    example_1$true_par_list,
    list(bio1, bio12),
    return_prob = TRUE
  )

  vals <- terra::values(out)
  expect_true(all(vals[!is.na(vals)] >= 0 & vals[!is.na(vals)] <= 1))
})

# Test 4: write to file works
test_that("habitat_suitability writes to file", {
  skip_if_not_installed("terra")
  bio1  <- terra::unwrap(example_1$bio01) / 100
  bio12 <- terra::unwrap(example_1$bio12) / 100

  tmp <- tempfile(fileext = ".tif")
  habitat_suitability(
    example_1$true_par_list,
    list(bio1, bio12),
    output    = tmp,
    overwrite = TRUE
  )

  expect_true(file.exists(tmp))
  from_disk <- terra::rast(tmp)
  expect_equal(terra::nlyr(from_disk), 1L)
  unlink(tmp)
})

# Test 5: geometry mismatch is caught
test_that("habitat_suitability errors on geometry mismatch", {
  skip_if_not_installed("terra")
  bio1  <- terra::unwrap(example_1$bio01) / 100
  bio12 <- terra::unwrap(example_1$bio12) / 100
  # crop bio12 to a different extent
  bio12_crop <- terra::crop(bio12, terra::ext(bio12) * 0.5)

  expect_error(
    habitat_suitability(example_1$true_par_list, list(bio1, bio12_crop)),
    regexp = "geometry"
  )
})

# Test 6: log and probability outputs are consistent via exp()
test_that("habitat_suitability return_prob=FALSE is log of return_prob=TRUE", {
  skip_if_not_installed("terra")
  bio1  <- terra::unwrap(example_1$bio01) / 100
  bio12 <- terra::unwrap(example_1$bio12) / 100

  hs_log  <- habitat_suitability(example_1$true_par_list, list(bio1, bio12),
                                 return_prob = FALSE)
  hs_prob <- habitat_suitability(example_1$true_par_list, list(bio1, bio12),
                                 return_prob = TRUE)

  v_log  <- as.vector(terra::values(hs_log))
  v_prob <- as.vector(terra::values(hs_prob))
  keep <- !is.na(v_log) & !is.na(v_prob)
  expect_equal(exp(v_log[keep]), v_prob[keep], tolerance = 1e-10)
})

# Test 7: output layer name reflects return_prob
test_that("habitat_suitability layer name reflects return_prob", {
  skip_if_not_installed("terra")
  bio1  <- terra::unwrap(example_1$bio01) / 100
  bio12 <- terra::unwrap(example_1$bio12) / 100

  r_prob <- habitat_suitability(example_1$true_par_list, list(bio1, bio12),
                                return_prob = TRUE)
  r_log  <- habitat_suitability(example_1$true_par_list, list(bio1, bio12),
                                return_prob = FALSE)

  expect_equal(terra::names(r_prob), "habitat_suitability")
  expect_equal(terra::names(r_log),  "log_prob_detect")
})

# Test 8: output geometry matches input
test_that("habitat_suitability output geometry matches input", {
  skip_if_not_installed("terra")
  bio1  <- terra::unwrap(example_1$bio01) / 100
  bio12 <- terra::unwrap(example_1$bio12) / 100
  env_list <- list(bio1, bio12)

  result <- habitat_suitability(example_1$true_par_list, env_list)
  ref    <- env_list[[1]]

  expect_equal(terra::nrow(result), terra::nrow(ref))
  expect_equal(terra::ncol(result), terra::ncol(ref))
  expect_equal(as.vector(terra::ext(result)),  as.vector(terra::ext(ref)))
  expect_equal(terra::crs(result),  terra::crs(ref))
})

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

test_that("habitat_suitability fails with non-SpatRaster env_list", {
  expect_error(
    habitat_suitability(
      param_list = list(mu = 1, sigltil = 1, sigrtil = 1,
                        o_mat = matrix(1, 1, 1), ctil = 0, pd = 0.5),
      env_list   = list(matrix(1:9, 3, 3))
    ),
    regexp = "SpatRaster"
  )
})

test_that("habitat_suitability fails with missing param_list keys", {
  skip_if_not_installed("terra")
  bio1 <- terra::unwrap(example_1$bio01) / 100
  bad_params <- list(mu = 1)
  expect_error(
    habitat_suitability(bad_params, list(bio1)),
    regexp = "param_list must contain"
  )
})
