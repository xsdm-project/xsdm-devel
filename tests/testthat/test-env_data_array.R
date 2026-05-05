library(testthat)
library(terra)

test_that("returns array when occ is NULL and one env layer", {
  # Create mock environmental data
  r1 <- rast(nrows = 2, ncols = 2)
  values(r1) <- 1:4
  r2 <- rast(nrows = 2, ncols = 2)
  values(r2) <- 5:8

  env_data <- list(bio1 = r1, bio12 = r2)

  # Create mock occurrence data
  occ <- data.frame(
    name = c("site1", "site2"),
    x = c(0.5, 1.5),
    y = c(0.5, 1.5),
    presence = c(1, 0)
  )
  result <- env_data_array(list(bio1 = r1), occ = NULL)
  expect_true(is.array(result))
  expect_equal(nrow(result), ncell(r1))
})

test_that("returns array when occ is NULL and multiple env layers", {
  # Create mock environmental data
  r1 <- rast(nrows = 2, ncols = 2)
  values(r1) <- 1:4
  r2 <- rast(nrows = 2, ncols = 2)
  values(r2) <- 5:8

  env_data <- list(bio1 = r1, bio12 = r2)

  # Create mock occurrence data
  occ <- data.frame(
    name = c("site1", "site2"),
    x = c(0.5, 1.5),
    y = c(0.5, 1.5),
    presence = c(1, 0)
  )
  result <- env_data_array(env_data, occ = NULL)
  expect_true(is.array(result))
  expect_equal(length(dim(result)), 3) # 3D array
})

test_that("returns array when occ is provided and one env layer", {
  # Create mock environmental data
  r1 <- rast(nrows = 2, ncols = 2)
  values(r1) <- 1:4
  r2 <- rast(nrows = 2, ncols = 2)
  values(r2) <- 5:8

  env_data <- list(bio1 = r1, bio12 = r2)

  # Create mock occurrence data
  occ <- data.frame(
    name = c("site1", "site2"),
    lon = c(0.5, 1.5),
    lat = c(0.5, 1.5),
    presence = c(1, 0)
  )
  result <- env_data_array(list(bio1 = r1), occ = occ)
  expect_true(is.array(result))
  expect_equal(nrow(result), nrow(occ))
})

test_that("returns array when occ is provided and multiple env layers", {
  # Create mock environmental data
  r1 <- rast(nrows = 2, ncols = 2)
  values(r1) <- 1:4
  r2 <- rast(nrows = 2, ncols = 2)
  values(r2) <- 5:8

  env_data <- list(bio1 = r1, bio12 = r2)

  # Create mock occurrence data
  occ <- data.frame(
    name = c("site1", "site2"),
    lon = c(0.5, 1.5),
    lat = c(0.5, 1.5),
    presence = c(1, 0)
  )
  result <- env_data_array(env_data, occ = occ)
  expect_true(is.array(result))
  expect_equal(length(dim(result)), 3)
})

test_that("fails when occ is missing required columns", {
  # Create mock environmental data
  r1 <- rast(nrows = 2, ncols = 2)
  values(r1) <- 1:4
  r2 <- rast(nrows = 2, ncols = 2)
  values(r2) <- 5:8

  env_data <- list(bio1 = r1, bio12 = r2)

  # Create mock occurrence data
  occ <- data.frame(
    name = c("site1", "site2"),
    lon = c(0.5, 1.5),
    lat = c(0.5, 1.5),
    presence = c(1, 0)
  )
  bad_occ <- data.frame(longitude = 1, latitude = 2)
  expect_error(env_data_array(env_data, bad_occ), "must.include")
})

test_that("fails when env_data is not a list", {
  # Create mock environmental data
  r1 <- rast(nrows = 2, ncols = 2)
  values(r1) <- 1:4
  r2 <- rast(nrows = 2, ncols = 2)
  values(r2) <- 5:8

  env_data <- list(bio1 = r1, bio12 = r2)

  # Create mock occurrence data
  occ <- data.frame(
    name = c("site1", "site2"),
    lon = c(0.5, 1.5),
    lat = c(0.5, 1.5),
    presence = c(1, 0)
  )
  expect_error(env_data_array(r1, occ), "list")
})
