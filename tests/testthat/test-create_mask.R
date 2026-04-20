library(testthat)

test_that("create_mask returns empty mask when params is NULL", {
  result <- create_mask(p = 2)

  # Expected names from make_mask_names(2)
  expected_names <- c(
    "mu1", "mu2", "sigltil1", "sigltil2",
    "sigrtil1", "sigrtil2", "ctil", "pd", "o_par1"
  )

  expect_type(result, "double")
  expect_equal(names(result), expected_names)
  expect_true(all(is.na(result)))
})


test_that("create_mask fills provided mask correctly", {
  mask <- c(mu1 = 11, mu2 = 5, pd = 0.8, ctil = -2)
  result <- create_mask(mask = mask, p = 2)

  # Use [[ to drop names
  expect_equal(result[["mu1"]], 11)
  expect_equal(result[["mu2"]], 5)
  expect_equal(result[["pd"]], 0.8)
  expect_equal(result[["ctil"]], -2)

  # Others remain NA
  expect_true(is.na(result[["sigltil1"]]))
  expect_true(is.na(result[["o_par1"]]))
})


test_that("create_mask handles full parameter set", {
  full_params <- c(
    mu1 = 10, mu2 = 20,
    sigltil1 = 1, sigltil2 = 2,
    sigrtil1 = 3, sigrtil2 = 4,
    o_par1 = 0.5, ctil = -1, pd = 0.9
  )

  result <- create_mask(mask = full_params, p = 2)

  expect_equal(result[names(full_params)], full_params)
})
