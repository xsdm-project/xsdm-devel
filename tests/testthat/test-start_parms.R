# tests/testthat/test-start_parms.R
library(testthat)
test_that("start_parms() returns expected shape and names with no mask", {
  # Subset environmental array to presences, as in the examples
  env_dat_pres <- example_1_env_array[example_1_occurrence_vector == 1, , ]

  # Small number of starts for speed/determinism
  set.seed(123)
  num_starts <- 7L
  sp <- start_parms(
    env_dat = env_dat_pres,
    mask = NULL,
    breadth = 1,
    num_starts = num_starts
  )

  # Basic structure
  expect_s3_class(sp, "tbl_df")
  expect_equal(nrow(sp), num_starts)

  # Column names should match all "free" parameters (i.e., all rows in get_range_df)
  range_df <- get_range_df(env_dat_pres, 1)
  expect_setequal(colnames(sp), rownames(range_df))

  # All values must lie within [lower, upper] of range_df
  for (nm in colnames(sp)) {
    lower <- range_df[nm, "lower", drop = TRUE]
    upper <- range_df[nm, "upper", drop = TRUE]
    expect_true(all(sp[[nm]] >= lower & sp[[nm]] <= upper),
      info = paste("Parameter", nm, "violates bounds")
    )
  }
})

test_that("start_parms() respects mask by excluding masked parameters", {
  env_dat_pres <- example_1_env_array[example_1_occurrence_vector == 1, , ]

  # For p = 2, typical canonical names include:
  # mu1, mu2, sigltil1, sigltil2, sigrtil1, sigrtil2, ctil, pd, o_par1
  mask <- c(mu2 = 5, pd = 1) # fixed values => should be excluded from output

  set.seed(456)
  out <- start_parms(
    env_dat = env_dat_pres,
    mask = mask,
    breadth = 1,
    num_starts = 5
  )

  expect_false("mu2" %in% names(out))
  expect_false("pd" %in% names(out))

  # Everything else should still be present
  range_df <- get_range_df(env_dat = env_dat_pres, breadth = 1)
  expected_names <- setdiff(rownames(range_df), names(mask))
  expect_setequal(names(out), expected_names)
})

test_that("start_parms() is (roughly) deterministic with a fixed seed", {
  env_dat_pres <- example_1_env_array[example_1_occurrence_vector == 1, , ]

  set.seed(42)
  a <- start_parms(env_dat_pres, num_starts = 5)

  set.seed(42)
  b <- start_parms(env_dat_pres, num_starts = 5)

  # Because Sobol' design is deterministic under a fixed seed + inputs,
  # we expect identical frames.
  expect_identical(a, b)
})

test_that("start_parms() input validation errors", {
  # wrong env_dat type
  expect_error(
    start_parms(env_dat = NULL),
    regexp = "Must be of type 'array'", # from checkmate::assert_array
    fixed = FALSE
  )

  # wrong quant_vec length
  env_dat_pres <- example_1_env_array[example_1_occurrence_vector == 1, , ]
  expect_error(
    start_parms(env_dat_pres, breadth  = c(0.1, 0.9)),
    regexp = "Must have length 1", # from checkmate::assert_numeric(len = 3)
    fixed = FALSE
  )

  # non-scalar num_starts
  expect_error(
    start_parms(env_dat_pres, num_starts = c(5, 10)),
    regexp = "Must have length 1", # from checkmate::assert_number()
    fixed = FALSE
  )
})
