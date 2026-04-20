library(testthat)

test_that("get_start_parms returns correct tibble with expected structure", {
  # Mock ranges data frame
  ranges <- data.frame(
    lower = c(-1, 0.1),
    center = c(0, 0.5),
    upper = c(1, 1.0),
    row.names = c("param1", "param2")
  )

  # Call the function
  result <- get_start_parms(ranges, numstarts = 10)

  # Check output type
  expect_s3_class(result, "tbl_df")

  # Check number of rows: numstarts
  expect_equal(nrow(result), 10)

  # Check column names match row names of ranges
  expect_true(all(colnames(result) %in% rownames(ranges)))

  # Check values are within bounds
  for (param in colnames(result)) {
    expect_true(all(result[[param]] >= ranges[param, "lower"]))
    expect_true(all(result[[param]] <= ranges[param, "upper"]))
  }
})

test_that("get_start_parms fails with invalid inputs", {
  # Wrong column names
  bad_ranges <- data.frame(a = 1, b = 2, c = 3)
  expect_error(get_start_parms(bad_ranges, numstarts = 10))

  # Wrong type for numstarts
  expect_error(get_start_parms(ranges, numstarts = "ten"))
})
