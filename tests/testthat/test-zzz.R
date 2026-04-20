library(testthat)
test_that(".onLoad sets future.globals.maxSize correctly", {
  # Save current option to restore later
  old_option <- getOption("future.globals.maxSize")

  # Call .onLoad
  .onLoad(libname = NULL, pkgname = NULL)

  # Check if the option was set to 2 GB
  expect_equal(getOption("future.globals.maxSize"), 8.0 * 1024^3)

  # Restore old option
  options(future.globals.maxSize = old_option)
})
