library(testthat)

test_that("xsdm:::get_range_df_ returns correct structure and values", {
  # Create a small environmental array: p = 2 env_var, 3 time steps, 2 locations
  set.seed(123)
  env_dat <- array(runif(2 * 3 * 2, min = 0.1, max = 1), dim = c(2, 3, 2))

  # Call the function
  result <- xsdm:::get_range_df_(env_dat)

  # Check class and dimensions
  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 3) # lower, center, upper
  expect_true(all(c("lower", "center", "upper") %in% colnames(result)))

  # Check row names include mu, sigl, sigr, ctil, pd, o_par
  rn <- rownames(result)
  expect_true(any(grepl("^mu", rn)))
  expect_true(any(grepl("^sigl", rn)))
  expect_true(any(grepl("^sigr", rn)))
  expect_true("ctil" %in% rn)
  expect_true("pd" %in% rn)
  expect_true(any(grepl("^o_par", rn)))

  # check that the "lower" column is always less than the "center" column, and the
  # "center" column is less than the "upper" column
  expect_true(all(result$lower <= result$center))
  expect_true(all(result$center <= result$upper))

  # Check that no NA remains in rows for mu, sigl, sigr
  mu_rows <- grep("^mu", rn)
  expect_false(any(is.na(result[mu_rows, ])))

  
})

test_that("breadth = 1 reproduces legacy quant_vec = c(0.1, 0.5, 0.9)", {
  set.seed(1)
  env <- array(rnorm(10 * 5 * 2), dim = c(10, 5, 2))
  legacy <- suppressWarnings(
    xsdm:::get_range_df_(env, quant_vec = c(0.1, 0.5, 0.9))
  )
  new <- xsdm:::get_range_df_(env)              # default breadth = 1
  expect_equal(new, legacy, tolerance = 1e-12)
})

test_that("breadth = 0 yields a near-point range (all three columns nearly equal)", {
  env <- array(rnorm(10 * 5 * 2), dim = c(10, 5, 2))
  r <- xsdm:::get_range_df_(env, breadth = 0)
  expect_true(all(abs(r$upper[-9] - r$lower[-9]) < 1e-5))
})

test_that("breadth is monotonic in width", {
  env <- array(rnorm(10 * 5 * 2), dim = c(10, 5, 2))
  r0  <- xsdm:::get_range_df_(env, breadth = 0.1)
  r1  <- xsdm:::get_range_df_(env, breadth = 0.5)
  r2  <- xsdm:::get_range_df_(env, breadth = 1.0)
  expect_true(all((r1$upper - r1$lower) - (r0$upper - r0$lower) >= -1e-12))
  expect_true(all((r2$upper - r2$lower) - (r1$upper - r1$lower) >= -1e-12))
})

test_that("breadth bounds are enforced", {
  env <- array(rnorm(10 * 5 * 2), dim = c(10, 5, 2))
  expect_error(xsdm:::get_range_df_(env, breadth = -0.01))
  expect_error(xsdm:::get_range_df_(env, breadth = 1.01))
})

test_that("quant_vec still works but warns", {
  env <- array(rnorm(10 * 5 * 2), dim = c(10, 5, 2))
  expect_warning(xsdm:::get_range_df_(env, quant_vec = c(0.2, 0.5, 0.8)),
                 "deprecated")
})

