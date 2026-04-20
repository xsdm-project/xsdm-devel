library(testthat)

test_that("make_mask_names returns a named numeric vector", {
  res <- make_mask_names(2)
  expect_type(res, "double")
  expect_true(!is.null(names(res)))
})

test_that("names order is correct for p = 1 (no o_par)", {
  p <- 1
  res <- make_mask_names(p)
  expected_names <- c("mu1", "sigltil1", "sigrtil1", "ctil", "pd")
  expect_identical(names(res), expected_names)
  expect_false(any(grepl("^o_par", names(res))))
})

test_that("names order is correct for p = 2", {
  p <- 2
  res <- make_mask_names(p)
  expected_names <- c(
    "mu1", "mu2",
    "sigltil1", "sigltil2",
    "sigrtil1", "sigrtil2",
    "ctil", "pd",
    "o_par1"
  )
  expect_identical(names(res), expected_names)
})

test_that("names order is correct for p = 3", {
  p <- 3
  res <- make_mask_names(p)
  expected_names <- c(
    paste0("mu", 1:3),
    paste0("sigltil", 1:3),
    paste0("sigrtil", 1:3),
    "ctil", "pd",
    paste0("o_par", 1:3)
  )
  expect_identical(names(res), expected_names)
})

test_that("all values are NA_real_", {
  for (p in c(1, 2, 4)) {
    res <- make_mask_names(p)
    expect_true(all(is.na(res)))
    expect_true(is.double(res))
  }
})

# ---------- Error cases (robust) ----------

test_that("invalid p raises informative messages", {
  expect_error(make_mask_names(0), "Must be >= 1")
  expect_error(make_mask_names(-1), "Must be >= 1")
  expect_error(make_mask_names(1.5), "Must be of type 'count'")
  expect_error(make_mask_names(NA_integer_), "May not be NA")
})


test_that("internal length check equals num_par(p)", {
  for (p in 1:6) {
    res <- make_mask_names(p)
    expect_equal(length(res), num_par(p))
  }
})
