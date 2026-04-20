library(testthat)

# If running outside the package test environment, optionally stub create_mask/make_mask_names/num_par.
# In your package, remove the stubs below.

test_that("p must be a positive count", {
  expect_error(create_param_vector_masked(param_vector = c(mu1 = 0), p = 0))
  expect_error(create_param_vector_masked(param_vector = c(mu1 = 0), p = -1))
  expect_error(create_param_vector_masked(param_vector = c(mu1 = 0), p = NA_integer_))
})

test_that("param_vector is required, named, numeric, and contains no NAs", {
  expect_error(
    create_param_vector_masked(param_vector = NULL, p = 1),
    "`param_vector` must not be NULL"
  )
  # Unnamed vector
  expect_error(
    create_param_vector_masked(param_vector = c(1, 2), p = 1),
    "Assertion on 'param_vector' failed: Must have names"
  )
  # NA inside param_vector
  expect_error(create_param_vector_masked(param_vector = c(mu1 = NA_real_), p = 1),
    "Assertion on 'param_vector' failed: Contains missing values (element 1).",
    fixed = TRUE
  )
})

test_that("errors on unexpected names in param_vector", {
  expect_error(
    create_param_vector_masked(param_vector = c(foo = 1), p = 1),
    "Unexpected name\\(s\\) in `param_vector`"
  )
})

test_that("mask optional; when present must use canonical names and no NAs", {
  p <- 2
  canon <- names(create_mask(NULL, p))
  good_mask <- setNames(c(0.1, 0.1, 0.05), c("mu1", "mu2", "pd"))
  out <- create_param_vector_masked(
    param_vector = setNames(
      c(1.0, 1.0, 2.0, 2.0, 0.3, 0.0),
      c(
        "sigltil1",
        "sigltil2",
        "sigrtil1",
        "sigrtil2",
        "ctil",
        "o_par1"
      )
    ),
    mask = good_mask,
    p = p
  )
  expect_true(all(names(out) == canon))
  expect_false(anyNA(out))

  bad_mask <- c(extra = 99)
  expect_error(
    create_param_vector_masked(
      param_vector = setNames(
        c(1.0, 2.0, 0.3, 0.0),
        c("sigltil1", "sigrtil1", "ctil", "o_par1")
      ),
      mask = bad_mask, p = p
    ),
    "Unexpected parameter name\\(s\\)"
  )

  mask_with_na <- c(mu1 = NA_real_)
  expect_error(
    create_param_vector_masked(
      param_vector = setNames(
        c(1.0, 2.0, 0.3, 0.0),
        c("sigltil1", "sigrtil1", "ctil", "o_par1")
      ),
      mask = mask_with_na, p = p
    )
  )
})

test_that("final output has canonical names, canonical length, and no NAs (p = 1)", {
  p <- 1
  canon <- create_mask(NULL, p)
  # Names likely: mu1, sigltil1, sigrtil1, ctil, pd (length 5)
  mask <- c(mu1 = -1, pd = 0.5)
  pv <- c(sigltil1 = 1.0, sigrtil1 = 2.0, ctil = 0.2)
  out <- create_param_vector_masked(param_vector = pv, mask = mask, p = p)

  expect_identical(names(out), names(canon))
  expect_equal(length(out), length(canon))
  expect_false(anyNA(out))
  # Precedence: param_vector overrides mask where overlapping
  expect_equal(as.numeric(out["mu1"]), -1) # from mask (no override in pv)
  expect_equal(as.numeric(out["ctil"]), 0.2) # from pv
})

test_that("p = 2: must include o_par1; precedence and fill checks", {
  p <- 2
  canon <- create_mask(NULL, p)
  # Fill all entries: mask first, then param overrides
  mask <- c(mu1 = 0.1, mu2 = 0.2, pd = 0.05)
  pv <- c(
    sigltil1 = 1.0, sigltil2 = 1.1,
    sigrtil1 = 2.0, sigrtil2 = 2.2,
    ctil = 0.3, o_par1 = 0.0
  )
  out <- create_param_vector_masked(param_vector = pv, mask = mask, p = p)

  expect_identical(names(out), names(canon))
  expect_equal(length(out), length(canon))
  expect_false(anyNA(out))
  expect_equal(as.numeric(out["mu1"]), 0.1)
  expect_equal(as.numeric(out["o_par1"]), 0.0)

  # If we omit o_par1, should error (NA remains)
  pv_incomplete <- c(sigltil1 = 1, sigltil2 = 1.1, sigrtil1 = 2, sigrtil2 = 2.2, ctil = 0.3)
  expect_error(
    create_param_vector_masked(param_vector = pv_incomplete, mask = mask, p = p),
    "contains NA values"
  )
})

test_that("p = 3: must include o_par1..o_par3; no NA remains", {
  p <- 3
  canon <- create_mask(NULL, p)

  mask <- c(mu1 = 0.1, mu2 = 0.2, mu3 = 0.3, pd = 0.01)
  pv <- c(
    sigltil1 = 1.0, sigltil2 = 1.1, sigltil3 = 1.2,
    sigrtil1 = 2.0, sigrtil2 = 2.1, sigrtil3 = 2.2,
    ctil = 0.4, o_par1 = -0.2, o_par2 = 0.0, o_par3 = 0.15
  )
  out <- create_param_vector_masked(param_vector = pv, mask = mask, p = p)

  expect_identical(names(out), names(canon))
  expect_equal(length(out), length(canon))
  expect_false(anyNA(out))
})


test_that("complementary names succeed and produce full, non-NA, canonical output", {
  p <- 3
  canon_names <- names(create_mask(NULL, p))

  # No overlap between mask and param_vector
  mask <- c(mu1 = 0.1, mu2 = 0.2, mu3 = 0.3, pd = 0.01)
  pv <- c(
    sigltil1 = 1.0, sigltil2 = 1.1, sigltil3 = 1.2,
    sigrtil1 = 2.0, sigrtil2 = 2.1, sigrtil3 = 2.2,
    ctil = 0.4, o_par1 = -0.2, o_par2 = 0.0, o_par3 = 0.15
  )

  out <- create_param_vector_masked(param_vector = pv, mask = mask, p = p)

  expect_identical(names(out), canon_names)
  expect_equal(length(out), length(canon_names))
  expect_false(anyNA(out))

  # Check a few values from each source
  expect_equal(as.numeric(out["mu1"]), 0.1) # from mask
  expect_equal(as.numeric(out["ctil"]), 0.4) # from param_vector
  expect_equal(as.numeric(out["o_par3"]), 0.15) # from param_vector
})
test_that("errors list which names are missing when NA remains", {
  p <- 1
  # Provide only mu1; others should be missing -> error lists names
  expect_error(
    create_param_vector_masked(param_vector = c(mu1 = 3), mask = NULL, p = p),
    regexp = "The final parameter vector contains NA values for:"
  )
})


test_that("param_vector and mask names must be complementary (no overlap)", {
  p <- 3

  # Overlapping case -> should error listing the overlapping names
  mask_overlap <- c(mu1 = 0.1, mu2 = 0.2, mu3 = 0.3, pd = 0.01)
  pv_overlap <- c(
    mu1 = 2, mu2 = 0.2, # overlaps with mask -> must fail
    sigltil1 = 1.0, sigltil2 = 1.1, sigltil3 = 1.2,
    sigrtil1 = 2.0, sigrtil2 = 2.1, sigrtil3 = 2.2,
    ctil = 0.4, o_par1 = -0.2, o_par2 = 0.0, o_par3 = 0.15
  )

  expect_error(
    create_param_vector_masked(param_vector = pv_overlap, mask = mask_overlap, p = p),
    regexp = "must be complementary \\(disjoint\\).*Overlapping names:.*mu1.*mu2"
  )

  # Disjoint case -> should succeed and produce full, non-NA, canonical output
  mask_disjoint <- c(mu1 = 0.1, mu2 = 0.2, mu3 = 0.3, pd = 0.01)
  pv_disjoint <- c(
    sigltil1 = 1.0, sigltil2 = 1.1, sigltil3 = 1.2,
    sigrtil1 = 2.0, sigrtil2 = 2.1, sigrtil3 = 2.2,
    ctil = 0.4, o_par1 = -0.2, o_par2 = 0.0, o_par3 = 0.15
  )

  out <- create_param_vector_masked(param_vector = pv_disjoint, mask = mask_disjoint, p = p)
  canon_names <- names(create_mask(mask = NULL, p = p))
  expect_identical(names(out), canon_names)
  expect_equal(length(out), length(canon_names))
  expect_false(anyNA(out))

  # spot-check values from each source to confirm application order
  expect_equal(as.numeric(out["mu1"]), 0.1) # from mask
  expect_equal(as.numeric(out["ctil"]), 0.4) # from param_vector
  expect_equal(as.numeric(out["o_par3"]), 0.15) # from param_vector
})
