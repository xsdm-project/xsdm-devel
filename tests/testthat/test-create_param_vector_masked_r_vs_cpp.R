# Parity test: pure-R create_param_vector_masked_r vs canonical
# create_param_vector_masked (C++ wrapper).
# Tolerance is 1e-6 (ecological precision); typical observed diff is at
# the level of floating-point round-off (~1e-15).

tol <- 1e-6

canonical_for_p <- function(p) names(make_mask_names(p))

test_that("create_param_vector_masked: example_1 fixture parity (p = 2)", {
  full_names <- canonical_for_p(2L)
  vals <- as.numeric(example_1$par_vec)
  names(vals) <- full_names
  pv   <- vals[setdiff(full_names, c("mu1", "pd"))]
  mask <- vals[c("mu1", "pd")]
  cpp <- create_param_vector_masked(pv, mask = mask, p = 2L)
  r   <- xsdm:::create_param_vector_masked_r(pv, mask = mask, p = 2L)
  expect_equal(unname(cpp), unname(r), tolerance = tol)
  expect_equal(names(cpp), names(r))
})

test_that("create_param_vector_masked: random sweep p in {1, 2, 3}", {
  for (p in 1:3) {
    full_names <- canonical_for_p(p)
    for (seed in 1:5) {
      set.seed(seed)
      vals <- stats::runif(length(full_names), -1, 1)
      names(vals) <- full_names
      # Random partition of names into mask vs param_vector.
      n <- length(full_names)
      k <- max(1L, n %/% 3L)
      mask_idx <- sample.int(n, k)
      mask <- vals[mask_idx]
      pv   <- vals[-mask_idx]
      cpp <- create_param_vector_masked(pv, mask = mask, p = p)
      r   <- xsdm:::create_param_vector_masked_r(pv, mask = mask, p = p)
      expect_equal(unname(cpp), unname(r), tolerance = tol)
      expect_equal(names(cpp), names(r))
    }
  }
})

test_that("create_param_vector_masked: NULL mask path", {
  p <- 2L
  full_names <- canonical_for_p(p)
  set.seed(101)
  vals <- stats::runif(length(full_names), -1, 1)
  names(vals) <- full_names
  cpp <- create_param_vector_masked(vals, mask = NULL, p = p)
  r   <- xsdm:::create_param_vector_masked_r(vals, mask = NULL, p = p)
  expect_equal(unname(cpp), unname(r), tolerance = tol)
  expect_equal(names(cpp), names(r))
})

test_that("create_param_vector_masked: rejects overlapping names", {
  pv   <- c(mu1 = 1.0, sigltil1 = 0.5, sigrtil1 = 0.4, ctil = 0.1, pd = 0.0)
  mask <- c(mu1 = 2.0)
  expect_error(create_param_vector_masked(pv, mask = mask, p = 1L))
  expect_error(xsdm:::create_param_vector_masked_r(pv, mask = mask, p = 1L))
})

test_that("create_param_vector_masked: rejects unknown names", {
  pv   <- c(zzz = 1.0)
  expect_error(create_param_vector_masked(pv, mask = NULL, p = 1L))
  expect_error(xsdm:::create_param_vector_masked_r(pv, mask = NULL, p = 1L))
})

test_that("create_param_vector_masked: rejects missing slots", {
  pv <- c(mu1 = 0.0)  # missing sigltil1, sigrtil1, ctil, pd
  expect_error(create_param_vector_masked(pv, mask = NULL, p = 1L))
  expect_error(xsdm:::create_param_vector_masked_r(pv, mask = NULL, p = 1L))
})
