library(testthat)
library(xsdm)

# ---------------------------------------------------------------------------
# Validates the num_starts >= 3 contract enforced by optimize_likelihood(),
# start_parms(), and get_start_parms_().
#
# Below 3, sobol::sobol_design either segfaults (nseq = 0) or returns a bare
# numeric vector (nseq = 1). Catching this at the public API boundary
# produces an actionable checkmate error instead of a downstream crash.
# ---------------------------------------------------------------------------

test_that("optimize_likelihood() rejects num_starts < 3", {
  set.seed(1)
  env_dat <- array(stats::runif(10 * 5 * 2), dim = c(10, 5, 2))
  occ <- rep(c(1L, 0L), length.out = 10)

  for (bad in c(0L, 1L, 2L)) {
    expect_error(
      optimize_likelihood(env_dat, occ, num_starts = bad, num_threads = 1L),
      regexp = "num_starts",
      info   = sprintf("num_starts = %d should fail validation", bad)
    )
  }
})

test_that("start_parms() rejects num_starts < 3", {
  set.seed(2)
  env_dat <- array(stats::runif(8 * 4 * 2), dim = c(8, 4, 2))

  for (bad in c(0L, 1L, 2L)) {
    expect_error(
      start_parms(env_dat, num_starts = bad),
      regexp = "num_starts",
      info   = sprintf("num_starts = %d should fail validation", bad)
    )
  }
})

test_that("start_parms() accepts num_starts = 3 (the boundary)", {
  set.seed(3)
  env_dat <- array(stats::runif(8 * 4 * 2), dim = c(8, 4, 2))
  out <- start_parms(env_dat, num_starts = 3L)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 3L)
})

test_that("get_start_parms_() rejects numstarts < 3", {
  ranges <- data.frame(
    lower  = c(-1, -2, -3),
    center = c( 0,  0,  0),
    upper  = c( 1,  2,  3),
    row.names = c("a", "b", "c")
  )
  for (bad in c(0L, 1L, 2L)) {
    expect_error(
      xsdm:::get_start_parms_(ranges, numstarts = bad),
      regexp = "numstarts",
      info   = sprintf("numstarts = %d should fail validation", bad)
    )
  }
})
