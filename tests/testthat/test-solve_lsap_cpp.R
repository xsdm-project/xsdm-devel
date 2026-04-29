library(testthat)

# ---------------------------------------------------------------------------
# Parity tests for the in-package C++ Hungarian solver `.solve_lsap_cpp`.
#
# `R/dist_between_params.R` calls `.solve_lsap_cpp` to solve the linear sum
# assignment problem on the cost matrix it builds in R. The reference R-level
# implementation that we must match is `clue::solve_LSAP`. These tests cover:
#
#   * trivial 1x1 / 2x2 hand-rolled matrices
#   * a known small example with a unique optimum
#   * random non-negative cost matrices of several sizes
#
# All tests are gated on `clue` being installed (it sits in Suggests:).
# ---------------------------------------------------------------------------

test_that(".solve_lsap_cpp matches clue::solve_LSAP on tiny hand-rolled cases", {
  skip_if_not_installed("clue")

  # 1x1: trivially the identity assignment.
  expect_equal(as.integer(xsdm:::.solve_lsap_cpp(matrix(7.5, 1, 1))),
               as.integer(clue::solve_LSAP(matrix(7.5, 1, 1))))

  # 2x2 with a unique optimum (off-diagonal pairing).
  cost_2x2 <- matrix(c(10, 1,
                       1, 10), 2, 2, byrow = TRUE)
  expect_equal(as.integer(xsdm:::.solve_lsap_cpp(cost_2x2)),
               as.integer(clue::solve_LSAP(cost_2x2)))

  # 3x3 small example with a unique optimum (assigns row i to column i+1 cyc.).
  cost_3x3 <- matrix(c(8, 1, 9,
                       9, 8, 1,
                       1, 9, 8), 3, 3, byrow = TRUE)
  expect_equal(as.integer(xsdm:::.solve_lsap_cpp(cost_3x3)),
               as.integer(clue::solve_LSAP(cost_3x3)))
})

test_that(".solve_lsap_cpp matches clue::solve_LSAP on random cost matrices", {
  skip_if_not_installed("clue")

  set.seed(20240429)
  for (n in c(2L, 3L, 5L, 8L, 12L)) {
    for (rep in seq_len(5L)) {
      cost <- matrix(runif(n * n, min = 0, max = 100), n, n)
      perm_cpp  <- as.integer(xsdm:::.solve_lsap_cpp(cost))
      perm_clue <- as.integer(clue::solve_LSAP(cost))

      # Optima may not be unique; instead of comparing permutations, compare
      # the achieved cost (this is what dist_between_params actually uses).
      pairing_cpp  <- cbind(seq_len(n), perm_cpp)
      pairing_clue <- cbind(seq_len(n), perm_clue)
      expect_equal(sum(cost[pairing_cpp]),
                   sum(cost[pairing_clue]),
                   tolerance = 1e-10,
                   info = sprintf("n=%d, rep=%d", n, rep))
    }
  }
})

test_that(".solve_lsap_cpp rejects non-square / non-finite inputs", {
  expect_error(xsdm:::.solve_lsap_cpp(matrix(0, 2, 3)),
               "square")
  expect_error(xsdm:::.solve_lsap_cpp(matrix(c(1, NA_real_, 3, 4), 2, 2)),
               "finite")
  expect_error(xsdm:::.solve_lsap_cpp(matrix(c(1, Inf, 3, 4), 2, 2)),
               "finite")
})
