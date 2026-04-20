testthat::test_that("build_orthogonal_matrix accepts only numeric values", {
  testthat::expect_error(build_orthogonal_matrix("a"), "numeric")
})
testthat::test_that(
  "build_orthogonal_matrix accepts NULL value and return a matrix with
                    entry 1 a dimenssion 1",
  {
    testthat::expect_equal(
      build_orthogonal_matrix(NULL),
      matrix(1, 1, 1)
    )
  }
)


test_that("build_orthogonal_matrix builds correct matrix for valid entries", {
  # Choose entries so that k is an integer (f(n) = 0.5 * (1 + sqrt(8n + 1)))
  # For n = 3, k = 3
  entries <- c(0.1, 0.2, 0.3)

  result <- build_orthogonal_matrix(entries)

  # Check dimensions
  expect_equal(dim(result), c(3, 3))

  # Check matrix properties: expm of skew-symmetric matrix should be orthogonal
  expect_equal(round(t(result) %*% result, 6), diag(3)) # Orthogonality check

  # Ensure no NA values
  expect_false(any(is.na(result)))
})
