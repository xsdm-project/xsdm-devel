# tests/testthat/test-expm_skew_cpp.R
#
# Parity test: xsdm:::.build_orthogonal_matrix_cpp must agree with
# build_orthogonal_matrix() to near machine precision.

test_that(".build_orthogonal_matrix_cpp: p = 1 returns 1x1 identity", {
  o_cpp <- xsdm:::.build_orthogonal_matrix_cpp(numeric(0))
  expect_equal(dim(o_cpp), c(1L, 1L))
  expect_equal(o_cpp[1, 1], 1)
})

test_that(".build_orthogonal_matrix_cpp: p = 2 matches closed-form rotation", {
  thetas <- c(0, 0.1, -0.7, 1.5, pi / 3)
  for (th in thetas) {
    o_ref <- build_orthogonal_matrix(th)
    o_cpp <- xsdm:::.build_orthogonal_matrix_cpp(th)
    expect_equal(o_cpp, o_ref, tolerance = 1e-14)
    # orthogonality
    expect_equal(t(o_cpp) %*% o_cpp, diag(2), tolerance = 1e-14)
  }
})

test_that(".build_orthogonal_matrix_cpp: p = 3 matches R expm", {
  set.seed(17)
  for (trial in 1:10) {
    entries <- rnorm(3, sd = 1.0)
    o_ref <- build_orthogonal_matrix(entries)
    o_cpp <- xsdm:::.build_orthogonal_matrix_cpp(entries)
    expect_equal(o_cpp, o_ref, tolerance = 1e-10)
    expect_equal(t(o_cpp) %*% o_cpp, diag(3), tolerance = 1e-10)
  }
})

test_that(".build_orthogonal_matrix_cpp: p = 4 and p = 5 match R expm", {
  set.seed(29)
  for (p in 4:5) {
    q <- p * (p - 1) / 2
    for (trial in 1:5) {
      entries <- rnorm(q, sd = 0.8)
      o_ref <- build_orthogonal_matrix(entries)
      o_cpp <- xsdm:::.build_orthogonal_matrix_cpp(entries)
      expect_equal(o_cpp, o_ref, tolerance = 1e-10)
      expect_equal(t(o_cpp) %*% o_cpp, diag(p), tolerance = 1e-10)
    }
  }
})

test_that(".build_orthogonal_matrix_cpp: large-norm entries still match R expm", {
  # Entries large enough to exercise the scaling-and-squaring branch.
  set.seed(101)
  entries <- rnorm(3, sd = 5.0)   # norm well above theta13 = 5.37
  o_ref <- build_orthogonal_matrix(entries)
  o_cpp <- xsdm:::.build_orthogonal_matrix_cpp(entries)
  expect_equal(o_cpp, o_ref, tolerance = 1e-10)
})

test_that(".build_orthogonal_matrix_cpp: invalid length errors", {
  # Length 2 is not a triangular number
  expect_error(
    xsdm:::.build_orthogonal_matrix_cpp(c(0.1, 0.2)),
    "triangular number"
  )
})
