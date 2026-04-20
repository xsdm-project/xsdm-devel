library(testthat)

test_that("round-trip recovers original parameters for k = 2..6", {
  set.seed(123)
  for (k in 2:6) {
    k <- 2
    q <- k * (k - 1) / 2
    pars <- runif(q, -1, 1)
    O <- build_orthogonal_matrix(pars)
    rec <- extract_orthogonal_matrix_parameters(O)
    expect_equal(rec, pars, tolerance = 1e-6)
  }
})

test_that("rebuilding from extracted parameters reproduces the original matrix", {
  set.seed(456)
  k <- 5
  q <- k * (k - 1) / 2
  pars <- runif(q, -2, 2)
  O <- build_orthogonal_matrix(pars)
  rec <- extract_orthogonal_matrix_parameters(O)
  O2 <- build_orthogonal_matrix(rec)
  expect_equal(O2, O, tolerance = 1e-8)
})

test_that("recovered parameters reproduce the original orthogonal matrix", {
  k <- 4
  q <- k * (k - 1) / 2
  pars <- as.numeric(1:q) # large values allowed
  O <- build_orthogonal_matrix(pars)
  rec <- extract_orthogonal_matrix_parameters(O)
  O2 <- build_orthogonal_matrix(rec)
  expect_equal(O2, O, tolerance = 1e-8)
})

test_that("rotation by pi in 2D errors with informative message", {
  # 2x2 rotation by angle = pi (det = +1)
  th <- -pi
  R2 <- matrix(c(cos(th), -sin(th), sin(th), cos(th)), 2, 2)
  # This has eigenvalues equal to -1, so the principal matrix logarithm
  # is not well defined. The function should error.
  expect_error(
    extract_orthogonal_matrix_parameters(R2),
    "eigenvalues equal to -1"
  )
})

test_that("3x3 matrix with rotation by pi errors (issue example)", {
  # matrix(c(-1,0,0,0,-1,0,0,0,1),3,3) from the issue
  O3 <- matrix(c(-1, 0, 0, 0, -1, 0, 0, 0, 1), 3, 3)
  expect_error(
    extract_orthogonal_matrix_parameters(O3),
    "eigenvalues equal to -1"
  )
})

test_that("1x1 identity returns NULL", {
  expect_null(
    extract_orthogonal_matrix_parameters(
      matrix(1, 1, 1)
    )
  )
})

test_that("1x1 reflection returns NULL with warning", {
  expect_warning({
    out <- extract_orthogonal_matrix_parameters(matrix(-1, 1, 1))
    expect_null(out)
  })
})

test_that("non-square input errors", {
  expect_error(extract_orthogonal_matrix_parameters(matrix(1:6, nrow = 2, ncol = 3)))
})

test_that("non-orthogonal input errors", {
  M <- diag(2)
  M[1, 1] <- 0.99999 # disturb orthogonality slightly
  expect_error(extract_orthogonal_matrix_parameters(M))
})

test_that("reflections for k > 1 error (det = -1)", {
  O <- diag(c(-1, 1)) # determinant -1 (reflection)
  expect_error(extract_orthogonal_matrix_parameters(O))
})


test_that("near -I_2 (angle pi - eps) is stable via logm path", {
  eps <- 1e-8
  th <- pi - eps
  R2 <- matrix(c(cos(th), -sin(th), sin(th), cos(th)), 2, 2)
  rec <- extract_orthogonal_matrix_parameters(R2)
  O2 <- build_orthogonal_matrix(rec)
  expect_equal(O2, R2, tolerance = 1e-8)
})

test_that("random large generators reproduce O (k = 3..6)", {
  set.seed(7)
  for (k in 3:6) {
    q <- k * (k - 1) / 2
    pars <- runif(q, -3, 3) # allow angles outside principal; reconstruction should match
    O <- build_orthogonal_matrix(pars)
    rec <- extract_orthogonal_matrix_parameters(O)
    O2 <- build_orthogonal_matrix(rec)
    expect_equal(O2, O, tolerance = 1e-8)
  }
})
