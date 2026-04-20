library(testthat)


test_that("tests for the function convert_equivalence_class", {
  skip_on_cran()
  # set up a random example
  S3 <- matrix(rnorm(9), 3, 3)
  S3 <- S3 - t(S3)
  o_mat <- expm::expm(S3)
  sigltil <- c(1, 2, 3)
  sigrtil <- c(4, 5, 6)
  p <- list(o_mat = o_mat, sigltil = sigltil, sigrtil = sigrtil)
  res <- convert_equivalence_class(p, c(1, 0, 0), c(1, 2, 3))

  # test for correct format of output
  expect_equal(names(res), c("o_mat", "sigltil", "sigrtil"))
  expect_equal(class(res$o_mat), c("matrix", "array"))
  expect_equal(dim(res$o_mat), c(3, 3))
  expect_equal(class(res$sigltil), "numeric")
  expect_equal(length(res$sigltil), 3)
  expect_equal(class(res$sigrtil), "numeric")
  expect_equal(length(res$sigrtil), 3)

  # test for accuracy
  h <- o_mat
  h[, 1] <- -h[, 1]
  expect_equal(res$o_mat, h)
  expect_equal(res$sigltil[1], sigrtil[1])
  expect_equal(res$sigrtil[1], sigltil[1])

  # another case, accuracy
  res <- convert_equivalence_class(p, c(0, 0, 0), c(2, 3, 1))
  expect_equal(res$o_mat, o_mat[, c(2, 3, 1)])
  expect_equal(res$sigltil, sigltil[c(2, 3, 1)])
  expect_equal(res$sigrtil, sigrtil[c(2, 3, 1)])
})
