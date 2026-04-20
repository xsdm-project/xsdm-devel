library(testthat)

test_that("tests for the function distance_between_params", {
  # generate some example inputs
  ang <- -46 * pi / 180
  o_mat1 <- matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), 2, 2)
  sigltil1 <- c(1, 2)
  sigrtil1 <- c(3, 4)
  mu <- c(0, 0)
  ctil <- 10
  pd <- .98
  p1 <- list(mu = mu, ctil = ctil, pd = pd, o_mat = o_mat1, sigltil = sigltil1, sigrtil = sigrtil1)

  ang <- -44 * pi / 180
  o_mat2 <- matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), 2, 2)
  sigltil2 <- c(1, 2)
  sigrtil2 <- c(3, 4)
  p2 <- list(mu = mu, ctil = ctil, pd = pd, o_mat = o_mat2, sigltil = sigltil2, sigrtil = sigrtil2)

  # check format
  res1 <- distance_between_params(p1, p2, FALSE)
  res2 <- distance_between_params(p1, p2, TRUE)
  expect_equal(class(res1), "numeric")
  expect_equal(length(res1), 1)
  expect_equal(class(res2), "list")
  expect_equal(length(res2), 2)
  expect_equal(names(res2), c("distance", "representative"))
  expect_equal(res1, res2$distance)
  expect_equal(class(res2$representative), "list")
  expect_equal(length(res2$representative), 6)
  expect_equal(names(res2$representative), c("mu", "sigltil", "sigrtil", "ctil", "pd", "o_mat"))
  expect_equal(res2$representative$mu, mu)
  expect_equal(res2$representative$ctil, ctil)
  expect_equal(res2$representative$pd, pd)

  # now test accuracy
  expect_equal(res1, sqrt(sum((o_mat1 - o_mat2)^2)))

  p2$mu <- c(1, 1)
  p2$ctil <- 3
  p2$pd <- .1
  res2 <- distance_between_params(p1, p2, TRUE)
  expect_equal(res2$distance, sqrt(res1^2 + 1 + 1 + (10 - 3)^2 + (.98 - .1)^2))

  # ANGEL to consider adding some more tests
})
