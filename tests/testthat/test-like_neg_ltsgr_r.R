library(testthat)
test_that("like_neg_ltsgr_r test for correct format of output", {
  # an ordinary case
  o_mat <- matrix(c(-0.4443546, 0.8958510, -0.8958510, -0.4443546), ncol = 2)
  mu <- c(11.433373, 5.046939)
  sigltil <- c(1.036834, 1.556083)
  sigrtil <- c(1.538972, 1.458738)
  env_dat <- example_1$env_array
  M <- xsdm:::like_neg_ltsgr_r(env_dat, mu, sigltil, sigrtil, o_mat)
  testthat::expect_equal(class(M), "numeric")
  testthat::expect_equal(length(M), 2000)
  testthat::expect_true(all(is.finite(M)))

  # some cases with Inf entries for sigltil or sigrtil
  sigltil[1] <- Inf
  M <- xsdm:::like_neg_ltsgr_r(env_dat, mu, sigltil, sigrtil, o_mat)
  testthat::expect_equal(class(M), "numeric")
  testthat::expect_equal(length(M), 2000)
  testthat::expect_true(all(is.finite(M)))

  sigrtil[2] <- Inf
  M <- xsdm:::like_neg_ltsgr_r(env_dat, mu, sigltil, sigrtil, o_mat)
  testthat::expect_equal(class(M), "numeric")
  testthat::expect_equal(length(M), 2000)
  testthat::expect_true(all(is.finite(M)))

  sigltil <- c(Inf, 1.556083)
  sigrtil <- c(Inf, 1.458738)
  M <- xsdm:::like_neg_ltsgr_r(env_dat, mu, sigltil, sigrtil, o_mat)
  testthat::expect_equal(class(M), "numeric")
  testthat::expect_equal(length(M), 2000)
  testthat::expect_true(all(is.finite(M)))
})

test_that("like_neg_ltsgr_r test for correctness of output in simple cases", {
  # trivial case where the environment is always optimal
  o_mat <- matrix(c(-0.4443546, 0.8958510, -0.8958510, -0.4443546), ncol = 2)
  mu <- c(11.433373, 5.046939)
  sigltil <- c(1.036834, 1.556083)
  sigrtil <- c(1.538972, 1.458738)
  env_dat <- array(rep(mu, each = 3 * 30), c(3, 30, 2))
  M <- xsdm:::like_neg_ltsgr_r(env_dat, mu, sigltil, sigrtil, o_mat)
  testthat::expect_equal(class(M), "numeric")
  testthat::expect_equal(length(M), 3)
  testthat::expect_equal(M, rep(0, 3))

  # trivial case where the species is insensitive to the environment
  sigltil <- rep(Inf, 2)
  sigrtil <- rep(Inf, 2)
  env_dat <- array(rnorm(3 * 30 * 2), c(3, 30, 2))
  M <- xsdm:::like_neg_ltsgr_r(env_dat, mu, sigltil, sigrtil, o_mat)
  testthat::expect_equal(class(M), "numeric")
  testthat::expect_equal(length(M), 3)
  testthat::expect_equal(M, rep(0, 3))

  # 1d case where sigltil and sigrtil are both 1
  mu <- 0
  sigltil <- 1
  sigrtil <- 1
  env_dat <- array(rnorm(10 * 30), c(10, 30, 1))
  o_mat <- matrix(1, 1, 1)
  M <- xsdm:::like_neg_ltsgr_r(env_dat, mu, sigltil, sigrtil, o_mat)
  testthat::expect_equal(M, 0.5 * apply(FUN = mean, X = env_dat^2, MARGIN = 1))

  # similar but 2d
  mu <- rep(0, 2)
  sigltil <- rep(1, 2)
  sigrtil <- rep(1, 2)
  o_mat <- diag(2)
  env_dat <- array(rnorm(10 * 30 * 2), c(10, 30, 2))
  M <- xsdm:::like_neg_ltsgr_r(env_dat, mu, sigltil, sigrtil, o_mat)
  testthat::expect_equal(
    M,
    0.5 * apply(X = apply(FUN = sum, MARGIN = c(1, 2), X = env_dat^2), MARGIN = 1, FUN = mean)
  )

  # similar but o_mat no longer the identity
  o_mat <- matrix(c(-0.4443546, 0.8958510, -0.8958510, -0.4443546), ncol = 2)
  M <- xsdm:::like_neg_ltsgr_r(env_dat, mu, sigltil, sigrtil, o_mat)
  env_dat_mat <- matrix(aperm(env_dat, c(3, 2, 1)), 2, 30 * 10)
  testthat::expect_equal(
    M,
    0.5 * apply(
      FUN = mean, MARGIN = 2,
      matrix(apply(FUN = sum, MARGIN = 2, (t(o_mat) %*% env_dat_mat)^2), 30, 10)
    )
  )
})
