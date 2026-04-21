library(testthat)

test_that("math_to_bio converts p=1 correctly (no o_par)", {
  # Build canonical vector for p=1
  math_vec <- make_mask_names(1)
  math_vec[] <- c(
    mu1 = 11.0,
    sigltil1 = log(1.2),
    sigrtil1 = log(0.8),
    ctil = -3.0,
    pd = xsdm:::logit(0.85)
  )
  
  bio <- math_to_bio(math_vec)
  
  expect_type(bio, "list")
  expect_named(bio, c("mu", "sigltil", "sigrtil", "ctil", "pd", "o_mat"))
  expect_equal(bio$mu, 11.0)
  expect_equal(bio$sigltil, 1.2)
  expect_equal(bio$sigrtil, 0.8)
  expect_equal(bio$ctil, -3.0)
  expect_equal(bio$pd, 0.85)
  expect_equal(bio$o_mat, matrix(1, 1, 1))
})

test_that("math_to_bio converts p=2 correctly (with o_par1)", {
  math_vec <- examples$optim_par_vec
  bio <- math_to_bio(math_vec)
  
  expect_type(bio, "list")
  expect_named(bio, c("mu", "sigltil", "sigrtil", "ctil", "pd", "o_mat"))
  expect_length(bio$mu, 2)
  expect_length(bio$sigltil, 2)
  expect_length(bio$sigrtil, 2)
  expect_length(bio$ctil, 1)
  expect_length(bio$pd, 1)
  expect_equal(dim(bio$o_mat), c(2, 2))
  expect_true(max(abs(t(bio$o_mat) %*% bio$o_mat - diag(2))) < 1e-12)
})

test_that("math_to_bio round-trips with bio_to_math (p=1,2,3)", {
  set.seed(42)
  for (p in 1:3) {
    nm <- make_mask_names(p)
    math1 <- nm
    math1[] <- rnorm(length(math1), mean = 0, sd = 2)
    o_inds <- grep("^o_par", names(math1))
    if (length(o_inds) > 0) {
      math1[o_inds] <- runif(length(o_inds), -0.5, 0.5)
    }
    
    bio <- math_to_bio(math1)
    math2 <- bio_to_math(bio)
    
    common_names <- intersect(names(math1), names(math2))
    expect_equal(math2[common_names], math1[common_names], tolerance = 1e-8)
  }
})

test_that("math_to_bio handles extreme math-scale values", {
  p <- 1
  math_vec <- make_mask_names(1)
  math_vec[] <- c(mu1 = 0, sigltil1 = log(2), sigrtil1 = log(3), ctil = 0, pd = 100)
  bio <- math_to_bio(math_vec)
  expect_equal(bio$pd, 1.0, tolerance = 1e-10)
  
  math_vec2 <- make_mask_names(1)
  math_vec2[] <- c(mu1 = 0, sigltil1 = 10, sigrtil1 = 10, ctil = 0, pd = 0)
  bio2 <- math_to_bio(math_vec2)
  expect_true(bio2$sigltil > 1000)
  expect_true(bio2$sigrtil > 1000)
})

test_that("math_to_bio errors on malformed inputs", {
  # Unnamed vector
  expect_error(math_to_bio(c(1, 2, 3, 4, 5)), "Must have names")
  
  # Missing required name (sigrtil1)
  bad1 <- c(mu1 = 1, sigltil1 = 0, ctil = 0, pd = 0)  # wrong length and names
  expect_error(math_to_bio(bad1), " Must have length >= 5")
  
  # Correct length but wrong names
  bad2 <- c(mu1 = 1, mu2 = 2, sigltil1 = 0, sigltil2 = 0,
            sigrtil1 = 0, sigrtil2 = 0, ctil = 0, pd = 0,
            o_par_wrong = 0.5)
  expect_error(math_to_bio(bad2), "names do not match the canonical order")
  
  # Correct names but wrong order
  nm <- rep(1, 9)
  names(nm) <- names(make_mask_names(2))
  bad3 <- nm
  bad3 <- bad3[c(2,1,3:9)]  # swap mu1 and mu2
  expect_error(math_to_bio(bad3), "names do not match the canonical order")
})

test_that("math_to_bio works with create_param_vector_masked output", {
  p <- 2
  full <- create_param_vector_masked(
    param_vector = examples$optim_par_vec,
    mask = NULL,
    p = p
  )
  bio <- math_to_bio(full)
  expect_type(bio, "list")
  expect_true(all(c("mu","sigltil","sigrtil","ctil","pd","o_mat") %in% names(bio)))
})
