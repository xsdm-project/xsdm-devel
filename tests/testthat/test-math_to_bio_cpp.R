# tests/testthat/test-math_to_bio_cpp.R
#
# Parity: xsdm:::.math_to_bio_cpp must agree element-wise with math_to_bio().

test_that(".math_to_bio_cpp: matches R on p = 1 (no o_par)", {
  pv <- make_mask_names(1)
  pv[] <- c(mu1 = 0.7, sigltil1 = log(1.5), sigrtil1 = log(0.8),
            ctil = -1.3, pd = 0.2)
  ref <- math_to_bio(pv)
  cpp <- xsdm:::.math_to_bio_cpp(pv)

  expect_equal(cpp$mu,      ref$mu,      tolerance = 1e-14)
  expect_equal(cpp$sigltil, ref$sigltil, tolerance = 1e-14)
  expect_equal(cpp$sigrtil, ref$sigrtil, tolerance = 1e-14)
  expect_equal(cpp$ctil,    ref$ctil,    tolerance = 1e-14)
  expect_equal(cpp$pd,      ref$pd,      tolerance = 1e-14)
  expect_equal(cpp$o_mat,   ref$o_mat,   tolerance = 1e-14)
})

test_that(".math_to_bio_cpp: matches R on examples$par_vec (p = 2)", {
  pv  <- examples$par_vec
  ref <- math_to_bio(pv)
  cpp <- xsdm:::.math_to_bio_cpp(pv)

  expect_equal(cpp$mu,      ref$mu,      tolerance = 1e-12)
  expect_equal(cpp$sigltil, ref$sigltil, tolerance = 1e-12)
  expect_equal(cpp$sigrtil, ref$sigrtil, tolerance = 1e-12)
  expect_equal(cpp$ctil,    ref$ctil,    tolerance = 1e-12)
  expect_equal(cpp$pd,      ref$pd,      tolerance = 1e-12)
  expect_equal(cpp$o_mat,   ref$o_mat,   tolerance = 1e-12)
})

test_that(".math_to_bio_cpp: matches R on random p = 3 draws", {
  set.seed(42)
  for (trial in 1:5) {
    nm <- make_mask_names(3)
    nm[] <- rnorm(length(nm))
    ref <- math_to_bio(nm)
    cpp <- xsdm:::.math_to_bio_cpp(nm)
    expect_equal(cpp$mu,      ref$mu,      tolerance = 1e-10)
    expect_equal(cpp$sigltil, ref$sigltil, tolerance = 1e-10)
    expect_equal(cpp$sigrtil, ref$sigrtil, tolerance = 1e-10)
    expect_equal(cpp$ctil,    ref$ctil,    tolerance = 1e-10)
    expect_equal(cpp$pd,      ref$pd,      tolerance = 1e-10)
    expect_equal(cpp$o_mat,   ref$o_mat,   tolerance = 1e-10)
  }
})

test_that(".math_to_bio_cpp: errors on out-of-order names", {
  pv <- examples$par_vec
  perm <- pv[c(2, 1, 3:length(pv))]
  expect_error(xsdm:::.math_to_bio_cpp(perm), "canonical order")
})

test_that(".build_canonical_param_vector_cpp: matches create_param_vector_masked", {
  p2 <- 2
  pv2 <- c(sigltil1 = 1.0, sigltil2 = 1.1, sigrtil1 = 2.0, sigrtil2 = 2.2,
           ctil = 0.3, o_par1 = 0.0)
  mask2 <- c(mu1 = 0.1, mu2 = 0.2, pd = 0.05)

  ref <- create_param_vector_masked(pv2, mask2, p2)
  cpp <- xsdm:::.build_canonical_param_vector_cpp(pv2, mask2, p2)

  # same values, same canonical order
  expect_equal(unname(cpp), unname(ref), tolerance = 1e-14)
  expect_identical(names(cpp), names(ref))
})

test_that(".build_canonical_param_vector_cpp: rejects overlap and missing", {
  p2 <- 2
  pv2 <- c(mu1 = 0.1, sigltil1 = 1.0, sigltil2 = 1.1,
           sigrtil1 = 2.0, sigrtil2 = 2.2, ctil = 0.3,
           pd = 0.0, o_par1 = 0.0)
  mask2 <- c(mu1 = 0.5)  # overlaps with param_vector
  expect_error(
    xsdm:::.build_canonical_param_vector_cpp(pv2, mask2, p2),
    "overlap"
  )

  # missing mu2
  bad <- pv2   # no mu2
  expect_error(
    xsdm:::.build_canonical_param_vector_cpp(bad, NULL, p2),
    "Missing canonical"
  )
})
