# Parity test: pure-R math_to_bio_r vs canonical math_to_bio (C++ wrapper).
# Tolerance is 1e-6 (ecological precision); typical observed diff is at
# the level of floating-point round-off (~1e-15).

tol <- 1e-6

draw_canonical_pv <- function(p, seed = 1L) {
  set.seed(seed)
  q <- p * (p - 1L) / 2L
  total <- 3L * p + 2L + q
  vals <- stats::runif(total, min = -2, max = 2)
  vals[(p + 1L):(2L * p)]       <- abs(vals[(p + 1L):(2L * p)]) + 0.1   # log(sigltil)
  vals[(2L * p + 1L):(3L * p)]  <- abs(vals[(2L * p + 1L):(3L * p)]) + 0.1
  names(vals) <- names(make_mask_names(p))
  vals
}

expect_bio_equal <- function(a, b, tol) {
  expect_equal(a$mu,      b$mu,      tolerance = tol)
  expect_equal(a$sigltil, b$sigltil, tolerance = tol)
  expect_equal(a$sigrtil, b$sigrtil, tolerance = tol)
  expect_equal(a$ctil,    b$ctil,    tolerance = tol)
  expect_equal(a$pd,      b$pd,      tolerance = tol)
  expect_equal(unname(as.numeric(a$o_mat)),
               unname(as.numeric(b$o_mat)), tolerance = tol)
}

test_that("math_to_bio: examples fixture parity (p = 2)", {
  pv <- examples$par_vec
  cpp <- math_to_bio(pv)
  r   <- xsdm:::math_to_bio_r(pv)
  expect_bio_equal(cpp, r, tol)
})

test_that("math_to_bio: random sweep p in {1, 2, 3}", {
  for (p in 1:3) {
    for (seed in 1:5) {
      pv  <- draw_canonical_pv(p, seed)
      cpp <- math_to_bio(pv)
      r   <- xsdm:::math_to_bio_r(pv)
      expect_bio_equal(cpp, r, tol)
    }
  }
})

test_that("math_to_bio: identity mu and ctil are preserved", {
  pv <- draw_canonical_pv(2L, seed = 11L)
  cpp <- math_to_bio(pv)
  expect_equal(cpp$mu,   unname(pv[c("mu1", "mu2")]), tolerance = tol)
  expect_equal(cpp$ctil, unname(pv["ctil"]),          tolerance = tol)
})

test_that("math_to_bio: p = 1 returns 1x1 identity o_mat", {
  pv <- draw_canonical_pv(1L, seed = 42L)
  cpp <- math_to_bio(pv)
  expect_equal(dim(cpp$o_mat), c(1L, 1L))
  expect_equal(as.numeric(cpp$o_mat), 1.0, tolerance = tol)
})

test_that("math_to_bio: rejects unnamed vectors (R wrapper -> C++)", {
  pv <- draw_canonical_pv(2L, seed = 7L)
  expect_error(math_to_bio(unname(pv)))
})

test_that("math_to_bio: rejects mis-ordered names", {
  pv <- draw_canonical_pv(2L, seed = 9L)
  swapped <- pv[c(2, 1, 3:length(pv))]
  expect_error(math_to_bio(swapped))
})
