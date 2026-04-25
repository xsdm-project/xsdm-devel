# tests/testthat/test-canonicalize_solutions.R

test_that("canonicalize_param_vector: identity for p = 1", {
  par <- c(mu1 = 11, sigltil1 = log(1.2), sigrtil1 = log(0.8),
           ctil = -6.7, pd = 0.85)
  expect_identical(canonicalize_param_vector(par, par), par)
})

test_that("canonicalize_param_vector: equivalence-class images map to same loglik", {
  data("examples", package = "xsdm")
  ref_math <- examples$par_vec  # canonical truth, p = 2

  ll_at <- function(v) {
    -xsdm:::loglik_math(
      param_vector = v,
      env_dat      = examples$env_array,
      occ          = examples$occ_vec,
      mask         = NULL,
      negative     = TRUE
    )
  }

  ll_ref <- ll_at(ref_math)

  # Build EC representatives on the bio scale and round-trip back to math.
  # Only representatives with det(o_mat) = +1 are math-scale-representable
  # (the math-scale `o_par` parameterises SO(p), not O(p)). For p = 2 those
  # are: (flip = (0,0), perm = (1,2)) — identity — and
  # (flip = (1,1), perm = (2,1)) — full sign flip + swap.
  bio <- math_to_bio(ref_math)
  # Math-scale-representable EC ops for p = 2 (det(o_mat) = +1):
  #   (flip = (1, 1), perm = (1, 2)) -- flip both columns
  #   (flip = (1, 0), perm = (2, 1)) -- flip col 1 + swap
  #   (flip = (0, 1), perm = (2, 1)) -- flip col 2 + swap
  bio_v1 <- convert_equivalence_class(bio, flip = c(1, 1), perm = c(1, 2))
  bio_v2 <- convert_equivalence_class(bio, flip = c(1, 0), perm = c(2, 1))
  bio_v3 <- convert_equivalence_class(bio, flip = c(0, 1), perm = c(2, 1))

  fold_bio <- function(b_partial) list(
    mu = bio$mu, ctil = bio$ctil, pd = bio$pd,
    sigltil = b_partial$sigltil, sigrtil = b_partial$sigrtil,
    o_mat = b_partial$o_mat
  )
  m_v1 <- bio_to_math(fold_bio(bio_v1))
  m_v2 <- bio_to_math(fold_bio(bio_v2))
  m_v3 <- bio_to_math(fold_bio(bio_v3))

  for (cand in list(m_v1, m_v2, m_v3)) {
    expect_equal(ll_at(cand), ll_ref, tolerance = 1e-8)
    canon <- canonicalize_param_vector(cand, reference = ref_math)
    expect_equal(ll_at(canon), ll_ref, tolerance = 1e-8)
    d_in    <- xsdm:::distance_between_params(math_to_bio(cand),  bio)
    d_canon <- xsdm:::distance_between_params(math_to_bio(canon), bio)
    expect_lte(d_canon, d_in + 1e-10)
  }
})

test_that("canonicalize_solutions: aligns every solution to the best fit", {
  data("examples", package = "xsdm")

  # Stub a fit-shaped object with the truth and its non-trivial
  # math-scale-representable EC image.
  ref_math <- examples$par_vec
  bio <- math_to_bio(ref_math)
  bio_v1 <- convert_equivalence_class(bio, flip = c(1, 1), perm = c(1, 2))
  bio_v2 <- convert_equivalence_class(bio, flip = c(1, 0), perm = c(2, 1))
  fold_bio <- function(b_partial) list(
    mu = bio$mu, ctil = bio$ctil, pd = bio$pd,
    sigltil = b_partial$sigltil, sigrtil = b_partial$sigrtil,
    o_mat   = b_partial$o_mat
  )
  variants <- list(ref_math, bio_to_math(fold_bio(bio_v1)),
                   bio_to_math(fold_bio(bio_v2)))

  fit <- list(
    solutions = data.frame(
      start_id    = seq_along(variants),
      loglik      = c(-100, -100.0001, -100.0002),
      convergence = c(1L, 1L, 1L)
    ),
    best = list(par = ref_math, loglik = -100, convergence = 1L)
  )
  fit$solutions$full_par <- variants

  canon <- canonicalize_solutions(fit)

  for (i in seq_along(canon$solutions$full_par)) {
    expect_equal(
      unname(canon$solutions$full_par[[i]]),
      unname(canon$best$par),
      tolerance = 1e-6
    )
  }
})
