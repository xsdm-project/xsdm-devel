library(testthat)
library(xsdmMle)
test_that("bio_to_math is inverse of math_to_bio for small generators (principal branch)", {
  set.seed(1)
  p <- 3
  q <- p * (p - 1) / 2
  # Param vector en escala math dentro de la rama principal
  param_vector <- c(
    setNames(runif(p, -1, 1), paste0("mu", seq_len(p))),
    setNames(runif(p, -0.3, 0.3), paste0("sigltil", seq_len(p))),
    setNames(runif(p, -0.3, 0.3), paste0("sigrtil", seq_len(p))),
    c(ctil = -1.2),
    c(pd = -0.5), # expit(-0.5) ~ 0.377
    setNames(runif(q, -0.5, 0.5), paste0("o_par", seq_len(q)))
  )

  bio <- math_to_bio(param_vector)
  math_back <- bio_to_math(bio)

  # Igualdad exacta en la rama principal
  expect_equal(math_back[names(param_vector)], param_vector, tolerance = 1e-8)
})

test_that("bio_to_math + build_orthogonal_matrix reproduces o_mat for general
          angles", {
  set.seed(2)
  p <- 4
  q <- p * (p - 1) / 2

  # Build a list of parameters in bio scale una lista bio arbitraria
  parms_bio <- list(
    mu = c(10, 12, 5, 7),
    sigltil = c(2, 1.5, 3, 0.8),
    sigrtil = c(1.9, 1.1, 2.8, 0.7),
    ctil = -2,
    pd = 0.85,
    o_mat = {
      o_par <- c(pi/2, 0, 0, 0, 0, pi/2) # bloque pi en (2,1) y (4,3): -I_4
      build_orthogonal_matrix(o_par)
    }
  )

  math <- bio_to_math(parms_bio)
  # reconstruir o_mat desde o_par extraído
  o_par <- math[grep("^o_par\\d+$", names(math))]
  O2 <- build_orthogonal_matrix(o_par)
  expect_equal(O2, parms_bio$o_mat, tolerance = 1e-5)
})

test_that("p=1 case has no o_par and round-trip works", {
  parms_bio <- list(
    mu      = 11,
    sigltil = 2,
    sigrtil = 2.5,
    ctil    = -3,
    pd      = 0.9,
    o_mat   = matrix(1, 1, 1) # 1x1 identity
  )
  math <- bio_to_math(parms_bio)
  expect_true(all(grepl("^(mu1|sigltil1|sigrtil1|ctil|pd)$", names(math))))
  bio <- math_to_bio(math)
  expect_equal(bio$mu, parms_bio$mu)
  expect_equal(bio$sigltil, parms_bio$sigltil)
  expect_equal(bio$sigrtil, parms_bio$sigrtil)
  expect_equal(bio$ctil, parms_bio$ctil)
  expect_equal(bio$pd, parms_bio$pd)
  expect_equal(bio$o_mat, parms_bio$o_mat)
})


test_that("bio_to_math <-> math_to_bio round-trip works for p = 1 (no o_par)", {
  # Fixed values for determinism
  mu <- c(10)
  sigltil <- c(1.2)
  sigrtil <- c(0.8)
  ctil <- 0.3
  pd <- 0.85
  O <- matrix(1, 1, 1) # 1x1 orthogonal

  bio <- list(
    mu = mu, sigltil = sigltil, sigrtil = sigrtil,
    ctil = ctil, pd = pd, o_mat = O
  )

  math <- bio_to_math(bio)

  # Canonical names and length
  expect_identical(names(math), names(make_mask_names(1L)))
  expect_false(any(is.na(math)))

  # Round-trip
  bio_rt <- math_to_bio(math)
  expect_equal(bio_rt$mu, mu, tolerance = 0)
  expect_equal(bio_rt$sigltil, sigltil, tolerance = 0)
  expect_equal(bio_rt$sigrtil, sigrtil, tolerance = 0)
  expect_equal(bio_rt$ctil, ctil, tolerance = 0)
  expect_equal(bio_rt$pd, pd, tolerance = 0)

  # O equality for p=1 is exact
  expect_equal(bio_rt$o_mat, O, tolerance = 0)
})

test_that("bio_to_math <-> math_to_bio round-trip works for p = 3 (with o_par)", {
  # Fixed values for determinism
  p <- 3L
  mu <- c(0.5, -1.2, 2.0)
  sigltil <- c(1.1, 0.7, 2.0)
  sigrtil <- c(1.3, 0.9, 1.5)
  ctil <- -0.3
  pd <- 0.8
  # q = p*(p-1)/2 = 3; choose modest angles to avoid logm instability
  o_par <- c(0.20, -0.15, 0.05)
  O <- build_orthogonal_matrix(o_par)

  # Sanity: O is orthogonal
  expect_true(max(abs(t(O) %*% O - diag(p))) < 1e-12)

  bio <- list(
    mu = mu, sigltil = sigltil, sigrtil = sigrtil,
    ctil = ctil, pd = pd, o_mat = O
  )

  math <- bio_to_math(bio)

  # Canonical names and length
  expect_identical(names(math), names(make_mask_names(p)))
  expect_false(any(is.na(math)))

  # Round-trip
  bio_rt <- math_to_bio(math)
  expect_equal(bio_rt$mu, mu, tolerance = 0)
  expect_equal(bio_rt$sigltil, sigltil, tolerance = 0)
  expect_equal(bio_rt$sigrtil, sigrtil, tolerance = 0)
  expect_equal(bio_rt$ctil, ctil, tolerance = 0)
  expect_equal(bio_rt$pd, pd, tolerance = 0)

  # Orthogonal basis preserved within numerical tolerance
  expect_true(max(abs(bio_rt$o_mat - O)) < 1e-8)
  expect_true(max(abs(t(bio_rt$o_mat) %*% bio_rt$o_mat - diag(p))) < 1e-12)
})

