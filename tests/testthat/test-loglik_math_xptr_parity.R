# tests/testthat/test-loglik_math_xptr_parity.R
#
# End-to-end parity: the pure-C++ XPtr, when driven through
# ucminfcpp::ucminf_xptr, produces an objective value equal to
# R::loglik_math(par, env_dat, occ, negative = TRUE) evaluated at the
# same par. This is the defining guarantee of the new path:
# whatever the optimizer reports as `value`, it is exactly
# -loglik(par).
#
# Note: the legacy R-callback closure (make_loglik_math_xptr) is kept
# in src/loglik_math_xptr.cpp for backwards compatibility, but cannot
# be exercised directly through ucminf_xptr (ucminfcpp strips the
# names of `par` before invoking the ObjFun, and the R-callback
# closure forwards an unnamed vector to create_param_vector_masked,
# which requires names). The new C++ closure avoids this by taking
# free_names at construction and splicing positional values into
# canonical slots internally.

test_that("pure-C++ XPtr: ucminf_xptr objective equals R::loglik_math at par", {
  set.seed(20260201)
  n <- 10
  Tt <- 5
  p <- 2
  env_dat <- array(runif(n * Tt * p, -1, 1), dim = c(n, Tt, p))
  occ <- rep(c(1L, 0L), length.out = n)

  par0 <- c(mu1 = 0.0, mu2 = 0.0,
            sigltil1 = 0.0, sigltil2 = 0.0,
            sigrtil1 = 0.0, sigrtil2 = 0.0,
            ctil = 0.0, pd = 0.0,
            o_par1 = 0.0)

  xptr_cpp <- xsdm:::make_loglik_math_xptr_cpp(
    env_dat     = env_dat,
    occ         = as.integer(occ),
    mask        = NULL,
    free_names  = names(par0),
    num_threads = 1L,
    grad        = "central",
    gradstep    = c(1e-6, 1e-8)
  )

  res <- ucminfcpp::ucminf_xptr(
    par     = par0,
    xptr    = xptr_cpp,
    control = list(maxeval = 50, xtol = 1e-8, stepmax = 1),
    hessian = 0
  )

  # ucminf_xptr drops names; restore them for the reference call.
  par_final <- res$par
  names(par_final) <- names(par0)

  f_r <- loglik_math(par_final, env_dat, as.integer(occ), negative = TRUE)
  expect_equal(res$value, f_r, tolerance = 1e-8)
})

test_that("pure-C++ XPtr: initial objective agrees with loglik_math_cpp", {
  set.seed(20260203)
  n <- 6
  Tt <- 3
  p <- 2
  env_dat <- array(runif(n * Tt * p, -1, 1), dim = c(n, Tt, p))
  occ <- rep(c(1L, 0L), length.out = n)

  par0 <- c(mu1 = 0.1, mu2 = -0.1,
            sigltil1 = 0.2, sigltil2 = -0.2,
            sigrtil1 = 0.3, sigrtil2 = -0.3,
            ctil = -0.4, pd = 0.5,
            o_par1 = 0.1)

  # Single-iteration run gives us `res$value = nll` at near-par0.
  xptr_cpp <- xsdm:::make_loglik_math_xptr_cpp(
    env_dat     = env_dat,
    occ         = as.integer(occ),
    mask        = NULL,
    free_names  = names(par0),
    num_threads = 1L,
    grad        = "central",
    gradstep    = c(1e-6, 1e-8)
  )
  res <- ucminfcpp::ucminf_xptr(
    par     = par0,
    xptr    = xptr_cpp,
    control = list(maxeval = 1L, xtol = 1e-12, stepmax = 1),
    hessian = 0
  )
  par_final <- res$par
  names(par_final) <- names(par0)

  f_ref <- loglik_math_cpp(par_final, env_dat, as.integer(occ),
                           negative = TRUE)
  expect_equal(res$value, f_ref, tolerance = 1e-8)
})

test_that("pure-C++ XPtr: unnamed mask is rejected with a clear error", {
  set.seed(20260205)
  n <- 6; Tt <- 3; p <- 2
  env_dat <- array(runif(n * Tt * p, -1, 1), dim = c(n, Tt, p))
  occ <- rep(c(1L, 0L), length.out = n)
  unnamed_mask <- c(0.0, 0.0)  # no names

  expect_error(
    xsdm:::make_loglik_math_xptr_cpp(
      env_dat     = env_dat,
      occ         = as.integer(occ),
      mask        = unnamed_mask,
      free_names  = c("mu1", "mu2", "sigltil1", "sigltil2",
                      "sigrtil1", "sigrtil2", "ctil", "o_par1"),
      num_threads = 1L,
      grad        = "central",
      gradstep    = c(1e-6, 1e-8)
    ),
    "must be a named numeric vector"
  )
})

test_that("pure-C++ XPtr: invalid gradstep is rejected", {
  n <- 4; Tt <- 2; p <- 2
  env_dat <- array(0.0, dim = c(n, Tt, p))
  occ <- rep(c(1L, 0L), length.out = n)
  expect_error(
    xsdm:::make_loglik_math_xptr_cpp(
      env_dat     = env_dat,
      occ         = as.integer(occ),
      mask        = NULL,
      free_names  = c("mu1"),
      num_threads = 1L,
      grad        = "central",
      gradstep    = c(0, 1e-8)
    ),
    "gradstep.*strictly positive"
  )
})

test_that("pure-C++ XPtr: works with a non-empty mask", {
  set.seed(20260204)
  n <- 6
  Tt <- 3
  p <- 2
  env_dat <- array(runif(n * Tt * p, -1, 1), dim = c(n, Tt, p))
  occ <- rep(c(1L, 0L), length.out = n)

  # free = all except pd; mask sets pd = 0.
  free0 <- c(mu1 = 0.1, mu2 = -0.1,
             sigltil1 = 0.0, sigltil2 = 0.0,
             sigrtil1 = 0.0, sigrtil2 = 0.0,
             ctil = -0.5,
             o_par1 = 0.0)
  mask <- c(pd = 0.0)

  xptr_cpp <- xsdm:::make_loglik_math_xptr_cpp(
    env_dat     = env_dat,
    occ         = as.integer(occ),
    mask        = mask,
    free_names  = names(free0),
    num_threads = 1L,
    grad        = "central",
    gradstep    = c(1e-6, 1e-8)
  )
  res <- ucminfcpp::ucminf_xptr(
    par     = free0,
    xptr    = xptr_cpp,
    control = list(maxeval = 1L, xtol = 1e-12, stepmax = 1),
    hessian = 0
  )
  par_final <- res$par
  names(par_final) <- names(free0)

  f_ref <- loglik_math_cpp(par_final, env_dat, as.integer(occ),
                           mask = mask, negative = TRUE)
  expect_equal(res$value, f_ref, tolerance = 1e-8)
})

test_that("pure-C++ XPtr: closure-side errors surface as sentinel +Inf, not a crash", {
  # Simulate a failure inside the ObjFun by constructing an XPtr whose
  # canonical-slot coverage is incomplete. The lazy-init throw on the
  # first invocation must be caught at the closure boundary and reported
  # as +Inf rather than propagating through ucminfcpp's template stack.
  set.seed(20260207)
  n <- 6; Tt <- 3; p <- 2
  env_dat <- array(runif(n * Tt * p, -1, 1), dim = c(n, Tt, p))
  occ <- rep(c(1L, 0L), length.out = n)

  # Deliberately missing canonical slots — not all of 3p+2+q are covered.
  bad_free <- c("mu1", "mu2", "ctil")

  xptr <- xsdm:::make_loglik_math_xptr_cpp(
    env_dat     = env_dat,
    occ         = as.integer(occ),
    mask        = NULL,
    free_names  = bad_free,
    num_threads = 1L,
    grad        = "central",
    gradstep    = c(1e-6, 1e-8)
  )

  # ucminfcpp invokes the closure: lazy-init throws, closure catches
  # and sets f = +Inf. ucminfcpp returns gracefully.
  par0 <- setNames(rep(0.0, length(bad_free)), bad_free)
  res <- ucminfcpp::ucminf_xptr(
    par     = par0,
    xptr    = xptr,
    control = list(maxeval = 2, stepmax = 1, grtol = 1e-6, xtol = 1e-12),
    hessian = FALSE
  )
  # No crash. The objective at the stopping point should be +Inf and
  # convergence should be a non-success code.
  expect_true(is.infinite(res$value) || !is.finite(res$value))
})
