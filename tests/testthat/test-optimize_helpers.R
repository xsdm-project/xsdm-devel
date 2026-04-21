# tests/testthat/test-optimize_helpers.R

test_that("optimize_loglik_math_ injects invh_lt and restores names if optimizer drops them", {
  env_dat <- array(0, dim = c(2, 2, 1))
  occ     <- c(1L, 0L)
  par0    <- c(a = 0.1, b = -0.2, c = 0.3)
  invh    <- 1:3
  
  # Mock optimizer: check control has invhessian.lt, return par WITHOUT names
  fake_opt <- function(par, fn, env_dat, mask, occ, negative, num_threads, control, hessian) {
    expect_identical(control$invhessian.lt, invh)      # covers invh_lt injection
    list(par = unname(par), value = 123, convergence = 0L)
  }
  
  res <- xsdm:::optimize_loglik_math_(
    param_vector = par0,
    env_dat      = env_dat,
    occ          = occ,
    mask         = NULL,
    num_threads  = 1L,
    base_control = list(maxeval = 5),
    invh_lt      = invh,
    optimizer_fun = fake_opt
  )
  
  # Names must be restored
  expect_identical(names(res$par), names(par0))
  expect_equal(res$convergence, 0L)
})

test_that("optimize_loglik_math_ sanitizes non-finite par and sets convergence = -99", {
  env_dat <- array(0, dim = c(2, 2, 1))
  occ     <- c(1L, 0L)
  par0    <- c(a = 0.1, b = -0.2, c = 0.3)
  
  # Mock optimizer returning NaN/Inf; triggers sanitization branch
  fake_opt <- function(par, ...) {
    list(par = c(a = NaN, b = Inf, c = 1), value = 5, convergence = 1L)
  }
  
  res <- xsdm:::optimize_loglik_math_(
    param_vector = par0,
    env_dat      = env_dat,
    occ          = occ,
    mask         = NULL,
    num_threads  = 1L,
    base_control = list(maxeval = 5),
    optimizer_fun = fake_opt
  )
  
  # Non-finite entries replaced with starting values
  expect_equal(res$par["a"], par0["a"])
  expect_equal(res$par["b"], par0["b"])
  # Finite entry preserved
  expect_equal(unname(res$par["c"]), 1)
  # Convergence forced to -99L
  expect_identical(res$convergence, -99L)
})

test_that("optimize_loglik_math_ returns structured error object on optimizer failure", {
  env_dat <- array(0, dim = c(2, 2, 1))
  occ     <- c(1L, 0L)
  par0    <- c(a = 0.1, b = -0.2)
  
  # Mock optimizer that throws
  fake_opt <- function(...) stop("boom")
  
  out <- xsdm:::optimize_loglik_math_(
    param_vector = par0,
    env_dat      = env_dat,
    occ          = occ,
    mask         = NULL,
    num_threads  = 1L,
    base_control = list(maxeval = 5),
    optimizer_fun = fake_opt
  )
  
  expect_s3_class(out, "ucminf_error")
  expect_equal(out$par, par0)
  expect_equal(out$value, Inf)
  expect_identical(out$convergence, NA_integer_)
  expect_match(out$error, "boom")
})

