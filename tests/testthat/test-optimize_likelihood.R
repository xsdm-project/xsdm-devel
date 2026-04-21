# tests/testthat/test-optimize_likelihood.R
library(testthat)
library(xsdm)
test_that("optimize_likelihood() returns expected structure and sorting", {
  # Small synthetic env array: n locations, T time steps, p variables
  set.seed(1)
  n <- 12
  Tt <- 6
  p <- 2
  env_dat <- array(runif(n * Tt * p, min = -2, max = 2), dim = c(n, Tt, p))

  # Occ vector with at least some presences (1) and absences (0)
  occ <- rep(c(1L, 0L, 1L, 0L), length.out = n)

  # Keep it fast: few starts and low maxeval
  set.seed(123)
  res <- optimize_likelihood(
    env_dat = env_dat,
    occ = occ,
    mask = NULL,
    num_starts = 4L,
    breadth = 1,
    parallel = FALSE,
    num_threads = 1L,
    control = list(maxeval = 20, stepmax = 1, xtol = 1e-8),
    verbose = FALSE
  )

  # Structure: list with named elements
  expect_type(res, "list")
  expect_named(res, c("solutions", "best"))

  # 'solutions' is a data.frame with required columns
  expect_s3_class(res$solutions, "data.frame")
  expect_true(all(c("start_id", "loglik", "convergence", "full_par") %in% names(res$solutions)))

  # rows <= num_starts (all succeed on valid data); 'full_par' is a list-column
  expect_lte(nrow(res$solutions), 4L)
  expect_gte(nrow(res$solutions), 1L)
  expect_true(is.list(res$solutions$full_par))
  expect_equal(length(res$solutions$full_par), nrow(res$solutions))

  # sorted by decreasing loglik
  expect_true(all(diff(res$solutions$loglik) <= 0))

  # 'best' equals the first row of 'solutions'
  expect_equal(res$best$loglik, res$solutions$loglik[1])
  expect_equal(res$best$convergence, res$solutions$convergence[1])
  expect_identical(res$best$par, res$solutions$full_par[[1]])

  # Each full_par is a named numeric vector of canonical length with no NA
  canon_names <- names(make_mask_names(p))
  canon_len <- length(canon_names)
  for (i in seq_len(nrow(res$solutions))) {
    fp <- res$solutions$full_par[[i]]
    expect_true(is.numeric(fp))
    expect_setequal(names(fp), canon_names)
    expect_equal(length(fp), canon_len)
    expect_false(anyNA(fp))
  }
})


test_that("optimize_likelihood() is reproducible under fixed seed", {
  set.seed(3)
  n <- 8
  Tt <- 5
  p <- 2
  env_dat <- array(runif(n * Tt * p, min = -0.5, max = 1.5), dim = c(n, Tt, p))
  occ <- rep(c(1L, 0L), length.out = n)

  set.seed(2026)
  a <- optimize_likelihood(env_dat, occ,
    num_starts = 3L,
    num_threads = 1L, parallel = FALSE,
    control = list(maxeval = 30), verbose = FALSE
  )

  set.seed(2026)
  b <- optimize_likelihood(env_dat, occ,
    num_starts = 3L,
    num_threads = 1L, parallel = FALSE,
    control = list(maxeval = 30), verbose = FALSE
  )

  # Determinism: starting points (Sobol') + ucminf path should yield identical results
  # for identical seeds and inputs
  expect_identical(a$solutions$loglik, b$solutions$loglik)
  expect_identical(a$solutions$convergence, b$solutions$convergence)
  # 'full_par' entries should match exactly as well
  expect_identical(a$solutions$full_par, b$solutions$full_par)
  expect_identical(a$best, b$best)
})

test_that("optimize_likelihood() errors when there are no presences", {
  set.seed(4)
  n <- 10
  Tt <- 6
  p <- 2
  env_dat <- array(runif(n * Tt * p), dim = c(n, Tt, p))
  occ <- rep(0L, n) # no presences

  expect_error(
    optimize_likelihood(env_dat, occ, num_starts = 2L, num_threads = 1L),
    regexp = "No presences \\(occ==1\\) available",
    fixed = FALSE
  )
})

test_that("optimize_likelihood() prints the single-thread notice in parallel mode", {
  skip_on_cran() # future/callr parallelism can be flaky on CRAN
  set.seed(5)
  n <- 6
  Tt <- 4
  p <- 2
  env_dat <- array(runif(n * Tt * p), dim = c(n, Tt, p))
  occ <- rep(c(1L, 0L), length.out = n)

  expect_message(
    optimize_likelihood(env_dat,
                        occ,
                        num_starts = 4L,
                        parallel = TRUE,
                        num_threads = 2L,
                        verbose = TRUE,
      control = list(maxeval = 5)
    ),
    regexp = "parallel=TRUE: forcing num_threads=1",
    fixed = FALSE
  )
})

test_that("optimize_likelihood() respects mask and reconstructs full parameter vectors", {
  set.seed(2)
  n <- 10
  Tt <- 6
  p <- 2
  env_dat <- array(runif(n * Tt * p, min = -1, max = 2), dim = c(n, Tt, p))
  occ <- c(1L, 0L, 1L, rep(0L, n - 3))

  # Mask on the math scale: fix pd (math=0 => bio=0.5) and mu2 = 0.25
  mask <- c(pd = 0, mu2 = 0.25)

  set.seed(999)
  res <- optimize_likelihood(
    env_dat     = env_dat,
    occ         = occ,
    mask        = mask,
    num_starts  = 4L,
    breadth     = 1,
    parallel    = FALSE,
    num_threads = 1L,
    control     = list(maxeval = 5),
    verbose     = FALSE
  )

  # Full vectors must include all canonical names
  canon_names <- names(make_mask_names(p))
  expect_true(all(vapply(res$solutions$full_par, function(x) setequal(names(x), canon_names), logical(1))))

  # Masked values must be present and equal in the reconstructed full vectors
  for (fp in res$solutions$full_par) {
    expect_equal(unname(fp["pd"]), mask[["pd"]])
    expect_equal(unname(fp["mu2"]), mask[["mu2"]])
  }
})

test_that("logical occ is normalized to integer (coverage of is.logical(occ))", {
  set.seed(2)
  n <- 10
  Tt <- 6
  p <- 2
  env_dat <- array(runif(n * Tt * p, min = -1, max = 2), dim = c(n, Tt, p))
  occ <- c(TRUE, FALSE, TRUE, rep(FALSE, n - 3))
  # Mask on the math scale: fix pd (math=0 => bio=0.5) and mu2 = 0.25
  mask <- c(pd = 0, mu2 = 0.25)

  res <- optimize_likelihood(
    env_dat     = env_dat,
    occ         = occ,
    mask        = mask,
    num_starts  = 4L,
    breadth     = 1,
    parallel    = FALSE,
    num_threads = 1L,
    control     = list(maxeval = 5),
    verbose     = FALSE
  )

  expect_type(res, "list")
  expect_named(res, c("solutions", "best"))
})

test_that("start_parms() returning 0 rows triggers informative error", {
  set.seed(12)
  n <- 8
  Tt <- 4
  p <- 2
  env_dat <- array(runif(n * Tt * p), dim = c(n, Tt, p))
  occ <- rep(c(1L, 0L), length.out = n)

  # Mock start_parms() inside the package to return an empty data.frame
  testthat::with_mocked_bindings(
    start_parms = function(env_dat,
                           mask = NULL,
                           breadth = 1,
                           num_starts = 1L) {
      # 0 rows, 0 cols -> nrow(...) == 0 -> error path in optimize_likelihood()
      data.frame()[FALSE, ]
    },
    {
      expect_error(
        optimize_likelihood(env_dat, occ, num_starts = 1L, num_threads = 1L, control = list()),
        regexp = "start_parms\\(\\) returned no starting points\\.",
        fixed = FALSE
      )
    },
    .package = "xsdm" # ensure binding is replaced in the package namespace
  )
})

test_that("verbose messages are printed in non-parallel mode", {
  set.seed(13)
  n <- 9
  Tt <- 4
  p <- 2
  env_dat <- array(runif(n * Tt * p, -0.5, 0.5), dim = c(n, Tt, p))
  occ <- rep(c(1L, 0L, 1L), length.out = n)

  # Two separate runs to match both messages cleanly with small workloads
  expect_message(
    optimize_likelihood(
      env_dat, occ,
      num_starts = 4L,
      num_threads = 1L,
      control = list(maxeval = 20),
      verbose = TRUE, # covers: "Optimizing from %d starting points%s."
      parallel = FALSE
    ),
    regexp = "Optimizing from 4 starting points\\.",
    fixed = FALSE
  )

  expect_message(
    optimize_likelihood(
      env_dat, occ,
      num_starts = 4L,
      num_threads = 1L,
      control     = list(maxeval = 5),
      verbose = TRUE, # covers: "Best log-likelihood: ..."
      parallel = FALSE
    ),
    regexp = "Best log-likelihood:",
    fixed = FALSE
  )
})

test_that("parallel branch executes, forces single-thread, and prints messages", {
  # These packages are Suggested; skip if missing to keep CI robust.
  skip_if_not_installed("furrr")
  skip_if_not_installed("future.callr")

  set.seed(14)
  n <- 6
  Tt <- 4
  p <- 2
  env_dat <- array(runif(n * Tt * p, -1, 1), dim = c(n, Tt, p))
  occ <- rep(c(1L, 0L), length.out = n)

  # 1) Covers the 'num_threads <- 1L' line and its message
  expect_message(
    optimize_likelihood(
      env_dat, occ,
      num_starts = 4L,
      parallel = TRUE, # enter parallel branch
      num_threads = 4L, # triggers "forcing num_threads=1 ..." + assignment
      control = NULL, # covers: if (is.null(control)) control <- list()
      verbose = TRUE
    ),
    regexp =
      "parallel=TRUE: forcing num_threads=1 to avoid nested RcppParallel\\.",
    fixed = FALSE
  )

  # 2) Also ensure we print the "(parallel)" message in the
  # 'Optimizing from ...' line
  expect_message(
    optimize_likelihood(
      env_dat, occ,
      num_starts = 4L,
      parallel = TRUE,
      num_threads = 2L,
      control = list(maxeval = 5),
      verbose = TRUE
    ),
    regexp = "Optimizing from 4 starting points \\(parallel\\)\\.",
    fixed = FALSE
  )
})

# ── Fix 1 & 2: runner tryCatch + NULL filtering ───────────────────────────────

test_that("optimize_likelihood() completes and excludes failed starts when runner errors", {
  # Inject a runner that raises an error for every other start by mocking
  # optimize_loglik_math_ so that odd-indexed calls throw.
  set.seed(42)
  n <- 10; Tt <- 5; p <- 2
  env_dat <- array(runif(n * Tt * p, -1, 1), dim = c(n, Tt, p))
  occ <- rep(c(1L, 0L), length.out = n)

  call_count <- 0L
  fake_opt <- function(par, fn, env_dat, mask, occ, negative, num_threads, control, hessian) {
    call_count <<- call_count + 1L
    if (call_count %% 2L == 1L) stop("injected runner failure")
    # successful path: return a minimal valid ucminf result
    list(par = par, value = 10, convergence = 0L)
  }

  # Temporarily replace the internal optimizer so odd starts always fail
  # without the need to corrupt env_dat values.
  with_mocked_bindings <- testthat::with_mocked_bindings
  with_mocked_bindings(
    `optimize_loglik_math_` = function(param_vector, env_dat, occ, mask,
                                       num_threads, base_control, invh_lt) {
      call_count <<- call_count + 1L
      if (call_count %% 2L == 1L) stop("injected runner failure")
      list(par = param_vector, value = 10, convergence = 0L)
    },
    {
      res <- optimize_likelihood(
        env_dat     = env_dat,
        occ         = occ,
        num_starts  = 4L,
        num_threads = 1L,
        parallel    = FALSE,
        control     = list(maxeval = 5),
        verbose     = FALSE
      )

      # Must return a valid result with fewer rows than num_starts
      expect_type(res, "list")
      expect_named(res, c("solutions", "best"))
      expect_s3_class(res$solutions, "data.frame")
      expect_lt(nrow(res$solutions), 4L)
      expect_gte(nrow(res$solutions), 1L)
    },
    .package = "xsdm"
  )
})

test_that("optimize_likelihood() reports failed-start count with verbose=TRUE", {
  set.seed(77)
  n <- 8; Tt <- 4; p <- 2
  env_dat <- array(runif(n * Tt * p, -1, 1), dim = c(n, Tt, p))
  occ <- rep(c(1L, 0L), length.out = n)

  call_count <- 0L
  testthat::with_mocked_bindings(
    `optimize_loglik_math_` = function(param_vector, env_dat, occ, mask,
                                       num_threads, base_control, invh_lt) {
      call_count <<- call_count + 1L
      if (call_count == 1L) stop("first start fails")
      list(par = param_vector, value = 10, convergence = 0L)
    },
    {
      expect_message(
        optimize_likelihood(
          env_dat     = env_dat,
          occ         = occ,
          num_starts  = 3L,
          num_threads = 1L,
          parallel    = FALSE,
          control     = list(maxeval = 5),
          verbose     = TRUE
        ),
        regexp = "1 of 3 starting point failed and were excluded\\.",
        fixed = FALSE
      )
    },
    .package = "xsdm"
  )
})

