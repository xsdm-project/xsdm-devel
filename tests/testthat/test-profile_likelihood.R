library(testthat)

test_that("profile_likelihood:
          minimal run returns expected structure (1 left + 1 right)", {
  # Small synthetic data (p = 2) to keep test fast and independent of large
  # datasets
  env_dat <- array(c(
    13, 13, 12, 2, 3, 3,
    15, 16, 15, 3, 5, 3,
    14, 14, 13, 2, 3, 3,
    13, 13, 12, 2, 3, 4
  ), dim = c(4, 3, 2))
  occ <- c(1, 0, 1, 0)

  # Use packaged named vector on math scale (p = 2)
  optim_param_vector <- examples$optim_par_vec

  inc <- 0.2

  expect_silent(
    res <- profile_likelihood(
      profile_parameter = "mu1",
      increment_left = inc,
      increment_right = inc,
      num_steps_left = 1L, # one step left
      num_steps_right = 1L, # one step right
      alpha = 0.95,
      optim_param_vector = optim_param_vector,
      env_dat = env_dat,
      occ = occ,
      mask = NULL,
      num_threads = 1L, # deterministic + fast
      control = list(maxeval = 50),
      verbose = FALSE
    )
  )

  # ---- top-level return structure ----
  expect_true(is.list(res))
  expect_true(all(c("profile", "found_better", "threshold", "parameters") %in% names(res)))

  expect_type(res$found_better, "logical")
  expect_length(res$found_better, 1L)

  expect_type(res$threshold, "double")
  expect_length(res$threshold, 1L)
  expect_true(is.finite(res$threshold))

  # ---- profile data.frame structure ----
  expect_true(is.data.frame(res$profile))
  expect_true(all(c("param", "value_math", "loglik", "convergence") %in%
                    names(res$profile)))

  # Exactly 3 rows: MLE + 1 left step + 1 right step
  expect_equal(nrow(res$profile), 3L)

  # param column should repeat the profiled parameter name
  expect_true(all(res$profile$param == "mu1"))

  # The profiled parameter values should be:
  # mu1, mu1 - inc, mu1 + inc (in that order)
  mu1_0 <- unname(optim_param_vector[["mu1"]])
  expect_equal(
    res$profile$value_math,
    c(
      mu1_0 - inc,
      mu1_0,
      mu1_0 + inc
    )
  )

  # loglik should be finite numeric
  expect_true(is.numeric(res$profile$loglik))
  expect_true(all(is.finite(res$profile$loglik)))

  # convergence: first row is NA (seed MLE row), remaining are integerish or NA
  expect_true(is.na(res$profile$convergence[2]))
  expect_true(all(is.na(res$profile$convergence[-2]) |
                    is.integer(res$profile$convergence[-1]) |
                    is.numeric(res$profile$convergence[-1])))

  # ---- parameters list of full vectors ----
  expect_true(is.list(res$parameters))
  expect_equal(nrow(res$parameters), 3L)

  # Each entry should be a named numeric vector with same names as
  # optim_param_vector
  for (k in seq_along(res$parameters)) {
    expect_true(!is.null(names(res$parameters)[k]))
  }
  # Expect equal names
  expect_true(setequal(
    names(res$parameters),
    names(optim_param_vector)
  ))

  # Check that mu1 is fixed correctly in the stored full vectors:
  expect_equal(unname(res$parameters[1, 1]), mu1_0 - inc)
  expect_equal(unname(res$parameters[2, 1]), mu1_0)
  expect_equal(unname(res$parameters[3, 1]), mu1_0 + inc)
})


test_that("profile_likelihood: MLE row matches direct loglik_math evaluation", {
  env_dat <- array(c(
    13, 13, 12, 2, 3, 3,
    15, 16, 15, 3, 5, 3,
    14, 14, 13, 2, 3, 3,
    13, 13, 12, 2, 3, 4
  ), dim = c(4, 3, 2))
  occ <- c(1, 0, 1, 0)

  opt_vec <- examples$optim_par_vec

  res <- profile_likelihood(
    profile_parameter = "mu1",
    increment_left = 0.2,
    increment_right = 0.2,
    num_steps_left = 1L,
    num_steps_right = 1L,
    alpha = 0.95,
    optim_param_vector = opt_vec,
    env_dat = env_dat,
    occ = occ,
    mask = NULL,
    num_threads = 1L,
    control = list(maxeval = 50),
    verbose = FALSE
  )

  ll0 <- loglik_math(
    param_vector = opt_vec,
    env_dat = env_dat,
    occ = occ,
    mask = NULL,
    negative = FALSE,
    num_threads = 1L
  )

  # check the likelihood of the profiled with the known loglik likelihood
  expect_equal(res$profile$loglik[2], ll0, tolerance = 1e-10)

  # Check that the profiled parameters check with the optim vector
  expect_equal(unlist(res$parameters[2, ]), opt_vec, tolerance = 0)
})


test_that("profile_likelihood: parameters length matches profile rows and
          preserves names", {
  env_dat <- array(0, dim = c(4, 3, 2))
  occ <- c(1, 0, 1, 0)

  res <- profile_likelihood(
    profile_parameter = "mu1",
    increment_left = 0.1,
    increment_right = 0.1,
    num_steps_left = 1L,
    num_steps_right = 1L,
    alpha = 0.95,
    optim_param_vector = examples$optim_par_vec,
    env_dat = env_dat,
    occ = occ,
    mask = NULL,
    num_threads = 1L,
    control = list(maxeval = 20),
    verbose = FALSE
  )

  expect_equal(nrow(res$parameters), nrow(res$profile))

  # Names should be identical to canonical names
  nm <- names(examples$optim_par_vec)
  for (k in 1:nrow(res$parameters)) {
    # thest that each row in data frame is numeris
    expect_true(is.numeric(unlist(res$parameters[k, ])))

    # test the names in the dataframe checks with the optim vector names
    expect_identical(sort(names(res$parameters[k, ])), sort(nm))
  }
})


test_that("profile_likelihood: verbose emits progress messages", {
  env_dat <- array(0, dim = c(4, 3, 2))
  occ <- c(1, 0, 1, 0)

  expect_message(
    profile_likelihood(
      profile_parameter = "mu1",
      increment_left = 0.1,
      increment_right = 0.1,
      num_steps_left = 1L,
      num_steps_right = 1L,
      alpha = 0.95,
      optim_param_vector = examples$optim_par_vec,
      env_dat = env_dat,
      occ = occ,
      mask = NULL,
      num_threads = 1L,
      control = list(maxeval = 10),
      verbose = TRUE
    ),
    regexp = "Start (left|right) side",
    fixed = FALSE
  )
})


test_that("profile_likelihood: errors when no free parameters remain", {
  env_dat <- array(0, dim = c(4, 3, 2))
  occ <- c(1, 0, 1, 0)

  # Mask all parameters except the profiled one -> zero free parameters
  all_names <- names(examples$optim_par_vec)
  mask_names <- setdiff(all_names, "mu1")
  mask <- examples$optim_par_vec[mask_names]

  expect_error(
    profile_likelihood(
      profile_parameter = "mu1",
      increment_left = 0.1,
      num_steps_left = 1L,
      alpha = 0.95,
      optim_param_vector = examples$optim_par_vec,
      env_dat = env_dat,
      occ = occ,
      mask = mask,
      num_threads = 1L,
      control = list(maxeval = 10),
      verbose = FALSE
    ),
    regexp = "No free parameters left to optimize",
    fixed = FALSE
  )
})


test_that("profile_likelihood: input validation errors (bad profile_parameter,
          increments, alpha, env_dat, occ)", {
  env_dat <- array(0, dim = c(4, 3, 2))
  occ <- c(1, 0, 1, 0)

  # profile_parameter not in names
  expect_error(
    profile_likelihood(
      profile_parameter = "does_not_exist",
      increment_left = 0.1,
      num_steps_left = 1L,
      alpha = 0.95,
      optim_param_vector = examples$optim_par_vec,
      env_dat = env_dat,
      occ = occ,
      mask = NULL,
      num_threads = 1L,
      control = list(maxeval = 10)
    ),
    regexp = "profile_parameter|choice",
    fixed = FALSE
  )

  # increment_left must be > 0
  expect_error(
    profile_likelihood(
      profile_parameter = "mu1",
      increment_left = 0,
      num_steps_left = 1L,
      alpha = 0.95,
      optim_param_vector = examples$optim_par_vec,
      env_dat = env_dat,
      occ = occ
    ),
    regexp = "increment_left",
    fixed = FALSE
  )

  # alpha must be in (0,1)
  expect_error(
    profile_likelihood(
      profile_parameter = "mu1",
      increment_left = 0.1,
      num_steps_left = 1L,
      alpha = 1,
      optim_param_vector = examples$optim_par_vec,
      env_dat = env_dat,
      occ = occ
    ),
    regexp = "alpha",
    fixed = FALSE
  )

  # env_dat must be 3D array
  expect_error(
    profile_likelihood(
      profile_parameter = "mu1",
      increment_left = 0.1,
      num_steps_left = 1L,
      alpha = 0.95,
      optim_param_vector = examples$optim_par_vec,
      env_dat = 1:10,
      occ = occ
    ),
    regexp = "env_dat",
    fixed = FALSE
  )

  # occ length mismatch
  expect_error(
    profile_likelihood(
      profile_parameter = "mu1",
      increment_left = 0.1,
      num_steps_left = 1L,
      alpha = 0.95,
      optim_param_vector = examples$optim_par_vec,
      env_dat = env_dat,
      occ = c(1, 0) # wrong length
    ),
    regexp = "occ",
    fixed = FALSE
  )
})


test_that("profile_likelihood: mask != NULL currently errors at baseline if
          overlaps with full optim_param_vector", {
  env_dat <- array(0, dim = c(4, 3, 2))
  occ <- c(1, 0, 1, 0)

  # Any overlap between mask names and full optim_param_vector can break if
  # create_param_vector_masked enforces disjointness.
  # This documents current behavior.
  expect_error(
    profile_likelihood(
      profile_parameter = "mu1",
      increment_left = 0.1,
      increment_right = 0.1,
      num_steps_left = 1L,
      num_steps_right = 1L,
      alpha = 0.95,
      optim_param_vector = examples$optim_par_vec,
      env_dat = env_dat,
      occ = occ,
      mask = c(pd = 0), # overlaps name with full vector
      num_threads = 1L,
      control = list(maxeval = 10),
      verbose = FALSE
    ),
    regexp = "disjoint|complementary|Overlapping names|param_vector.*mask",
    fixed = FALSE
  )
})
