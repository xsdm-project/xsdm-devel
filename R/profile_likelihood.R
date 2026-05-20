#' Helper. Profile one side of a likelihood profile (internal)
#'
#' Fixes one parameter on the math scale and re-optimizes the remaining
#' parameters along a single direction (left/right) until the LR threshold
#' is reached or a step cap is hit.
#'
#' @inheritParams profile_likelihood
#' @param direction Integer. -1 (left) or +1 (right).
#' @param increment Numeric. Step size on math scale for this side.
#' @param max_steps Integer. Maximum iterations for this side.
#' @param base_control Named list. Control passed to \code{ucminfcpp::ucminf_xptr(control = ...)}.
#'   User-specified entries should be merged in the caller (see \code{profile_likelihood}).
#' @param start_full Named numeric. Full warm-start parameter vector.
#'   Defaults to \code{optim_param_vector}.
#' @param invh_lt Optional numeric. Lower triangle of the inverse Hessian for warm-start.
#' @param verbose Logical. If \code{TRUE}, prints compact progress messages; otherwise silent.
#'
#' @return A list with elements \code{ll}, \code{vals}, \code{fulls}, \code{conv},
#'   \code{last_full}, \code{last_invh}, \code{steps}, and \code{crossed}.
#' @keywords internal
#' @seealso \code{\link{profile_likelihood}}
profile_one_side_ <- function(
    direction, # -1 for left, +1 for right
    increment, # step size (math scale)
    max_steps, # iteration cap
    profile_parameter, # name of the profiled parameter
    optim_param_vector, # full parameter template (named)
    env_dat, occ, mask, # data & fixed params
    num_threads, # threads for loglik
    optim_ll, thresh, # MLE loglik & LR threshold
    base_control, # list passed to ucminfcpp
    start_full = optim_param_vector, # warm-start full vector
    invh_lt = NULL, # warm-start inverse Hessian (lower-tri vector)
    verbose = FALSE
) {
  cur_ll <- optim_ll
  nsteps <- 0L
  
  # storage
  ll_vec <- list()
  val_vec <- list()
  full_lst <- list()
  conv_vec <- list()
  
  if (verbose) {
    side_label <- if (direction < 0) "left" else "right"
    message(sprintf("Start %s side", side_label))
  }
  
  # Embedded text update cadence: fine-grained for short runs, coarser for long
  every <- if (max_steps <= 20L) 1L else ceiling(max_steps / 20L)
  
  while (cur_ll > thresh && nsteps < max_steps) {
    nsteps <- nsteps + 1L
    
    # math-scale value for the profiled parameter at this step
    new_val <- optim_param_vector[profile_parameter] + direction * nsteps * increment
    new_mask <- c(stats::setNames(new_val, profile_parameter), mask)
    
    # free parameters = all minus masked
    free_names <- setdiff(names(optim_param_vector), names(new_mask))
    par_init <- start_full[free_names]
    
    # reuse invH if available
    ctrl <- base_control
    if (!is.null(invh_lt)) ctrl$invhessian.lt <- invh_lt

    # Use the pure-C++ XPtr factory (same one optimize_likelihood() uses).
    # The legacy R-callback factory `make_loglik_math_xptr()` does not
    # propagate `free_names` through the ucminfcpp -> loglik_math() boundary,
    # which makes every profile_likelihood() step error out with
    # "param_vector must have names". Switching to make_loglik_math_xptr_cpp
    # both fixes the bug and removes per-evaluation R-callback overhead.
    grad_ctrl <- resolve_xptr_grad_control_(ctrl)
    occ_i <- as.integer(occ)
    loglik_xptr <- make_loglik_math_xptr_cpp(
      env_dat     = env_dat,
      occ         = occ_i,
      mask        = new_mask,
      free_names  = free_names,
      num_threads = num_threads,
      grad        = grad_ctrl$grad,
      gradstep    = grad_ctrl$gradstep
    )

    res <- ucminfcpp::ucminf_xptr(
      par = par_init,
      xptr = loglik_xptr,
      control = ctrl,
      hessian = FALSE
    )
    
    cur_ll <- -res$value
    
    # reconstruct a complete, ordered vector for storage & next warm-start
    full_next <- optim_param_vector
    full_next[names(new_mask)] <- unname(new_mask)
    full_next[free_names] <- res$par
    
    # store
    ll_vec <- c(ll_vec, cur_ll)
    val_vec <- c(val_vec, unname(new_val))
    full_lst <- c(full_lst, list(full_next))
    conv_vec <- c(conv_vec, if (!is.null(res$convergence)) res$convergence else NA_integer_)
    
    # advance warm-starts
    start_full <- full_next
    if (!is.null(res$invhessian.lt)) invh_lt <- res$invhessian.lt
    
    # lightweight text progress (no progress bar)
    if (verbose && (nsteps %% every == 0L || cur_ll <= thresh)) {
      pct <- (cur_ll - optim_ll) / (thresh - optim_ll)
      pct <- max(0, min(1, pct))
      message(sprintf(
        "Progress: %d%% (iter %d, conv=%s)",
        round(pct * 100),
        nsteps,
        if (is.null(res$convergence)) "NA" else as.character(res$convergence)
      ))
      if (cur_ll <= thresh) message("Reached likelihood ratio threshold on this side.")
    }
  }
  
  if (verbose && nsteps >= max_steps && cur_ll > thresh) {
    side_label <- if (direction < 0) "left" else "right"
    message(sprintf("%s side: reached maximum number of steps without crossing threshold", side_label))
  }
  
  list(
    ll        = ll_vec,
    vals      = val_vec,
    fulls     = full_lst,
    conv      = conv_vec,
    last_full = start_full,
    last_invh = invh_lt,
    steps     = nsteps,
    crossed   = (cur_ll <= thresh)
  )
}

#' Basic (non-adaptive) tool for profiling the likelihood
#'
#' @param profile_parameter Character. Name of the parameter to profile. Profiles
#' are done on the math scale.
#' @param increment_left Numeric. Step size (math scale) when moving to the left,
#' from the start point of the parameter point estimate, to construct the profile.
#' @param increment_right Numeric. Step size (math scale) when moving to the right,
#' from the start point of the paraneter point estimate, to construct the profile.
#' @param num_steps_left Integer. Maximum number of steps to take to the left.
#' @param num_steps_right Integer. Maximum number of steps to take to the right.
#' @param alpha Numeric value between 0 and 1. Confidence level used for the
#' likelihood ratio (LR) threshold: threshold = MLE_loglik - qchisq(alpha, 1)/2.
#' @param optim_param_vector Named numeric. MLE parameters on math scale.
#' @param env_dat 3D array (locations x time x variables).
#' @param occ Logical, either 0 or 1, vector (length = number of locations).
#' @param mask Named numeric or NULL. Parameters kept fixed (math scale).
#' @param num_threads Integer. Threads used internally by log-likelihood.
#' @param control Named list. Control passed to
#' \code{ucminfcpp::ucminf_xptr(control = ...)}.
#'   User-specified entries override defaults:
#'   \itemize{
#'     \item \code{grad = "central"}
#'     \item \code{gradstep = c(1e-6, 1e-8)}
#'     \item \code{grtol = 1e-5}
#'     \item \code{xtol = 1e-12}
#'     \item \code{stepmax = 5}
#'     \item \code{maxeval = 2000}
#'   }
#'   If you want optimizer iteration trace, set \code{control$trace > 0}.
#' @param verbose Logical. If \code{TRUE}, prints compact progress messages; otherwise silent.
#'
#' @returns A list with:
#' \itemize{
#'   \item \code{profile}: data.frame with columns \code{param}, \code{value_math},
#'         \code{loglik}, \code{convergence}, and a list-column \code{full_par}.
#'         (No side/step columns.)
#'   \item \code{found_better}: logical; TRUE if any profiled point exceeds the MLE log-likelihood.
#'   \item \code{threshold}: numeric; LR threshold used.
#'   \item \code{parameters}: data.frame with the parameters found in each step
#'   of the profiling
#' }
#' @export
#'
#' @examples
#' ## Minimal profiling example (fast): 1 step left + 1 step right
#' res <- profile_likelihood(
#'   profile_parameter = "mu1",
#'  increment_left = 0.2,
#'  increment_right = 0.2,
#'  num_steps_left = 1L, # one iteration on the left
#'  num_steps_right = 1L, # one iteration on the right
#'  alpha = 0.95,
#'  optim_param_vector = example_1$optim_par_vec,
#'  env_dat = example_1$env_array,
#'  occ = example_1$occ_df$presence,
#'  num_threads = 1L, # keep it fast and deterministic
#'  control = list(maxeval = 20),
#'  verbose = FALSE
#')
#' # Check the structure of the output:
#' res$profile
#' res$threshold
#' res$found_better
#' ## Full math-scale parameter vectors used at each evaluated point:
#' res$parameter_df
profile_likelihood <- function(
    profile_parameter = "mu1",
    increment_left = 0.1,
    increment_right = increment_left,
    num_steps_left = 20L,
    num_steps_right = num_steps_left,
    alpha = 0.95,
    optim_param_vector,
    env_dat,
    occ,
    mask = NULL,
    num_threads = RcppParallel::defaultNumThreads(),
    control = list(),
    verbose = FALSE
) {
  # --------------------------
  # Parameter validation (checkmate)
  # --------------------------
  # env_dat: 3D array with no NAs
  check_env_array(env_dat)
  
  # occ: logical (no NAs) OR integerish in {0,1} (no NAs); and length matches locations
  n_loc <- dim(env_dat)[1]
  checkmate::assert(
    checkmate::check_logical(occ, any.missing = FALSE, len = n_loc),
    checkmate::check_integerish(occ, lower = 0, upper = 1, any.missing = FALSE, len = n_loc),
    .var.name = "occ"
  )
  
  # optim_param_vector: named numeric, finite, no NAs, unique names, non-empty
  checkmate::assert_numeric(optim_param_vector, any.missing = FALSE, finite = TRUE, min.len = 1)
  checkmate::assert_character(names(optim_param_vector), any.missing = FALSE, min.len = length(optim_param_vector), unique = TRUE)
  
  
  # profile_parameter: single, non-missing character; must be in names(optim_param_vector)
  checkmate::assert_character(profile_parameter, any.missing = FALSE, len = 1)
  checkmate::assert_choice(profile_parameter, choices = names(optim_param_vector))
  
  # mask: NULL or named numeric; names subset of optim_param_vector names; no NAs; unique names
  checkmate::assert_numeric(mask, any.missing = FALSE, null.ok = TRUE)
  # if (!is.null(mask)) {
  #   checkmate::assert_character(names(mask), any.missing = FALSE, min.len = length(mask), unique = TRUE)
  #   checkmate::assert_subset(names(mask), choices = names(optim_param_vector))
  #   # avoid conflicts: profiled parameter must not be pre-fixed by user mask
  #   if (profile_parameter %in% names(mask)) {
  #     stop("`mask` must not contain the `profile_parameter` you are profiling; remove it from `mask`.", call. = FALSE)
  #   }
  # }
  
  
  # increments: numeric scalars > 0 (finite)
  
  
  checkmate::assert_number(increment_left, lower = 0, finite = TRUE)
  checkmate::assert_true(increment_left > 0)
  
  checkmate::assert_number(increment_right, lower = 0, finite = TRUE)
  checkmate::assert_true(increment_right > 0)
  
  # alpha in (0,1)
  checkmate::assert_number(alpha, lower = 0, upper = 1, finite = TRUE)
  checkmate::assert_true(alpha > 0 && alpha < 1)
  
  
  # step caps: integerish scalars >= 1masked_now
  checkmate::assert_integerish(num_steps_left, lower = 1, any.missing = FALSE, len = 1)
  checkmate::assert_integerish(num_steps_right, lower = 1, any.missing = FALSE, len = 1)
  
  
  # num_threads: integerish scalar >= 1
  checkmate::assert_integerish(num_threads, lower = 1, any.missing = FALSE, len = 1)
  
  # control: list (content validated by ucminfcpp)
  checkmate::assert_list(control, any.missing = FALSE, null.ok = TRUE)
  
  # Ensure there is at least one free parameter to optimize at the MLE point
  masked_now <- unique(c(profile_parameter, if (!is.null(mask)) names(mask) else character(0)))
  n_free_now <- length(setdiff(names(optim_param_vector), masked_now))
  if (n_free_now < 1L) {
    stop("No free parameters left to optimize after fixing the profile parameter and the `mask`. ",
         "Either remove some names from `mask` or profile a different parameter.",
         call. = FALSE
    )
  }
  
  
  if (!is.null(mask)) {
    overlap <- intersect(names(optim_param_vector), names(mask))
    if (length(overlap) > 0) {
      stop(
        "`optim_param_vector` (free parameters) and `mask` (fixed parameters) ",
        "must be disjoint. Overlapping names: ", paste(overlap, collapse = ", "), ".\n",
        "Pass the optimizer output vector (free parameters only) together with the same `mask` ",
        "used in that optimization.",
        call. = FALSE
      )
    }
  }
  
  
  # --------------------------
  # Merge user control with defaults (user wins)
  # --------------------------
  default_ctrl <- list(
    grad     = "central",
    gradstep = c(1e-6, 1e-8),
    grtol    = 1e-5,
    xtol     = 1e-12,
    stepmax  = 5,
    maxeval  = 1000
  )
  base_control <- utils::modifyList(default_ctrl, control)
  
  # --------------------------
  # MLE log-likelihood at the provided optimum
  # --------------------------
  optim_ll <- loglik_math(
    param_vector = optim_param_vector,
    env_dat = env_dat,
    occ = occ,
    negative = FALSE,
    mask = mask,
    num_threads = num_threads
  )
  
  # Likelihood-ratio threshold (df = 1)
  thresh <- optim_ll - stats::qchisq(alpha, df = 1L) / 2
  
  # --------------------------
  # Seed storages with the MLE row
  # --------------------------
  storage_ll <- list(optim_ll)
  storage_val <- list(unname(optim_param_vector[profile_parameter]))
  storage_full <- list(optim_param_vector)
  storage_conv <- list(NA_integer_)
  
  # --------------------------
  # LEFT
  # --------------------------
  left <- profile_one_side_(
    direction = -1,
    increment = increment_left,
    max_steps = num_steps_left,
    profile_parameter = profile_parameter,
    optim_param_vector = optim_param_vector,
    env_dat = env_dat, occ = occ,
    mask = mask,
    num_threads = num_threads,
    optim_ll = optim_ll, thresh = thresh,
    base_control = base_control,
    start_full = optim_param_vector,
    invh_lt = NULL,
    verbose = verbose
  )
  
  storage_ll <- c(storage_ll, left$ll)
  storage_val <- c(storage_val, left$vals)
  storage_full <- c(storage_full, left$fulls)
  storage_conv <- c(storage_conv, left$conv)
  
  # --------------------------
  # RIGHT
  # --------------------------
  right <- profile_one_side_(
    direction = +1,
    increment = increment_right,
    max_steps = num_steps_right,
    profile_parameter = profile_parameter,
    optim_param_vector = optim_param_vector,
    env_dat = env_dat, occ = occ, mask = mask,
    num_threads = num_threads,
    optim_ll = optim_ll, thresh = thresh,
    base_control = base_control,
    start_full = optim_param_vector,
    invh_lt = NULL,
    verbose = verbose
  )
  
  storage_ll <- c(storage_ll, right$ll)
  storage_val <- c(storage_val, right$vals)
  storage_full <- c(storage_full, right$fulls)
  storage_conv <- c(storage_conv, right$conv)
  
  # --------------------------
  # Assemble output (no side/step columns)
  # --------------------------
  profile_df <- data.frame(
    param = profile_parameter,
    value_math = unlist(storage_val, use.names = FALSE),
    loglik = unlist(storage_ll, use.names = FALSE),
    convergence = unlist(storage_conv, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  # --- Adding mask to the optim_vector
  # --- Build storage_full with same lenght
  if (!is.null(mask)) {
    storage_full[[1]] <- c(storage_full[[1]], mask)
  }
  # --- Sort by value_math (increasing), and apply the same order to storage_full ---
  ord <- order(profile_df$value_math, na.last = TRUE)
  profile_df <- profile_df[ord, , drop = FALSE]
  storage_full <- storage_full[ord]
  
  
  # --- Convert storage_full (list of named numeric vectors) to a data.frame ---
  parameter_df <- do.call(
    rbind,
    lapply(storage_full, function(v) {
      out <- as.data.frame(t(v), stringsAsFactors = FALSE)
      rownames(out) <- NULL
      out
    })
  )
  found_better <- any(profile_df$loglik > optim_ll + .Machine$double.eps)
  
  list(
    profile      = profile_df,
    found_better = found_better,
    threshold    = thresh,
    parameters   = parameter_df
  )
}
