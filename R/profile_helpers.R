#' Profile one side of a likelihood profile (internal)
#'
#' Fixes one parameter on the math scale and re-optimizes the remaining
#' parameters along a single direction (left/right) until the LR threshold
#' is reached or a step cap is hit.
#'
#' @inheritParams profile_likelihood
#' @param direction Integer. -1 (left) or +1 (right).
#' @param increment Numeric. Step size on math scale for this side.
#' @param max_steps Integer. Maximum iterations for this side.
#' @param base_control Named list. Control passed to \code{ucminf::ucminf(control = ...)}.
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
.profile_one_side <- function(
  direction, # -1 for left, +1 for right
  increment, # step size (math scale)
  max_steps, # iteration cap
  profile_parameter, # name of the profiled parameter
  optim_param_vector, # full parameter template (named)
  env_dat, occ, mask, # data & fixed params
  num_threads, # threads for loglik
  optim_ll, thresh, # MLE loglik & LR threshold
  base_control, # list passed to ucminf
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

    res <- ucminf::ucminf(
      par = par_init,
      fn = loglik_math,
      env_dat = env_dat,
      mask = new_mask,
      occ = occ,
      negative = TRUE,
      num_threads = num_threads,
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
