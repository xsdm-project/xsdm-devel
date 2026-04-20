#' Automatically derive relative plot limits for `interpret_parameters()`
#'
#' Computes variable-specific plotting limits for use with
#' \code{\link{interpret_parameters}()}, based on the full observed
#' environmental range. Limits are expressed *relative to* the species'
#' optimal environmental values \code{mu}, i.e. each returned limit is
#' reported as \code{(abs_limit - mu_i)}.
#'
#' Width is controlled by a single scalar \code{breadth} with the same
#' semantics used elsewhere in the package (see
#' \code{\link{get_range_df}()}):
#' \itemize{
#'   \item \code{breadth = 1} (default) gives the widest view: the full
#'         min–max range of the environmental data (across *all* locations
#'         and time steps, not filtered by \code{occ}) plus a symmetric
#'         margin of \code{margin * range} on each side. This ensures
#'         that every scatter point is comfortably inside the plotting
#'         window.
#'   \item \code{breadth = 0} collapses every limit pair to essentially
#'         a single point around \code{mu}.
#'   \item Values in between interpolate linearly.
#' }
#'
#' @param env_dat A 3D numeric array of environmental data with dimensions
#'   \code{(locations) x (time) x (variables)}. Must not contain missing
#'   values. The full array is used (no filtering by \code{occ}).
#' @param param_list Named list of parameters such as returned by
#'   \code{\link{math_to_bio}()}. Must contain entry \code{mu}.
#' @param indices Integer vector specifying which environmental variables
#'   (elements of \code{mu}) to compute limits for. Must be a subset of
#'   \code{seq_along(param_list$mu)}.
#' @param breadth Scalar in \code{[0, 1]} controlling the width of the
#'   plotting window around \code{mu}. Default \code{1}.
#' @param margin Non-negative scalar. Fraction of the observed data range
#'   by which to expand the limits symmetrically on each side when
#'   \code{breadth = 1}. Default \code{0.1} (10\% of the data range on
#'   each side). The margin scales linearly with \code{breadth}.
#'
#' @returns A list of length \code{length(indices)}. Each element is a
#'   numeric length-2 vector \code{c(lower, upper)} of limits *relative to*
#'   the corresponding \code{mu[indices[i]]}, suitable for direct use as
#'   the \code{plot_lims} argument of \code{\link{interpret_parameters}()}.
#'
#' @examples
#' \dontrun{
#'   lims <- auto_plot_lims(
#'     env_dat    = example_1_env_array,
#'     param_list = example_1_param_list_example,
#'     indices    = c(1, 2)
#'   )
#'   interpret_parameters(
#'     example_1_param_list_example,
#'     plot_indices = c(1, 2),
#'     plot_lims    = lims,
#'     env_dat      = example_1_env_array,
#'     occ          = example_1_occurrence_vector
#'   )
#' }
#' @export
auto_plot_lims <- function(env_dat,
                           param_list,
                           indices,
                           breadth = 1,
                           margin  = 0.1) {
  # ---- Validation --------------------------------------------------------
  check_env_array(env_dat)
  checkmate::assert_list(param_list, any.missing = FALSE)
  checkmate::assert_true("mu" %in% names(param_list))
  mu <- param_list$mu
  checkmate::assert_numeric(mu, any.missing = FALSE, finite = TRUE)
  p <- length(mu)
  checkmate::assert_integerish(indices,
                               lower = 1, upper = p,
                               any.missing = FALSE
  )
  checkmate::assert_number(
    breadth,
    lower = 0, upper = 1, finite = TRUE, na.ok = FALSE
  )
  checkmate::assert_number(margin, lower = 0, finite = TRUE, na.ok = FALSE)
  
  # ---- Full env range per variable (no occ filter) -----------------------
  env_sub  <- env_dat[, , indices, drop = FALSE]
  abs_lo   <- apply(env_sub, 3, min, na.rm = TRUE)
  abs_hi   <- apply(env_sub, 3, max, na.rm = TRUE)
  
  # Breathing room: fraction of the observed range on each side
  data_range <- abs_hi - abs_lo
  dev        <- margin * data_range
  
  # ---- Relative to mu at breadth = 1 (widest) ---------------------------
  rel_lo_full <- (abs_lo - dev) - mu[indices]
  rel_hi_full <- (abs_hi + dev) - mu[indices]
  
  # ---- Linear interpolation from pinprick (breadth = 0) to full ----------
  eps    <- 1e-6
  rel_lo <- breadth * rel_lo_full + (1 - breadth) * (-eps)
  rel_hi <- breadth * rel_hi_full + (1 - breadth) * ( eps)
  
  Map(c, rel_lo, rel_hi)
}
