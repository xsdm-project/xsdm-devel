#' Helper (internal): automatically derive relative plot limits for
#' `interpret_parameters()`
#'
#' Computes variable-specific plotting limits **for internal use only**.
#'
#' ## Details
#'
#' Limits are expressed *relative to* the species' optimal environmental
#' values `mu`, i.e. each returned limit is reported as
#' `(abs_limit - mu_i)`.
#'
#' Width is controlled by a single scalar `breadth`:
#'
#' - `breadth = 1` (default) gives the widest view: the full min–max range
#'   plus a symmetric margin of `margin * range` on each side.
#' - `breadth = 0` collapses every limit pair to essentially a single point
#'   around `mu`.
#' - Values in between interpolate linearly.
#'
#' @param env_dat A 3D numeric array. **Must not contain missing values.**
#' @param param_list Named list, must contain entry `mu`.
#' @param indices Integer vector specifying which environmental variables
#'   to compute limits for.
#' @param breadth Scalar in \eqn{[0, 1]}. Default `1`.
#' @param margin Non-negative scalar. Fraction of observed range to expand.
#' Default is 0.1, i.e. 10% of the observed range on each side.
#'
#' @returns A list of length `length(indices)`. Each element is a numeric
#'   length-2 vector `c(lower, upper)` of limits *relative to* `mu`.
#'
#' @noRd
#' @keywords internal
auto_plot_lims_ <- function(env_dat,
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


#' Tool to help interpret xsdm model parameters
#'
#' Due to the parameter reduction step which was carried out to eliminate
#' structural non-identifiability in the xsdm model, parameter interpretation
#' is more difficult. This function helps with that difficulty, displaying
#' contours for the inferred log growth-environment function. The shapes of
#' these contours are determined by inference, even though their levels are
#' not; and the shapes are generally more informative anyway. See the manual
#' documents ``The xsdm model'' and ``How to fit xsdm models with species
#' occurrence data using xsdm'' for additional details.
#'
#' If \code{env_dat} and \code{occ} are provided, two panels are drawn side
#' by side: on the left the growth-environment function is shown together
#' with the environmental values at presence locations (\code{occ == 1});
#' on the right the same function is shown together with the environmental
#' values at non-detections (\code{occ == 0}). Both panels share identical
#' contour breaks (bivariate case) or identical axes (univariate case), so
#' the two are directly comparable.
#'
#' @param param_list A named list of xsdm model parameters such as returned
#'   by \code{math_to_bio}. Must contain elements \code{mu}, \code{sigltil},
#'   \code{sigrtil}, \code{ctil}, \code{pd}, and \code{o_mat}.
#' @param plot_indices A length-1 or length-2 integer vector of indices of
#'   environmental variables against which the growth-environment function
#'   is to be plotted. For a length-2 vector, the first index is the
#'   horizontal axis, the second the vertical. Other environmental
#'   variables are held at their values in \code{param_list$mu}.
#' @param plot_lims Optional list of the same length as \code{plot_indices},
#'   each element a 2-vector giving the plotting extent *relative to*
#'   \code{mu}. If \code{NULL} (the default) and \code{env_dat} is
#'   supplied, limits are auto-derived via \code{auto_plot_lims_()} using
#'   the \code{breadth} argument. The auto-derived limits cover the
#'   full observed environmental range plus a symmetric margin on each side.
#' @param env_dat Optional 3D numeric array of environmental data with
#'   dimensions \code{(locations) x (time) x (variables)}. Required for the
#'   two-panel (presence vs non-detection) display and for auto-derived
#'   \code{plot_lims}. If \code{NULL}, a single-panel legacy plot is drawn
#'   and \code{plot_lims} must be supplied.
#' @param occ Optional length-\code{(locations)} logical or 0/1 vector of
#'   presence/absence. Required together with \code{env_dat} for the
#'   two-panel display.
#' @param breadth Scalar in \code{[0, 1]} controlling how wide the
#'   auto-derived plotting window is around \code{mu}:
#'   \code{breadth = 1} (default) shows
#'   the full min-max environmental range plus a 10\% margin on each side;
#'   \code{breadth = 0} collapses to essentially a single point.
#'   Ignored when \code{plot_lims} is supplied.
#' @param ... Additional graphical arguments passed to \code{plot} (1D case)
#'   or \code{image} (2D case).
#'
#' @returns Invisibly returns the (possibly auto-derived) \code{plot_lims}
#'   list, so downstream code can reuse the same limits. The main purpose
#'   of the function is its side effect: plots are sent to the default
#'   graphics device.
#'
#' @details
#' The log growth-environment function is determined by inference only up
#' to an affine transformation \eqn{g = a f(e) + b} with \eqn{a > 0}. Its
#' contours are therefore unlabelled in the output; their shape is what is
#' interpretively meaningful. In code the function is
#' \deqn{y(e) = -\sum_i \left( \frac{[u_i]_+}{\sigma^R_i}
#'                            + \frac{[u_i]_-}{\sigma^L_i} \right)^2 ,
#'       \quad u = O^{T} (e - \mu),}
#' which is always \eqn{\le 0}, attains its maximum 0 at \eqn{e = \mu},
#' and decreases without bound as \eqn{e} moves away from \eqn{\mu}.
#' Consequently the numeric values on the y-axis of the univariate plot
#' and the numeric values of the image colors in the bivariate plot carry
#' no units of their own.
#'
#' @importFrom graphics abline image mtext par points
#' @importFrom grDevices hcl.colors adjustcolor
#' @export
#'
#' @examples
#' \dontrun{
#'   # Two-panel (presence vs non-detection) plot with auto-derived limits
#'   interpret_parameters(
#'     examples$par_list,
#'     plot_indices = c(1, 2),
#'     env_dat      = examples$env_array,
#'     occ          = examples$occ_vec
#'   )
#'
#'   # Narrower auto-derived window
#'   interpret_parameters(
#'     examples$par_list,
#'     plot_indices = c(1, 2),
#'     env_dat      = examples$env_array,
#'     occ          = examples$occ_vec,
#'     breadth      = 0.7
#'   )
#'
#' }
interpret_parameters <- function(param_list,
                                 plot_indices,
                                 plot_lims = NULL,
                                 env_dat   = NULL,
                                 occ       = NULL,
                                 breadth   = 1,
                                 ...) {
  
  # ---- Validate param_list ------------------------------------------------
  checkmate::assert_true(is.list(param_list))
  checkmate::assert_true(setequal(
    names(param_list),
    c("mu", "sigltil", "sigrtil", "ctil", "pd", "o_mat")
  ))
  checkmate::assert_numeric(param_list$mu,
                            finite = TRUE, any.missing = FALSE
  )
  checkmate::assert_numeric(param_list$sigltil,
                            lower = 0, any.missing = FALSE
  )
  checkmate::assert_numeric(param_list$sigrtil,
                            lower = 0, any.missing = FALSE
  )
  checkmate::assert_numeric(param_list$o_mat,
                            finite = TRUE, any.missing = FALSE
  )
  
  p <- length(param_list$mu)
  checkmate::assert_true(length(param_list$sigltil) == p)
  checkmate::assert_true(length(param_list$sigrtil) == p)
  checkmate::assert_true(all(dim(param_list$o_mat) == c(p, p)))
  
  # ---- Validate plot_indices ---------------------------------------------
  checkmate::assert_numeric(plot_indices,
                            finite = TRUE, any.missing = FALSE
  )
  checkmate::assert_true(all(plot_indices %in% seq_len(p)))
  checkmate::assert_true(length(plot_indices) %in% c(1L, 2L))
  
  # ---- Validate env_dat / occ consistency --------------------------------
  have_data <- !is.null(env_dat) && !is.null(occ)
  if (!is.null(occ) && is.null(env_dat)) {
    stop("`occ` was supplied without `env_dat`.", call. = FALSE)
  }
  if (have_data) {
    check_env_array(env_dat)
    checkmate::assert(
      checkmate::check_logical(occ,
                               any.missing = FALSE, len = dim(env_dat)[1]
      ),
      checkmate::check_integerish(occ,
                                  lower = 0, upper = 1, any.missing = FALSE,
                                  len = dim(env_dat)[1]
      )
    )
    occ <- as.integer(occ)
  }
  
  checkmate::assert_number(
    breadth,
    lower = 0, upper = 1, finite = TRUE, na.ok = FALSE
  )
  
  # ---- Resolve plot_lims -------------------------------------------------
  if (is.null(plot_lims)) {
    if (is.null(env_dat)) {
      stop(
        "`plot_lims` was not supplied and cannot be auto-derived because ",
        "`env_dat` is also missing. Provide either `plot_lims` or ",
        "`env_dat`.",
        call. = FALSE
      )
    }
    plot_lims <- auto_plot_lims_(
      env_dat    = env_dat,
      param_list = param_list,
      indices    = plot_indices,
      breadth    = breadth
    )
  }
  
  checkmate::assert_true(is.list(plot_lims))
  checkmate::assert_true(length(plot_lims) == length(plot_indices))
  for (i in seq_along(plot_indices)) {
    checkmate::assert_numeric(plot_lims[[i]],
                              finite = TRUE, any.missing = FALSE, len = 2L
    )
    checkmate::assert_true(plot_lims[[i]][1] < plot_lims[[i]][2])
  }
  
  # ---- Internal helper: log growth-environment function ------------------
  compute_y <- function(e_m_mu) {
    n <- ncol(e_m_mu)
    u <- t(param_list$o_mat) %*% e_m_mu
    inv_r <- matrix(1 / param_list$sigrtil, nrow = p, ncol = n)
    inv_l <- matrix(1 / param_list$sigltil, nrow = p, ncol = n)
    -colSums((pmax(u, 0) * inv_r + pmin(u, 0) * inv_l)^2)
  }
  
  # ---- Flatten (location x time) env values for scatter overlays ---------
  if (have_data) {
    n_time   <- dim(env_dat)[2]
    occ_flat <- rep(occ, times = n_time)
    env_flat <- vector("list", length(plot_indices))
    for (k in seq_along(plot_indices)) {
      env_flat[[k]] <- as.numeric(env_dat[, , plot_indices[k]])
    }
  }
  
  # ---- Univariate branch --------------------------------------------------
  if (length(plot_indices) == 1L) {
    idx <- plot_indices[1L]
    
    len  <- 500L
    x    <- seq(
      from = param_list$mu[idx] + plot_lims[[1]][1],
      to   = param_list$mu[idx] + plot_lims[[1]][2],
      length.out = len
    )
    e_m_mu <- matrix(0, nrow = p, ncol = len)
    e_m_mu[idx, ] <- x - param_list$mu[idx]
    y <- compute_y(e_m_mu)
    
    if (have_data) {
      op <- par(mfrow = c(1L, 2L))
      on.exit(par(op), add = TRUE)
      plot_one_1d <- function(title_tag, keep, ...) {
        plot(x, y,
             type = "l",
             xlab = paste("Environmental variable", idx),
             ylab = "", yaxt = "n",
             main = title_tag,
             ...
        )
        mtext("Log growth-environment function", side = 2, line = 1.1)
        abline(v = param_list$mu[idx], lty = "dashed")
        ymin  <- min(y)
        ypad  <- 0.10 * diff(range(y))
        env_k <- env_flat[[1]][keep]
        if (length(env_k)) {
          points(
            env_k,
            stats::runif(length(env_k), ymin - ypad, ymin - 0.2 * ypad),
            pch = 16, cex = 0.5,
            col = grDevices::adjustcolor("black", alpha.f = 0.3)
          )
        }
      }
      plot_one_1d("Presences (occ == 1)",      occ_flat == 1L, ...)
      plot_one_1d("Non-detections (occ == 0)", occ_flat == 0L, ...)
    } else {
      plot(x, y,
           type = "l",
           xlab = paste("Environmental variable", idx),
           ylab = "", yaxt = "n",
           ...
      )
      mtext("Log growth-environment function", side = 2, line = 1.1)
      abline(v = param_list$mu[idx], lty = "dashed")
    }
  }
  
  # ---- Bivariate branch ---------------------------------------------------
  if (length(plot_indices) == 2L) {
    len  <- 500L
    idx1 <- plot_indices[1L]
    idx2 <- plot_indices[2L]
    
    x_1 <- seq(
      from = param_list$mu[idx1] + plot_lims[[1]][1],
      to   = param_list$mu[idx1] + plot_lims[[1]][2],
      length.out = len
    )
    x_2 <- seq(
      from = param_list$mu[idx2] + plot_lims[[2]][1],
      to   = param_list$mu[idx2] + plot_lims[[2]][2],
      length.out = len
    )
    grid   <- t(expand.grid(x_1, x_2))
    e_m_mu <- matrix(0, nrow = p, ncol = len^2)
    e_m_mu[plot_indices, ] <- grid - matrix(
      param_list$mu[plot_indices], nrow = 2L, ncol = len^2
    )
    z <- matrix(compute_y(e_m_mu), nrow = len, ncol = len)
    
    dots <- list(...)
    if (is.null(dots$zlim)) dots$zlim <- range(z, finite = TRUE)
    if (is.null(dots$col))  dots$col  <- grDevices::hcl.colors(64, "YlOrRd",
                                                               rev = TRUE)
    
    draw_one_2d <- function(title_tag, keep, ...) {
      do.call(image, c(
        list(x = x_1, y = x_2, z = z, axes = TRUE,
             xlab = paste("Environmental variable", idx1),
             ylab = paste("Environmental variable", idx2),
             main = title_tag),
        dots,
        list(...)
      ))
      points(
        param_list$mu[idx1], param_list$mu[idx2],
        pch = 16, col = "green", cex = 1
      )
      if (have_data) {
        xk <- env_flat[[1]][keep]
        yk <- env_flat[[2]][keep]
        if (length(xk)) {
          points(
            xk, yk,
            pch = 16, cex = 0.4,
            col = grDevices::adjustcolor("black", alpha.f = 0.4)
          )
        }
      }
    }
    
    if (have_data) {
      op <- par(mfrow = c(1L, 2L))
      on.exit(par(op), add = TRUE)
      draw_one_2d("Presences (occ == 1)",      occ_flat == 1L, ...)
      draw_one_2d("Non-detections (occ == 0)", occ_flat == 0L, ...)
    } else {
      draw_one_2d("", logical(0), ...)
    }
  }
  
  invisible(plot_lims)
}