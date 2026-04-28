#' Helper. Reasonable ranges for initial conditions for optimizations seeking to
#' maximize the likelihood of the xsdm model
#'
#' Given environmental data, constructs reasonable ranges for each xsdm
#' model parameter within which to select initial conditions for the
#' optimizer.
#'
#' @param env_dat A 3D numeric array of environmental time series data
#'   with dimensions \code{(locations) x (time) x (environmental variables)}.
#'   Missing values are not allowed. Assumed to pertain only to locations
#'   where the species was observed.
#' @param breadth Scalar in \code{[0, 1]} controlling how wide the search
#'   ranges are around their (fixed, data-driven) center. \code{breadth = 1}
#'   (the default) reproduces the pre-v0.3 behaviour that corresponded to
#'   \code{quant_vec = c(0.1, 0.5, 0.9)}; \code{breadth = 0} collapses
#'   every range to essentially a single point, equivalent to
#'   \code{quant_vec = c(0.5 - 1e-6, 0.5, 0.5 + 1e-6)}. Values in between
#'   interpolate linearly.
#' @param quant_vec Deprecated. If supplied, a deprecation warning is
#'   emitted and \code{breadth} is set to
#'   \code{(quant_vec[3] - quant_vec[1]) / 0.8}. This argument will be
#'   removed in a future release.
#'
#' @returns A data.frame with three columns (lower bound, center, upper
#'   bound) giving the search range for each parameter.
#'
#' @details The center of every range is fixed and data-driven: it is
#'   the empirical median for \code{mu}-type parameters, \code{0} for
#'   orthogonal-matrix angles, \code{0.5} for \code{pd}, and the median
#'   of the evaluated log-likelihood surface for \code{ctil}. Only the
#'   half-width around that center is user-controllable, via
#'   \code{breadth}.
#'
#'   Internally, \code{breadth} is mapped to two monotonic quantities:
#'   \itemize{
#'     \item \code{half_width = 1e-6 + breadth * (0.4 - 1e-6)} — the
#'           half-width in probability space used for quantiles of
#'           \code{mu}, \code{pd}, and \code{ctil}. Quantile arguments
#'           are clamped so that \code{logit} stays finite.
#'     \item \code{fact = 1 + breadth} — the multiplicative factor on
#'           the log scale for \code{sigltil} and \code{sigrtil} ranges.
#'           \code{breadth = 1} gives \code{fact = 2} (the legacy value);
#'           \code{breadth = 0} gives \code{fact = 1}, i.e. a
#'           degenerate single-point range.
#'   }
#' @keywords internal
#' @examples
#' set.seed(1)
#' env <- array(rnorm(10 * 5 * 2), dim = c(10, 5, 2))
#' xsdm:::get_range_df_(env)                 # default breadth = 1
#' xsdm:::get_range_df_(env, breadth = 0.3) 
get_range_df_ <- function(env_dat,
                         breadth = 1,
                         quant_vec = NULL) {
  # Validate env_dat
  check_env_array(env_dat)
  
  # --- Soft deprecation of the old quant_vec argument ---------------------
  if (!is.null(quant_vec)) {
    .Deprecated(
      msg = paste0(
        "`quant_vec` is deprecated in `get_range_df_()`; use the scalar ",
        "`breadth` argument instead (0 <= breadth <= 1, default 1 = old ",
        "quant_vec = c(0.1, 0.5, 0.9))."
      )
    )
    checkmate::assert_numeric(
      quant_vec,
      any.missing = FALSE,
      finite = TRUE,
      len = 3
    )
    checkmate::assert(
      all(quant_vec > 0 & quant_vec < 1),
      msg = "All values in quant_vec must be strictly between 0 and 1 (0 < x < 1)."
    )
    checkmate::assert(
      if (all(diff(quant_vec) > 0)) TRUE else "not strictly increasing",
      msg = "quant_vec must be strictly increasing (no ties)."
    )
    # Inverse of the new symmetric mapping. Asymmetry in the supplied
    # quant_vec is ignored on purpose: the new parametrization is
    # symmetric around 0.5 by design.
    breadth <- (quant_vec[3] - quant_vec[1]) / 0.8
  }
  
  # Validate breadth
  checkmate::assert_number(
    breadth,
    lower = 0,
    upper = 1,
    finite = TRUE,
    na.ok = FALSE
  )
  
  # --- Derived quantities used throughout the function --------------------
  half_width <- 1e-6 + breadth * (0.4 - 1e-6)
  probs <- c(0.5 - half_width, 0.5, 0.5 + half_width)
  fact <- 1 + breadth
  
  # Number of environmental variables
  p <- dim(env_dat)[3]
  
  # Receptacle for results
  ranges <- data.frame(
    lower = NA * numeric(num_par(p)),
    center = NA * numeric(num_par(p)),
    upper = NA * numeric(num_par(p))
  )
  
  # Number of parameters needed for the orthogonal matrix
  q <- (p^2 - p) / 2
  
  # Name parameters in the receptacle for results
  if (q != 0) {
    o_inds <- seq_len(q)
    rownames(ranges)[o_inds] <- paste0("o_par", 1:q)
  }
  mu_inds <- (1 + q):(q + p)
  rownames(ranges)[mu_inds] <- paste0("mu", 1:p)
  sigl_inds <- (1 + q + p):(q + 2 * p)
  rownames(ranges)[sigl_inds] <- paste0("sigltil", 1:p)
  sigr_inds <- (1 + q + 2 * p):(q + 3 * p)
  rownames(ranges)[sigr_inds] <- paste0("sigrtil", 1:p)
  pd_inds <- (1 + q + 3 * p):(num_par(p) - 1)
  rownames(ranges)[pd_inds] <- "pd"
  ctil_inds <- num_par(p)
  rownames(ranges)[ctil_inds] <- "ctil"
  
  # Infer a reasonable range for start guesses for mu, sigltil and sigrtil
  for (counter in 1:p) {
    h <- as.numeric(env_dat[, , counter])
    
    # get mu range
    ranges[mu_inds[counter], ] <- unname(stats::quantile(h, probs = probs))
    
    # get sigL range
    mu_center <- ranges[mu_inds[counter], 2]
    h2 <- sqrt(mean((h[h < mu_center] - mu_center)^2))
    # the log is because we want ranges on the math scale
    ranges[sigl_inds[counter], ] <- log(c(h2 / fact, h2, fact * h2))
    
    # get sig_r range
    h2 <- sqrt(mean((h[h > mu_center] - mu_center)^2))
    ranges[sigr_inds[counter], ] <- log(c(h2 / fact, h2, fact * h2))
  }
  
  # Start range for orthogonal matrix parameters intended to blanket the
  # space of (special) orthogonal matrices
  if (q != 0) {
    ranges[o_inds, 1] <- -3 * pi
    ranges[o_inds, 2] <- 0
    ranges[o_inds, 3] <- 3 * pi
    o_star <- build_orthogonal_matrix(ranges[o_inds, 2])
  } else if (q == 0) {
    o_star <- build_orthogonal_matrix(entries = NULL)
  }
  
  # pd parameter
  ranges[pd_inds, ] <- logit(probs)
  
  # ctil range
  #
  # 1. Evaluate loglik in  central values of all the other parameters
  # 2. Pick the central ctil to be the opposite of the median of those values
  # 2.1 For those central parameters half the probabilities of detection
  # will be < 0.5 and the other half will be > 0.5,
  # 3. Pick the lower and upper ctil based on quantiles of the distribution
  mu_star <- ranges[mu_inds, 2]
  sigl_star <- exp(ranges[sigl_inds, 2])
  sigr_star <- exp(ranges[sigr_inds, 2])
  h_star <- like_neg_ltsgr(
    env_dat = env_dat,
    mu = mu_star,
    sigltil = sigl_star,
    sigrtil = sigr_star,
    o_mat = o_star
  )
  ranges[ctil_inds, ] <- -stats::quantile(h_star, rev(probs))
  if (q != 0) {
    ranges <- ranges[c(
      mu_inds, sigl_inds, sigr_inds,
      ctil_inds, pd_inds, o_inds
    ), ]
  }
  
  return(ranges)
}

#' get_start_parms_. Generates a Latin hypercube design for the parameters
#' based on the Sobol' low-discrepancy sequence.
#' Given a set of ranges of environmental variables create a sample
#' of parameters.
#'
#' @param ranges A data frame with ranges to generate the parameter hypercube
#' of parameters
#' @param numstarts The number of require samples
#'
#' @keywords internal
get_start_parms_ <- function(ranges, numstarts = 100) {
  if (!requireNamespace("pomp", quietly = TRUE)) {
    stop("Package 'pomp' is required for this function. Install it with install.packages('pomp').")
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required for this function. Install it with install.packages('tibble').")
  }
  # Check parameters. numstarts must be >= 3: we call
  # pomp::sobol_design(..., nseq = numstarts - 1) and then append the
  # data-driven center row. nseq = 0 segfaults pomp, nseq = 1 returns a
  # vector instead of a matrix and downstream colnames<- fails.
  checkmate::assert_data_frame(ranges, any.missing = FALSE, ncols = 3)
  checkmate::assert_integerish(numstarts, lower = 3, any.missing = FALSE,
                               len = 1,
                               .var.name = "numstarts")
  checkmate::assert_names(names(ranges),
                          must.include = c("lower", "center", "upper")
  )
  # Set limits on the parameters
  lower <- ranges[, 1]
  names(lower) <- rownames(ranges)
  center <- ranges[, 2]
  names(center) <- rownames(ranges)
  upper <- ranges[, 3]
  names(upper) <- rownames(ranges)
  
  # Get the actual start parameters, math scale
  if (!requireNamespace("pomp", quietly = TRUE)) {
    stop("Package 'pomp' is required for get_start_parms_(). Install it with: install.packages('pomp')")
  }
  startparms_math <- pomp::sobol_design(
    lower = lower,
    upper = upper,
    nseq = numstarts - 1
  )
  # Defensive coercion: pomp::sobol_design has historically returned a
  # bare numeric vector (instead of a data.frame) for very small nseq.
  # The validation above already rules out nseq < 2, but the cast keeps
  # us safe against future pomp behavior changes and produces a clearer
  # error if the contract ever drifts.
  if (!is.data.frame(startparms_math) && !is.matrix(startparms_math)) {
    startparms_math <- as.data.frame(
      matrix(startparms_math, nrow = numstarts - 1L,
             dimnames = list(NULL, names(lower)))
    )
  }
  startparms_math <- rbind(startparms_math, center)
  tibble::as_tibble(startparms_math)
}


#' Starting parameters for the optimization
#'
#' Generate starting parameters for the optimization of the xsdm log-likelihood.
#' Starting points are constructed from the environmental conditions at observed
#' presences (where \code{occ == 1}) using a Latin hypercube design for the
#' parameters based on the Sobol' low-discrepancy sequence.
#'
#' @param env_dat The environmental array for only the observed occurrences
#' @param mask Either NULL or a named numeric vector. Names must be as specified
#' by calling \code{make_mask_names}. The NULL case means planned optimizations
#' will be over all model parameters, so start parameter sets should include all
#' model parameters. The non-NULL case means some parameters will not be
#' specified in the output of this function because the optimizations which
#' are planned will fix those parameters anyway. The most common case for most
#' applications will be \code{mask=NULL}.
#' @param breadth Scalar in \code{[0, 1]} controlling how wide the search
#'   ranges are around their (fixed, data-driven) center. \code{breadth = 1}
#'   (the default) reproduces the pre-v0.3 behaviour that corresponded to
#'   \code{quant_vec = c(0.1, 0.5, 0.9)}; \code{breadth = 0} collapses
#'   every range to essentially a single point, equivalent to
#'   \code{quant_vec = c(0.5 - 1e-6, 0.5, 0.5 + 1e-6)}. Values in between
#'   interpolate linearly.
#' @param num_starts The number of samples of the hypercube
#'
#' @returns A data frame with samples for each parameter to optimize
#' @export
#'
#' @details The bounds and center of the search range for a mu parameter are
#' based on the quantiles in quant_vec applied to all observations of that
#' environmental variable, over space and time. The relationship between
#' quant_vec and the width of the ranges selected for the other parameters
#' varies, but generally wider ranges in quant_vec produce wider ranges for
#' start parameters.
#'
#' @examples
#' env_dat <- examples$env_array[examples$occ_vec == 1, , ]
#' start_parms(env_dat)
#' start_parms(env_dat, mask = c(mu2 = 5, pd = 1))
start_parms <- function(env_dat,
                        mask = NULL,
                        breadth  = 1 ,
                        num_starts = 100) {
  check_env_array(env_dat)
  checkmate::assert_numeric(breadth ,  null.ok = FALSE, )
  # num_starts >= 3 is required: pomp::sobol_design (called by
  # get_start_parms_) needs nseq = num_starts - 1 >= 2.
  checkmate::assert_integerish(num_starts, lower = 3, any.missing = FALSE,
                               len = 1)
  
  range_df <- get_range_df_(env_dat, breadth)
  
  if (is.null(mask)) {
    out <- get_start_parms_(range_df, num_starts)
  } else {
    diff_names <- setdiff(rownames(range_df), names(mask))
    out <- get_start_parms_(range_df[diff_names, ], num_starts)
  }

  return(out)
}
