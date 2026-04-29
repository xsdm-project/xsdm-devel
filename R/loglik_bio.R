#' Log-likelihood function for the xsdm model, parameters on the biological
#' scale.
#'
#' Computes the log-likelihood for the xsdm model given environmental data, a
#' vector of occurrences and pseudo-absences, and model parameters on the
#' biological scale.
#'
#' This is a thin R wrapper around the C++ implementation \code{loglik_bio_cpp};
#' the optimizer hot path (\code{sum_log_p = TRUE}, \code{return_prob = FALSE})
#' is pure C++. The non-default flag combinations
#' (\code{sum_log_p = FALSE} or \code{return_prob = TRUE}) are computed by
#' delegating to the C++-backed \code{log_prob_detect} and reducing in R.
#' A pure-R reference implementation, \code{loglik_bio_r}, is kept internal
#' to the package and is used only by the parity tests in
#' \code{tests/testthat/test-loglik_bio_r_vs_cpp.R}.
#'
#' @param env_dat The environmental data array, dimensions
#' \code{n_loc x n_time x p} (number of locations x time-series length x number
#' of environmental variables). Must be a 3-dimensional array with no missing
#' values.
#' @param occ Presence/pseudo-absence binary vector. Same length as dimension 1 of
#' \code{env_dat}.
#' @param mu Vector of optimal environmental values. Length \code{p=dim(env_dat)[3]}.
#' Unconstrained real numbers.
#' @param sigltil Vector specifying width of the growth-environment function.
#' Length \code{p=dim(env_dat)[3]}. Positive real numbers, \code{Inf} entries also allowed.
#' @param sigrtil Vector specifying width of the growth-environment function.
#' Length \code{p=dim(env_dat)[3]}. Positive real numbers, \code{Inf} entries also allowed.
#' @param o_mat An orthogonal matrix, dimensions \code{p} by \code{p}.
#' @param ctil Scalar. Relates to the center of the detection-link function.
#' @param pd Maximum probability of detection of the species. Parameter between
#' 0 and 1.
#' @param return_prob Logical (default FALSE). Flag to return likelihood instead
#' of log-likelihood.
#' @param sum_log_p Logical (default TRUE). If FALSE, returns the individual
#' log-likelihoods (or likelihoods, if \code{return_prob} is TRUE) associated with the
#' individual locations, instead of their sum (product, if \code{return_prob} is TRUE).
#' @param num_threads Number of threads for parallel computation. Defaults to
#' \code{RcppParallel::defaultNumThreads()}.
#'
#' @returns A single value, the log-likelihood (or the likelihood, if \code{return_prob}
#' is TRUE); or a vector of location specific values of \code{sum_log_p} is FALSE.
#' @export
#'
#' @examples
#' ll <- loglik_bio(
#'   env_dat = examples$env_array,
#'   occ = examples$occ_vec,
#'   mu = examples$par_list$mu,
#'   sigltil = examples$par_list$sigltil,
#'   sigrtil = examples$par_list$sigrtil,
#'   o_mat = examples$par_list$o_mat,
#'   ctil = examples$par_list$ctil,
#'   pd = examples$par_list$pd
#' )
#' ll
loglik_bio <- function(env_dat,
                       occ,
                       mu,
                       sigltil,
                       sigrtil,
                       o_mat,
                       ctil,
                       pd,
                       return_prob = FALSE,
                       sum_log_p = TRUE,
                       num_threads = RcppParallel::defaultNumThreads()) {
  # Input validation (kept in R so the user-facing error messages and
  # checkmate-style asserts continue to match the historical API surface;
  # the math itself runs entirely in C++ via `loglik_bio_cpp`).
  checkmate::assert(
    checkmate::check_logical(occ, any.missing = FALSE),
    checkmate::check_integerish(occ, lower = 0, upper = 1, any.missing = FALSE),
    .var.name = "occ"
  )
  check_env_array(env_dat)
  checkmate::assert_numeric(mu, any.missing = FALSE, min.len = 1)
  checkmate::assert_numeric(sigltil, any.missing = FALSE, min.len = 1)
  checkmate::assert_numeric(sigrtil, any.missing = FALSE, min.len = 1)
  checkmate::assert_numeric(ctil, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(pd, any.missing = FALSE, len = 1)
  checkmate::assert_matrix(
    o_mat,
    min.rows = 1, min.cols = 1, any.missing = FALSE
  )
  dims <- as.integer(dim(env_dat))

  # Hot path: full log-likelihood scalar -> pure C++.
  if (isTRUE(sum_log_p) && !isTRUE(return_prob)) {
    return(loglik_bio_cpp(
      env_dat_vec  = as.numeric(env_dat),
      env_dat_dims = dims,
      occ          = as.integer(occ),
      mu           = as.numeric(mu),
      sigltil      = as.numeric(sigltil),
      sigrtil      = as.numeric(sigrtil),
      o_mat        = as.matrix(o_mat),
      ctil         = as.numeric(ctil),
      pd           = as.numeric(pd),
      num_threads  = as.integer(num_threads)
    ))
  }

  # Slow path (per-location output / linear-scale): use C++ log_prob_detect
  # for the inner xtensor kernel and reduce in R.
  log_p <- log_prob_detect(
    env_dat = env_dat,
    mu = mu,
    sigltil = sigltil,
    sigrtil = sigrtil,
    o_mat = o_mat,
    ctil = ctil,
    pd = pd,
    return_prob = FALSE,
    num_threads = num_threads
  )

  if (sum_log_p) {
    res <- sum(occ * log_p + (1 - occ) * log1mexp(-log_p))
  } else {
    res <- occ * log_p + (1 - occ) * log1mexp(-log_p)
  }

  if (return_prob) {
    res <- exp(res)
  }

  res
}
