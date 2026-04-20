#' Log-likelihood function for the xsdm model, parameters on the biological
#' scale.
#'
#' Computes the log-likelihood for the xsdm model given environmental data, a
#' vector of occurrences and pseudo-absences, and model parameters on the
#' biological scale.
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
#' mu <- param_list_example$mu
#' sigltil <- param_list_example$sigltil
#' sigrtil <- param_list_example$sigrtil
#' ctil <- param_list_example$ctil
#' pd <- param_list_example$pd
#' o_mat <- param_list_example$o_mat
#' env_dat <- example_1_env_array
#' occ <- example_1_occurrence_vector
#'
#' ll <- loglik_bio(
#'   env_dat,
#'   occ,
#'   mu = mu,
#'   sigltil = sigltil,
#'   sigrtil = sigrtil,
#'   o_mat = o_mat,
#'   ctil = ctil,
#'   pd = pd
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
  # Validate inputs for modeling function --------------------------------------
  # occ: must be either a logical vector (TRUE/FALSE) with no NAs or a numeric
  # or integer  vector containing only 0 and 1 with no NA

  #   Using a disjunctive assert so either condition is acceptable
  checkmate::assert(
    checkmate::check_logical(occ, any.missing = FALSE),
    checkmate::check_integerish(occ, lower = 0, upper = 1, any.missing = FALSE),
    .var.name = "occ"
  )

  # env_dat: must be a 3-dimensional array (n_loc x n_time x p) with no NAs.

  # This prevents passing a vector or 1D array by mistake.
  check_env_array(env_dat)

  # mu: numeric vector (length >= 1) with no missing values.
  checkmate::assert_numeric(mu, any.missing = FALSE, min.len = 1)

  # sigl: numeric vector (length >= 1) with no missing values. Left side scale
  # of asymmetrical long term stochastic growht function)
  checkmate::assert_numeric(sigltil, any.missing = FALSE, min.len = 1)

  # sigl: numeric vector (length >= 1) with no missing values. Right side scale
  # of asymmetrical long term stochastic growht function)
  checkmate::assert_numeric(sigrtil, any.missing = FALSE, min.len = 1)

  # ctil: single numeric scalar (len == 1) with no missing values. Threshold
  # parameter; enforcing scalar avoids vector to be automatically repeated
  # (recycled) to match the length of longer vectors in operations
  # without warning
  checkmate::assert_numeric(ctil, any.missing = FALSE, len = 1)

  # pd: single numeric scalar (len == 1) with no missing values. Penalty
  # in the probability of detection; enforcing scalar (see above ctil).
  checkmate::assert_numeric(pd, any.missing = FALSE, len = 1)

  # o_mat: numeric matrix with at least 1 row and 1 column and no NAs.
  # Observation/occurrence matrix; dimensions must be valid and no missing.
  checkmate::assert_matrix(
    o_mat,
    min.rows = 1, min.cols = 1, any.missing = FALSE
  )

  # establish the desired number of threads to use. Is set as defaultNumThreads
  RcppParallel::setThreadOptions(numThreads = num_threads)
  # Restore previous thread setting on exit (success or failure) so that the
  # calling session/worker thread state is not permanently altered (Fix 4).
  on.exit(
    RcppParallel::setThreadOptions(numThreads = RcppParallel::defaultNumThreads()),
    add = TRUE
  )
  
  # get the probability of detection for each location
  log_p <- log_prob_detect(
    env_dat = env_dat,
    mu = mu,
    sigltil = sigltil,
    sigrtil = sigrtil,
    o_mat = o_mat,
    ctil = ctil,
    pd = pd,
    return_prob = FALSE
  )
  
  # If sum_log_p is TRUE, the user wants the location-specific log-likelihoods
  # to be summed, otherwise they want them separately as a vector.
  if (sum_log_p) {
    res <- sum(occ * log_p + (1 - occ) * log1mexp(-log_p))
  } else {
    res <- occ * log_p + (1 - occ) * log1mexp(-log_p)
  }

  # If return_prob is TRUE the user wants linear-scale instead of log-scale
  # likelihoods.
  if (return_prob) {
    res <- exp(res)
  }
  
  res
}
