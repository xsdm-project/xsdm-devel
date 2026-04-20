#' Probability of detection of the species in each location
#'
#' Computes the probability of detection of the species in each location for the
#' xsdm model, given environmental data and model parameters.
#'
#' @param env_dat The environmental data array, dimensions
#' \code{n_loc x n_time x p} (number of locations x time-series length x number
#' of environmental variables). Must be a 3-dimensional array with no missing
#' values.
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
#' @param return_prob Logical (default FALSE). Flag to return probabilities of
#' detection instead their logs.
#' @param num_threads Number of threads for parallel computation. Defaults to
#' \code{RcppParallel::defaultNumThreads()}.
#'
#' @returns A vector of length equal to the number of locations, containing the
#' probabilities of detection (or their logs) of the species in each location.
#' @export
#'
#' @examples
#' mu <- c(-1, 5.046939)
#' sigltil <- c(1.036834, 1.556083)
#' sigrtil <- c(1.538972, 1.458738)
#' ctil <- -2
#' pd <- 0.9
#' o_mat <- matrix(c(-0.4443546, 0.8958510, -0.8958510, -0.4443546), ncol = 2)
log_prob_detect <- function(env_dat,
                            mu,
                            sigltil,
                            sigrtil,
                            o_mat,
                            ctil,
                            pd,
                            return_prob = FALSE,
                            num_threads = RcppParallel::defaultNumThreads()) {
  check_env_array(env_dat)
  old <- RcppParallel::defaultNumThreads()
  on.exit(RcppParallel::setThreadOptions(numThreads = old), add = TRUE)
  RcppParallel::setThreadOptions(numThreads = num_threads)

  h <- like_neg_ltsgr_cpp(
    env_dat = env_dat,
    mu = mu,
    sigltil = sigltil,
    sigrtil = sigrtil,
    o_mat = o_mat
  )


  # Get probability of detection
  logpdetect <- log(pd) - log1pexp(ctil + h)
  
  if (return_prob) {
    logpdetect <- exp(logpdetect)
  }
  return(logpdetect)
}
