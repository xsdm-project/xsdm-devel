#' Long-term stochastic growth rate worker function for the xsdm model
#'
#' Computes the negative of the long-term stochastic growth rate, plus
#' \code{log(lambda_max)}, for the xsdm model, for each location. This is the
#' fast version of a worker function that relies on C++ code, optimized using
#' \pkg{RcppParallel}. The legacy pure-R reference is preserved as
#' \code{xsdm:::like_neg_ltsgr_r} for testing and comparison.
#'
#' @param env_dat The environmental data array, dimensions
#' (number of locations) x (time series length) x (number of environmental
#' variables). Must not contain missing values.
#' @param mu Vector of optimal environmental values.
#' Length \code{p=dim(env_dat)[3]}. Unconstrained real numbers.
#' @param sigltil Vector specifying width of the growth-environment function.
#' Length \code{p=dim(env_dat)[3]}. Positive real numbers, \code{Inf} entries
#' also allowed.
#' @param sigrtil Vector specifying width of the growth-environment function.
#' Length \code{p=dim(env_dat)[3]}. Positive real numbers, \code{Inf} entries
#' also allowed.
#' @param o_mat An orthogonal matrix, dimensions \code{p} by \code{p}.
#' @param num_threads Number of threads for parallel computation. Defaults to
#' \code{RcppParallel::defaultNumThreads()}.
#'
#' @returns A vector of length equal to the number of locations, as described
#' above.
#' @export
#'
#' @details
#' Internally, this function:
#' \enumerate{
#'   \item Reshapes the environmental data into a matrix.
#'   \item Computes inverse (diagonal) matrices for asymmetry adjustments.
#'   \item Calls the C++ function \code{like_ltsg()} for efficient likelihood
#'   computation.
#' }
#'
#' @note
#' Ensure that \code{env_dat} has no missing values. The parameter vectors
#' \code{mu}, \code{sigltil}, and \code{sigrtil} must have length equal to the
#' number of environmental variables (\code{p}). The argument \code{o_mat} must
#' be a \code{p x p} orthogonal matrix (i.e., \code{o_mat \%*\% t(o_mat)} is
#' the identity).
#'
#' @examples
#' # Example usage:
#' like_neg_ltsgr(env_dat = example_1$env_array,
#'                mu      = example_1$par_list$mu,
#'                sigltil = example_1$par_list$sigltil,
#'                sigrtil = example_1$par_list$sigrtil,
#'                o_mat   = example_1$par_list$o_mat)
like_neg_ltsgr <- function(env_dat,
                           mu,
                           sigltil,
                           sigrtil,
                           o_mat,
                           num_threads = RcppParallel::defaultNumThreads()) {
  # ---- Assertions ----
  checkmate::assert_array(env_dat,
    min.d = 3, any.missing = FALSE,
    null.ok = FALSE
  )
  checkmate::assert_numeric(mu,
    any.missing = FALSE,
    min.len = 1
  )
  checkmate::assert_numeric(sigltil,
    lower = 0,
    any.missing = FALSE,
    len = length(mu)
  )
  checkmate::assert_numeric(sigrtil,
    lower = 0, any.missing = FALSE,
    len = length(mu)
  )
  checkmate::assert_matrix(o_mat,
    nrows = length(mu), ncols = length(mu),
    any.missing = FALSE
  )
  checkmate::assert_number(num_threads, lower = 1)

  # ---- Set threads ----
  RcppParallel::setThreadOptions(numThreads = num_threads)

  # ---- Dimensions ----
  n <- dim(env_dat)[1] # number of locations
  ts_length <- dim(env_dat)[2] # time steps
  p <- length(mu) # number of environmental variables

  # ---- Permutate the array env_dat ----
  env_dat <- aperm(env_dat, c(3, 2, 1))

  # ---- Reshape the array env_dat to a matrix ----
  # env_dat_mat is column-major: time varies fastest, then location.
  # Each column j corresponds to (location l, time t) via j = (l-1)*ts_length + t.
  env_dat_mat <- matrix(env_dat, nrow = p, ncol = ts_length * n)

  # Verify that the column-major reshape is consistent with the dimensions
  # passed to the C++ core (q = ts_length, r = n), matching the guard in like_ltsg.cpp.
  if (ts_length * n != ncol(env_dat_mat)) {
    stop(
      "Dimension mismatch: ts_length * n (", ts_length * n,
      ") must equal ncol(env_dat_mat) (", ncol(env_dat_mat),
      "). env_dat_mat is expected to be column-major with time varying fastest."
    )
  }

  # ---- Compute inverse matrices ----
  if (p == 1) {
    dl_inv <- matrix(1 / sigltil, 1, 1)
    dr_inv <- matrix(1 / sigrtil, 1, 1)
  } else {
    dl_inv <- diag(1 / sigltil)
    dr_inv <- diag(1 / sigrtil)
  }
  drl_inv <- dr_inv - dl_inv

  # ---- Call C++ core function ----
  res <- like_ltsg(
    env_m = env_dat_mat,
    mu = mu,
    dl_mat = dl_inv,
    drl_mat = drl_inv,
    ortho_m = t(o_mat),
    q = ts_length,
    r = n
  )

  # ---- Reset threads ----
  RcppParallel::setThreadOptions(numThreads = RcppParallel::defaultNumThreads())

  res
}
