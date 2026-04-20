#' Long-term stochastic growth rate worker function for the xsdm model, R version
#'
#' Computes the negative of the long-term stochastic growth rate, plus
#' log(lambda_max), for the xsdm model, for each location. This is the R version
#' of a worker function, see also the accompanying C version, which should
#' produce identical results but faster.
#'
#' @param env_dat The environmental data array, dimensions
#' (number of locations) x (time series length) x (number of environmental
#' variables). Must not contain missing values.
#' @param mu Vector of optimal environmental values. Length \code{p=dim(env_dat)[3]}.
#' Unconstrained real numbers.
#' @param sigltil Vector specifying width of the growth-environment function.
#' Length \code{p=dim(env_dat)[3]}. Positive real numbers, \code{Inf} entries also allowed.
#' @param sigrtil Vector specifying width of the growth-environment function.
#' Length \code{p=dim(env_dat)[3]}. Positive real numbers, \code{Inf} entries also allowed.
#' @param o_mat An orthogonal matrix, dimensions \code{p} by \code{p}.
#'
#' @returns A vector of length equal to the number of locations, as described
#' above.
#'
#' @details Being an internal function, there is no error checking. Note that
#' env_dat must be a 3d array (not a matrix or a vector) even if one of its
#' dimensions is 1. And \code{o_mat} must be a matrix even when \code{p} is 1
#' (in that case it's a 1 x 1 matrix, but not a scalar).
#'
#' @keywords internal
like_neg_ltsgr_r <- function(env_dat,
                             mu,
                             sigltil,
                             sigrtil,
                             o_mat) {
  # get various dimensions for convenience
  n <- dim(env_dat)[1] # number of locations
  tslen <- dim(env_dat)[2] # time series length
  p <- length(mu) # number of env vars

  # subtract mu and apply o_mat to get to u
  envdat_mat <- matrix(aperm(env_dat, c(3, 2, 1)), nrow = p, ncol = tslen * n)
  u <- t(o_mat) %*% (envdat_mat - matrix(mu, p, tslen * n))

  # apply the asymmetries
  if (p == 1) {
    DLinv <- matrix(1 / sigltil, 1, 1)
    DRinv <- matrix(1 / sigrtil, 1, 1)
  } else {
    DLinv <- diag(1 / sigltil)
    DRinv <- diag(1 / sigrtil)
  }
  uasymsq <- (DLinv %*% u + (DRinv - DLinv) %*% matrix(pmax(0, u), p, tslen * n))^2

  # compute the negative of the long-term stochastic growth rate, plus
  # log(lambda_max), in each location and return
  res <- matrix(apply(FUN = sum, X = uasymsq, MARGIN = 2), tslen, n)
  res <- 0.5 * apply(FUN = mean, X = res, MARGIN = 2)

  return(res)
}
