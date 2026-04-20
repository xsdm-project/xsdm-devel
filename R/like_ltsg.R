#' Compute likelihood for LTSG model
#'
#' This function calculates a likelihood-like measure using orthogonal matrices,
#' environmental data, and diagonal matrices, leveraging parallel computation.
#'
#' @name like_ltsg
#' @title Compute likelihood for LTSG model
#' @param mu Numeric vector of means (length equal to number of rows in `env_m`)
#' @param env_m Numeric matrix of environmental data. Must be column-major with
#'   time varying fastest: column \code{j} corresponds to \code{(location, time)}
#'   via index \code{j*q + i}. This matches the memory layout expected by the
#'   underlying C++ implementation.
#' @param dl_mat Diagonal matrix (as NumericMatrix)
#' @param drl_mat Diagonal matrix (as NumericMatrix)
#' @param ortho_m Numeric matrix (orthogonal basis)
#' @param q Integer, number of rows for reshaping
#' @param r Integer, number of columns for reshaping.
#'
#' @return A numeric vector of length `r` with computed sums.
#' @examples
#' mu <- c(1, 2)
#' ortho_m <- matrix(1:4, nrow = 2)
#' env_m <- matrix(1:4, nrow = 2)
#' dl_mat <- diag(2)
#' drl_mat <- diag(2)
#' like_ltsg(mu, env_m, dl_mat, drl_mat, ortho_m, q = 1, r = 2)
#'
#' @export
like_ltsg
