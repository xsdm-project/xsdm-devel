#' Example MLE-like parameter vector on the math scale (p = 2)
#'
#' A small named numeric vector of math-scale parameters used in examples and
#' tests for the New Mexico virtual species (p = 2). The names follow the
#' canonical conventions described in \code{loglik_math()} details.
#'
#' @format A named numeric vector of length 9 with the following entries:
#' \describe{
#'   \item{\code{mu1}, \code{mu2}}{Optimal environmental values (math scale).}
#'   \item{\code{sigltil1}, \code{sigltil2}}{Left-width parameters (math scale; exponentiated in \code{math_to_bio()}).}
#'   \item{\code{sigrtil1}, \code{sigrtil2}}{Right-width parameters (math scale; exponentiated in \code{math_to_bio()}).}
#'   \item{\code{ctil}}{Detection-link center parameter (math scale).}
#'   \item{\code{pd}}{Detection maximum parameter (math scale; expit-transformed in \code{math_to_bio()}).}
#'   \item{\code{o_par1}}{Orthogonal-matrix parameter (math scale) for p = 2 (lower-triangle of the skew-symmetric matrix).}
#' }
#'
#' @details
#' This object is intended for test the distance_between_params_Hungarian function
#' comparing example_1_optim_param_vector with this example_1_optim_param_vector_equivalent
"example_1_optim_param_vector_equivalent"
