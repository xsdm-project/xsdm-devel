#' Parameter list of optim conditions for the Example 1
#'
#' A small list of biological-scale parameters used in examples and tests.
#' Includes the optimal environmental values (`mu`), asymmetric width parameters
#' (`sigltil`, `sigrtil`), detection-link parameters (`ctil`, `pd`), and a
#' 2x2 orthogonal matrix (`o_mat`) defining the rotation in environmental space.
#'
#' @format A named list with the following elements:
#' \describe{
#'   \item{\code{mu}}{Numeric vector of length 2: \code{c(14.110506,  6.561854)}}
#'   \item{\code{sigltil}}{Numeric vector of length 2: \code{c(0.4627375, 1.0819189)}}
#'   \item{\code{sigrtil}}{Numeric vector of length 2: \code{c(0.1051458, 1.3447374)}}
#'   \item{\code{ctil}}{Numeric scalar: \code{-18.14259}}
#'   \item{\code{pd}}{Numeric scalar in (0,1): \code{0.8940732}}
#'   \item{\code{o_mat}}{2x2 orthogonal matrix:
#'     \deqn{
#'       \\begin{pmatrix}
#'       -0.1801624 & -0.9836369 \\\\
#'       0.9836369 &  -0.1801624
#'       \\end{pmatrix}
#'     }
#'   }
#' }
#'
#' @details
#' The matrix \code{o_mat} corresponds to a rotation by approximately 0.2 radians
#' (i.e., \code{cos(0.2)} and \code{sin(0.2)} appear in the entries).
#'
#' @examples
#' example_1_param_list_example
#' str(example_1_param_list_example)
"example_1_optim_param_list"
