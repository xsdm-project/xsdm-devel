#' Example parameter list for the New Mexico virtual species (p = 2)
#'
#' A small list of biological-scale parameters used in examples and tests.
#' Includes the optimal environmental values (`mu`), asymmetric width parameters
#' (`sigltil`, `sigrtil`), detection-link parameters (`ctil`, `pd`), and a
#' 2x2 orthogonal matrix (`o_mat`) defining the rotation in environmental space.
#'
#' @format A named list with the following elements:
#' \describe{
#'   \item{\code{mu}}{Numeric vector of length 2: \code{c(14.2, 6.8)}}
#'   \item{\code{sigltil}}{Numeric vector of length 2: \code{c(1.4, 0.5)}}
#'   \item{\code{sigrtil}}{Numeric vector of length 2: \code{c(1.35, 0.10)}}
#'   \item{\code{ctil}}{Numeric scalar: \code{-17.6}}
#'   \item{\code{pd}}{Numeric scalar in (0,1): \code{0.9}}
#'   \item{\code{o_mat}}{2x2 orthogonal matrix:
#'     \deqn{
#'       \\begin{pmatrix}
#'       0.9800666 & -0.1986693 \\\\
#'       0.1986693 &  0.9800666
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
"example_1_param_list_example"
