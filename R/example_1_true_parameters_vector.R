#' Example of vector of true parameters on the math scale for the New Mexico
#' virtual species and two bioclimatic variables (p = 2).
#'
#' A small named numeric vector of **math-scale** parameters used in examples and
#' tests. This object is the math-scale representation corresponding to the
#' biological-scale list \code{example_1_param_list_example}; i.e., it is
#' essentially \code{bio_to_math(example_1_param_list_example)} expressed as a
#' single named vector following the canonical parameter naming conventions.
#'
#' The entries include:
#' \itemize{
#'   \item \code{mu1}, \code{mu2} (unconstrained on the math scale),
#'   \item \code{sigltil1}, \code{sigltil2} and \code{sigrtil1}, \code{sigrtil2}
#'         (log-scale widths; these correspond to \code{log(sigltil)} and
#'         \code{log(sigrtil)} on the biological scale),
#'   \item \code{ctil} (unconstrained),
#'   \item \code{pd} (logit / expit scale; on the biological scale \code{pd = expit(pd)}),
#'   \item \code{o_par1} (orthogonal-matrix parameter for \code{p = 2}).
#' }
#'
#' @format A named numeric vector of length 9 with the following values:
#' \describe{
#'   \item{\code{mu1}}{14.2}
#'   \item{\code{mu2}}{6.8}
#'   \item{\code{sigltil1}}{0.3364722}
#'   \item{\code{sigltil2}}{-0.6931472}
#'   \item{\code{sigrtil1}}{0.3001046}
#'   \item{\code{sigrtil2}}{-2.example_1_true_parameters_vector3025851}
#'   \item{\code{ctil}}{-17.6}
#'   \item{\code{pd}}{2.1972246}
#'   \item{\code{o_par1}}{0.2}
#' }
#'
#' @details
#' This object is useful for reproducible examples that require a canonical
#' named vector on the math scale (e.g., \code{loglik_math()}, profiling utilities,
#' and tests). For the corresponding biological-scale parameters, see
#' \code{\link{example_1_param_list_example}}.
#'
#' @examples
#' example_1_true_parameters_vector
#' str(example_1_true_parameters_vector)
#'
#' ## Relationship to the biological-scale example (if bio_to_math() is available):
#' ## all.equal(example_1_true_parameters_vector, bio_to_math(example_1_param_list_example))
"example_1_true_parameters_vector"
