#' Get the number of parameters of the main xsdm model given the number
#' of environmental variables to be considered
#'
#' @param p A positive integer representing the number of environmental variables
#' to be used in the xsdm model, i.e., \code{dim(env_dat)[3]} for the
#' \code{env_dat} argument of the \code{loglik_math} function
#'
#' @returns An integer with the number of parameters
#'
#' @details
#' For instance, in the `p=1` case, the xsdm model parameters are `mu`,
#' `sigltil`, and `sigrtil` (which are scalars in the `p=1` case); `ctil`,
#' and `pd` (which are scalars for any value `p`). That makes 5 parameters,
#' so this function returns 5. In the `p=2` case, the parameters are `mu`,
#' `sigltil`, and `sigrtil` (each of which is now a length-2 vector); `ctil`,
#' and `pd` (again scalars); and the single parameter pertaining to `o_mat`;
#' for a total of 9.
#'
#' @export
#' @examples
#' num_par(2)
num_par <- function(p) {
  # ctil, plus the O params, plus the mu, sigltil and sigrtil params, plus pd
  1 + (p^2 - p) / 2 + 3 * p + 1
}
