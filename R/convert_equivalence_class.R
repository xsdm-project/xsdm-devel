#' Converts a set of parameters to other representatives of the same equivalence
#' class
#'
#' Model parameters in the biological scale are only determined up to an
#' equivalence class. This function converts a set of parameters to another set
#' of equivalent parameters.
#'
#' @param p Named list with entries mu, sigltil, sigrtil, ctil, pd, and o_mat
#' @param flip Vector of binaries corresponding to which columns of o_mat are to
#' have their sign change (which a concomitant switch of the corresponding
#' entries of sigltil and sigrtil). Length must equal the number of columns of
#' o_mat.
#' @param perm Permutation to be applied to the columns of o_mat.
#'
#' @return List with entries o_mat, sigltil, and sigrtil
#'
#' @export
#' @examples
#' convert_equivalence_class(
#'   p = examples$optim_par_list,
#'   flip = c(1, 0),
#'   perm = c(1, 2)
#' )
convert_equivalence_class <- function(p, flip, perm) {
  o_mat <- p$o_mat
  sigltil <- p$sigltil
  sigrtil <- p$sigrtil
  dd <- dim(o_mat)[1]

  # do the sign flipping
  for (cc in 1:dd) {
    if (flip[cc]) {
      o_mat[, cc] <- -o_mat[, cc]
      h <- sigrtil[cc]
      sigrtil[cc] <- sigltil[cc]
      sigltil[cc] <- h
    }
  }

  # do the permuting
  o_mat <- o_mat[, perm]
  sigltil <- sigltil[perm]
  sigrtil <- sigrtil[perm]
  output <- list(o_mat = o_mat, sigltil = sigltil, sigrtil = sigrtil)
  output
}
