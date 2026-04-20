#' Function to facilitate the creation of the argument \code{mask} to the function
#' \code{loglik_math}
#'
#' The argument \code{mask} to the function \code{loglik_math} is required to
#' follow some very specific conventions in order to reduce the risk of errors
#' coming from mismatched arguments. This function facilitates the creation of
#' such vectors.
#'
#' @param p A positive integer representing the number of environmental variables
#' to be used in the xsdm model, i.e., \code{dim(env_dat)[3]} for the
#' \code{env_dat} argument of the \code{loglik_math} function
#'
#' @returns A named numeric vector full of NAs, with the names generated according
#' to the conventions in Details of the function \code{loglik_math}. See also Details
#' below.
#' @export
#'
#' @details
#' The output has length \code{3*p+(p^2-p)/2+2}. The names of the entries are
#' \code{mu1}, \code{mu2}, \ldots, \code{mup}, \code{sigltil1}, \code{sigltil2},
#' \ldots, \code{sigltilp}, \code{sigrtil1}, \code{sigrtil2}, \ldots,
#' \code{sigrtilp}, \code{o_pari} for \code{i} ranging from 1 to \code{(p^2-p)/2},
#' \code{ctil}, and \code{pd}. All entries are \code{NA}.
#'
#' @examples
#' make_mask_names(2)
make_mask_names <- function(p) {
  checkmate::assert_count(p, positive = TRUE)

  q <- (p^2 - p) / 2

  # Build names in the documented order; omit o_par when q == 0
  o_names <- if (q > 0) paste0("o_par", seq_len(q)) else character(0)

  names_order <- c(
    paste0("mu", seq_len(p)),
    paste0("sigltil", seq_len(p)),
    paste0("sigrtil", seq_len(p)),
    "ctil",
    "pd",
    o_names
  )

  # Ensure names length matches num_par(p)
  expected_len <- num_par(p)
  if (length(names_order) != expected_len) {
    stop(
      "Internal error: names length (", length(names_order),
      ") does not match num_par(p) = ", expected_len, "."
    )
  }

  stats::setNames(rep(NA_real_, length(names_order)), names_order)
}
