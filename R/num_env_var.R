#' Get the number of environmental variables given the number of parameters
#'
#' Inverts \code{num_par(p)} to recover \code{p} from \code{n}, the number of
#' parameters of the main xsdm model. Uses the closed-form solution of the
#' quadratic: \eqn{2n = p^2 + 5p + 4}, i.e. \eqn{p = (-5 + \sqrt{9 + 8n})/2}.
#' Errors if \eqn{n} is not a valid value of \code{num_par(p)} for some integer
#' \eqn{p \ge 1}.
#'
#' @param n Integerish scalar: total number of parameters.
#'
#' @returns A single integer \code{p}, the number of environmental variables.
#'
#' @export
#' @examples
#' num_env_var(5) # -> 1  (since num_par(1) = 5)
#' num_env_var(9) # -> 2  (since num_par(2) = 9)
#' num_env_var(14) # -> 3  (since num_par(3) = 14)
#' # round-trip check:
#' p <- 4
#' stopifnot(num_env_var(num_par(p)) == p)
num_env_var <- function(n) {
  checkmate::assert_integerish(n, lower = 5, any.missing = FALSE, len = 1)

  disc <- 9 + 8 * as.numeric(n)
  sqrt_disc <- sqrt(disc)
  p <- (-5 + sqrt_disc) / 2

  # require integer p >= 1 and perfect-square discriminant
  if (!isTRUE(all.equal(sqrt_disc, round(sqrt_disc))) || p < 1 || !isTRUE(all.equal(p, round(p)))) {
    stop(
      "Invalid 'n': must equal num_par(p) for some integer p >= 1. 2n = p^2 + 5p + 4 must hold for a positive integer p."
    )
  }
  as.integer(round(p))
}
