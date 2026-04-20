#' Generalized expit function
#'
#' @param x A real number to ma
#' @param l_bound Lower bound
#' @param u_bound Upper bound
#' @param x0 Middle point
#'
#' @returns A re-scaled valued between lower bound and upper bound with
#' the middle point x0
#' @export
#'
#' @examples
#' expit_general(3, -pi, pi, 0)
expit_general <- function(x, l_bound, u_bound, x0) {
  checkmate::assert_numeric(x)
  checkmate::assert_numeric(l_bound)
  checkmate::assert_numeric(u_bound)
  out <- l_bound + (u_bound - l_bound) / (1 + exp(-(x - x0)))
  out[x > 100] <- u_bound
  out[x < -100] <- l_bound
  out
}
