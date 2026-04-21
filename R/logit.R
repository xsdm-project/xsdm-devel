#' Logit function. Inverse of the expit function.
#'
#' @param x A numeric value to scale from 0 to 1 to -inf to inf
#'
#' @returns A numeric
#' @keywords internal
#'
#' @examples
#' xsdm:::logit(0.5)
logit <- function(x) {
  no <- (x < 0) | (x > 1)
  out <- numeric(length(x))
  out[no] <- NaN
  out[!no] <- log(x[!no] / (1 - x[!no]))
  dim(out) <- dim(x)
  out
}
