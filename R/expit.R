#' Functions to take the expit of numerical vectors.
#' expit exp(x)/(1 + exp(x))
#'
#' @param x A numeric value
#'
#' @returns A real vector corresponding to the expits of x
#' @export
#'
#' @examples
#' expit(0)
#' expit(0.5)
#' expit(-1)

expit <- function(x) {
  # Allow NA values in the input
  checkmate::assert_numeric(x, any.missing = TRUE)
  
  # Preallocate output with NAs
  out <- rep(NA_real_, length(x))
  
  # Identify non-missing elements
  valid <- !is.na(x)
  if (any(valid)) {
    x_valid <- x[valid]
    idx_pos <- x_valid >= 0
    idx_neg <- x_valid < 0
    
    # Compute on valid entries only
    out[valid][idx_pos] <- 1 / (1 + exp(-x_valid[idx_pos]))
    out[valid][idx_neg] <- exp(x_valid[idx_neg]) / (1 + exp(x_valid[idx_neg]))
  }
  
  out
}