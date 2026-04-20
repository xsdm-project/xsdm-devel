#' Numerically stable `log(1 - exp(-a))`
#'
#' Computes \eqn{\log(1 - \exp(-a))} accurately for non-negative `a`,
#' using two different formulas depending on whether `a` is above or
#' below `log(2)`.
#'
#' @param a Numeric vector of non-negative values. `NA` values are
#'   preserved; negative values emit a warning and return `NaN`.
#' @param cutoff Positive numeric scalar. Threshold between the two
#'   formulas; `log(2)` is near-optimal.
#'
#' @return A numeric vector the same length as `a` with
#'   \eqn{\log(1 - \exp(-a))}.
#'
#' @references Mächler, M. (2012). *Accurately Computing
#'   log(1 − exp(− |a|)).* CRAN package `copula` vignette.
#' @seealso \code{\link{log1pexp}}, \code{\link{log1p}}, \code{\link{expm1}}
#' @examples
#' a <- 2^seq(-20, 5, length.out = 10)
#' cbind(a, log(1 - exp(-a)), log1mexp(a))
#' @export
log1mexp <- function(a, cutoff = log(2)) {
  if (has_na <- any(ina <- is.na(a))) {
    y <- a
    a <- a[ok <- !ina]
  }
  if (any(a < 0)) {
    warning("'a' >= 0 needed")
  }
  tst <- a <= cutoff
  r <- a
  r[tst] <- log(-expm1(-a[tst]))
  r[!tst] <- log1p(-exp(-a[!tst]))
  if (has_na) {
    y[ok] <- r
    y
  } else {
    r
  }
}
