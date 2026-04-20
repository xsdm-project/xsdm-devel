#' Numerically stable `log(1 + exp(x))`
#'
#' Computes \eqn{\log(1 + \exp(x))} accurately for any real `x`,
#' avoiding overflow as `x -> +Inf` and catastrophic cancellation as
#' `x -> -Inf`.
#'
#' @param x Numeric vector. `NA` values are preserved.
#' @param c0,c1,c2 Numeric scalars defining the switch points between
#'   four asymptotically optimal formulas. Defaults (-37, 18, 33.3)
#'   are from Mächler (2012) and should not normally be changed.
#'
#' @return A numeric vector the same length as `x` with
#'   \eqn{\log(1 + \exp(x))}.
#'
#' @references Mächler, M. (2012). *Accurately Computing
#'   log(1 − exp(− |a|)).* CRAN package `copula` vignette.
#' @seealso \code{\link{log1mexp}}, \code{\link{log1p}},  \code{\link{expm1}} 
#' @examples
#' x <- seq(-40, 40, by = 10)
#' cbind(x, log1p(exp(x)), log1pexp(x))
#' @export

log1pexp <- function(x, c0 = -37, c1 = 18, c2 = 33.3) {
  if (has_na <- any(ina <- is.na(x))) {
    y <- x
    x <- x[ok <- !ina]
  }
  r <- exp(x)
  if (any(i <- c0 < x & (i1 <- x <= c1))) {
    r[i] <- log1p(r[i])
  }
  if (any(i <- !i1 & (i2 <- x <= c2))) {
    r[i] <- x[i] + 1 / r[i]
  }
  if (any(i3 <- !i2)) {
    r[i3] <- x[i3]
  }
  if (has_na) {
    y[ok] <- r
    y
  } else {
    r
  }
}
