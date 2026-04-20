#' Enumerate Permutations of Vector Elements
#'
#' `permutations()` enumerates all permutations of size `r` drawn from a source
#' vector of length `n`. The result is a matrix with one permutation per row.
#'
#' **Caution:** The number of permutations grows rapidly with `n` and `r`.
#' For large `n`, you may need to increase R's expression/recursion limit;
#' see `options(expressions = ...)`.
#'
#' @aliases permutations combinations
#' @param n Integer. Size of the source vector.
#' @param r Integer. Size of the target vectors (permutation length).
#' @param v Vector. Source vector; defaults to `1:n`.
#' @param set Logical. If `TRUE`, duplicates in `v` are removed. Default `TRUE`.
#' @param repeats.allowed Logical. If `TRUE`, permutations may include repeated
#'   values. Default `FALSE`.
#'
#' @return A matrix with `choose(n, r) * r!` rows (or fewer with `set = TRUE`
#'   and duplicated `v`), and `r` columns, where each row is one permutation.
#'
#' @details
#' This interface is adapted from the implementation in the \pkg{gtools} package
#' (function `gtools::permutations`). The original function was derived from code
#' by Bill Venables, extended by Gregory R. Warnes to handle `repeats.allowed`.
#'
#' @seealso [base::choose()], [base::options()]
#'
#' @references
#' Venables, W. N. (2001). Programmer's Niche. *R News*, 1(1).
#'   <https://cran.r-project.org/doc/Rnews/>
#'
#' @note
#' If you copied or adapted code from \pkg{gtools}, ensure your package license
#' is GPL-2 compatible and that attribution is provided (see package `LICENSE`
#' and `Authors@R`).
#'
permutations <- function(n, r, v = 1:n, set = TRUE, repeats.allowed = FALSE) {
  if (mode(n) != "numeric" || length(n) != 1 || n < 1 || (n %% 1) !=
    0) {
    stop("bad value of n")
  }
  if (mode(r) != "numeric" || length(r) != 1 || r < 1 || (r %% 1) !=
    0) {
    stop("bad value of r")
  }
  if (!is.atomic(v) || length(v) < n) {
    stop("v is either non-atomic or too short")
  }
  if ((r > n) & repeats.allowed == FALSE) {
    stop("r > n and repeats.allowed=FALSE")
  }
  if (set) {
    v <- unique(sort(v))
    if (length(v) < n) {
      stop("too few different elements")
    }
  }
  v0 <- vector(mode(v), 0)
  if (repeats.allowed) {
    sub <- function(n, r, v) {
      if (r == 1) {
        matrix(v, n, 1)
      } else if (n == 1) {
        matrix(v, 1, r)
      } else {
        inner <- Recall(n, r - 1, v)
        cbind(rep(v, rep(nrow(inner), n)), matrix(t(inner),
          ncol = ncol(inner), nrow = nrow(inner) * n,
          byrow = TRUE
        ))
      }
    }
  } else {
    sub <- function(n, r, v) {
      if (r == 1) {
        matrix(v, n, 1)
      } else if (n == 1) {
        matrix(v, 1, r)
      } else {
        X <- NULL
        for (i in 1:n) {
          X <- rbind(X, cbind(v[i], Recall(n -
            1, r - 1, v[-i])))
        }
        X
      }
    }
  }
  sub(n, r, v[1:n])
}
