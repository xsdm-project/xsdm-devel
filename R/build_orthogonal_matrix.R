#' Build an orthogonal matrix from a real-parameter vector
#'
#' Constructs a \eqn{k \times k} orthogonal matrix \eqn{O} by exponentiating a
#' skew-symmetric matrix \eqn{S} built by assigning its lower-triangular entries
#' from the input vector. Specifically, the function sets \eqn{S_{ij}}
#' (for \eqn{i>j}) from `entries`, mirrors it to enforce \eqn{S = L - L^\top},
#' and returns `expm::expm(S)`, which is guaranteed orthogonal because
#' \eqn{\exp(S)} is orthogonal whenever \eqn{S} is real skew-symmetric.
#'
#' The dimension `k` is inferred from `length(entries)` via the relation
#' \eqn{n = k(k-1)/2}, so `length(entries)` must equal a triangular number.
#'
#' @param entries A numeric vector (possibly `NULL`). If `NULL`, returns the
#'   `1 x 1` identity. Otherwise, its length must be \eqn{n = k(k-1)/2} for some
#'   integer \eqn{k \ge 2}, supplying the strictly lower-triangular entries of a
#'   skew-symmetric generator.
#'
#' @returns A `k x k` orthogonal matrix. If `entries` is `NULL`, returns
#'   `matrix(1, 1, 1)` (identity).
#'
#' @details
#' - Dimension inference uses \eqn{k = \frac{1 + \sqrt{1 + 8n}}{2}} where
#'   \eqn{n = \text{length(entries)}}. If \eqn{k} is not an integer, the input
#'   is invalid and an error is thrown.
#' - Orthogonality follows from the fact that \eqn{S^\top = -S} implies
#'   \eqn{\exp(S)^\top \exp(S) = I}.
#' - Note the function actually returns a special orthogonal matrix, i.e., the
#'   determinant is +1.
#'
#' @examples
#' # 1x1 identity (NULL input)
#' build_orthogonal_matrix(NULL)
#'
#' # 2x2 orthogonal matrix from one parameter
#' O2 <- build_orthogonal_matrix(0.0)
#' all.equal(t(O2) %*% O2, diag(2)) # should be TRUE
#'
#' # 3x3 example: length(entries) = 3 (= 3*2/2), so k = 3
#' O3 <- build_orthogonal_matrix(c(0.1, -0.2, 0.3))
#' all.equal(t(O3) %*% O3, diag(3), tolerance = 1e-10)
#'
#' @export
build_orthogonal_matrix <- function(entries) {
  # Validate type: numeric vector or NULL; no NAs allowed.
  checkmate::assert_numeric(entries, any.missing = FALSE, null.ok = TRUE)

  # NULL => 1x1 identity
  if (is.null(entries)) {
    return(matrix(1, 1, 1))
  }

  # Infer k from n = length(entries) using n = k*(k-1)/2
  n <- length(entries)
  k_real <- 0.5 * (1 + sqrt(1 + 8 * n))
  k <- as.integer(round(k_real))

  # Ensure k is a valid integer solution
  if (k <= 1L || k * (k - 1L) / 2L != n) {
    stop(
      "Invalid 'entries' length: must be a triangular number n = k*(k-1)/2 ",
      "for some integer k >= 2. Got n = ", n, "."
    )
  }

  # Build a k x k skew-symmetric matrix s_matrix from lower-triangular entries.
  s_matrix <- matrix(0, nrow = k, ncol = k)
  s_matrix[lower.tri(s_matrix)] <- entries

  # enforce skew-symmetry: s_matrix^T = -s_matrix
  s_matrix <- (s_matrix - t(s_matrix))

  # Exponentiate to obtain an orthogonal matrix: O = expm(s_matrix)
  o_matrix <- expm::expm(s_matrix)
  o_matrix
}
