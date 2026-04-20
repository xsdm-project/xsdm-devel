#' Extract a math-scale real-parameter vector corresponding to a special orthogonal matrix
#'
#' Computes the principal matrix logarithm of a special orthogonal matrix `o_mat`,
#' then returns the strictly lower-triangular entries of the resulting skew-symmetric
#' matrix. This is a partial inverse of `build_orthogonal_matrix`, up to the periodicity
#' of the exponential map.
#'
#' @details
#' The matrix exponential \eqn{\exp: \mathfrak{so}(k) \to SO(k)} is surjective but not
#' injective: different skew-symmetric matrices can exponentiate to the same orthogonal
#' matrix. This function uses the **principal matrix logarithm** as implemented in
#' `expm::logm`. Consequently, it may fail (or produce complex results) for matrices
#' that have eigenvalues equal to \eqn{-1} (i.e., rotations by \eqn{\pi}). Such matrices
#' lie on the cut locus of the exponential map and do not possess a unique real logarithm.
#' If you encounter this, consider perturbing the matrix slightly away from the problematic
#' rotation.
#'
#' @param o_mat A \eqn{k \times k} special orthogonal matrix (\code{t(o_mat) \%*\% o_mat = I}
#'   and \code{det(o_mat) = 1}).
#'
#' @returns A numeric vector of length \eqn{k(k-1)/2} containing the strictly
#'   lower-triangular entries of the skew-symmetric generator. For \eqn{k=1},
#'   returns \code{NULL} (the identity matrix).
#'
#' @export
#' @examples
#' o_par2 <- 0.25
#' O2 <- build_orthogonal_matrix(o_par2)
#' extract_orthogonal_matrix_parameters(O2)
extract_orthogonal_matrix_parameters <- function(o_mat) {
  tol <- 1e-10

  # Basic validation
  checkmate::assert_matrix(o_mat, mode = "numeric", any.missing = FALSE)
  k <- nrow(o_mat)
  checkmate::assert_true(nrow(o_mat) == ncol(o_mat), .var.name = "o_mat must be square")

  # 1x1 special case
  if (k == 1L) {
    if (abs(o_mat[1, 1] - 1) < tol) {
      return(NULL)
    } else if (abs(o_mat[1, 1] + 1) < tol) {
      warning("Input is a reflection (det = -1). No skew-symmetric generator
              exists. Returning NULL.")
      return(NULL)
    } else {
      stop("1x1 matrix is not orthogonal: value must be +1 or -1")
    }
  }

  # Orthogonality and determinant checks
  I <- diag(k)
  err_orth <- max(abs(o_mat %*% t(o_mat) - I))
  checkmate::assert_true(err_orth < 1e-5, .var.name = "o_mat is not orthogonal")

  det_o <- as.numeric(det(o_mat))
  checkmate::assert_true(abs(det_o - 1) < 1e-5, .var.name = "det(o_mat) != 1")

  # Check for eigenvalues near -1 (cut locus) and provide informative error
  eig <- eigen(o_mat, symmetric = FALSE, only.values = TRUE)$values
  if (any(Re(eig) < -1 + tol & abs(Im(eig)) < tol)) {
    stop(
      "o_mat has eigenvalues equal to -1 (rotation by pi). ",
      "The principal matrix logarithm is not defined. ",
      "Consider perturbing the matrix slightly."
    )
  }

  # Compute principal logarithm
  S <- expm::logm(o_mat)

  # Enforce real and skew-symmetric (within tolerance)
  if (max(abs(Im(S))) > tol) {
    warning("Logarithm has non-negligible imaginary part; real part used.")
  }
  S <- Re(S)
  if (max(abs(S + t(S))) > tol) {
    warning("Logarithm is not perfectly skew-symmetric;
            using lower triangle anyway.")
  }

  S[lower.tri(S)]
}
