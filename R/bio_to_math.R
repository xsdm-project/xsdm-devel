#' Converts parameters from the biological scale to the math (unconstrained)
#' scale
#'
#' Given a list with biological-scale parameters (`mu`, `sigltil`, `sigrtil`,
#' `ctil`, `pd`, `o_mat`), returns a named numeric vector on the math scale,
#' in the canonical order produced by `make_mask_names(p)`:
#' `mu1..mup`, `sigltil1..p`, `sigrtil1..p`, `o_par1..q`, `ctil`, `pd`,
#' where `q = p*(p-1)/2` and `p = length(mu) = nrow(o_mat)`.
#'
#' Transformations:
#' - `mu`      : identity
#' - `sigltil` : `log()`
#' - `sigrtil` : `log()`
#' - `ctil`    : identity
#' - `pd`      : `logit()`
#' - `o_mat`   : lower-triangular parameters recovered via
#'               `extract_orthogonal_matrix_parameters()` (see Details)
#'
#' @param parms_bio A named list with entries: `mu`, `sigltil`, `sigrtil`,
#'   `ctil`, `pd`, `o_mat`.
#'
#' @returns A named numeric vector on the math scale, ordered per
#'   `make_mask_names(p)`.
#'
#' @details
#' The `o_mat` entries are mapped to a vector via the principal matrix
#' logarithm, i.e., one of the (skew-symmetric) matrices S such that
#' `o_mat = expm(S)`. The math-scale parameters `o_par` are then the
#' strictly lower-triangular elements of `S`. For `p = 1`, there are no
#' `o_par` entries. Note, however, that the principal logarithm is not
#' defined for all special orthogonal matrices (even though all such
#' matrices are in the image of the matrix exponential), so this function
#' may fail for some valid `o_mat` inputs.
#'
#' @seealso [math_to_bio()], [make_mask_names()], [build_orthogonal_matrix()]
#' @export
#' @examples
#' ## --- p = 1 (no o_par entries) ---
#' mu1 <- 10
#' sigltil1 <- 1.2
#' sigrtil1 <- 0.8
#' bio_parameters <- list(
#'   mu      = c(mu1),
#'   sigltil = c(sigltil1),
#'   sigrtil = c(sigrtil1),
#'   ctil    = 0.3,
#'   pd      = 0.85,
#'   o_mat   = matrix(1, 1, 1) # 1x1 orthogonal
#' )
#' math1 <- bio_to_math(bio_parameters)
#' # Canonical names
#' names(math1)
#' # Back to biological scale
#' math_parameters <- math_to_bio(math1)
#' all.equal(math_parameters$mu, bio_parameters$mu)
#' all.equal(math_parameters$sigltil, bio_parameters$sigltil)
#' all.equal(math_parameters$sigrtil, bio_parameters$sigrtil)
#' all.equal(math_parameters$ctil, bio_parameters$ctil)
#' all.equal(math_parameters$pd, bio_parameters$pd)
#'
#' ## --- p = 2 (includes one o_par) ---
#' mu2 <- c(11, 5)
#' sigltil2 <- c(1.1, 1.5)
#' sigrtil2 <- c(1.4, 1.3)
#' ctil2 <- -0.2
#' pd2 <- 0.9
#' o_par2 <- 0.25
#' O2 <- build_orthogonal_matrix(o_par2)
#' bio_parameters_2d <- list(
#'   mu      = mu2,
#'   sigltil = sigltil2,
#'   sigrtil = sigrtil2,
#'   ctil    = ctil2,
#'   pd      = pd2,
#'   o_mat   = O2
#' )
#' math_parameters_2d <- bio_to_math(bio_parameters_2d)
#' # check canonical name order produced by make_mask_names(2)
#' identical(names(math_parameters_2d), names(make_mask_names(2)))
bio_to_math <- function(parms_bio) {
  # --- check name assertions   ---
  checkmate::assert_list(parms_bio, names = "unique", any.missing = FALSE)
  req <- c("mu", "sigltil", "sigrtil", "ctil", "pd", "o_mat")
  if (!all(req %in% names(parms_bio))) {
    stop("parms_bio must contain: ", paste(req, collapse = ", "), ".")
  }

  mu <- parms_bio$mu
  sigltil <- parms_bio$sigltil
  sigrtil <- parms_bio$sigrtil
  ctil <- parms_bio$ctil
  pd <- parms_bio$pd
  o_mat <- parms_bio$o_mat

  checkmate::assert_numeric(mu, any.missing = FALSE, min.len = 1)
  checkmate::assert_numeric(sigltil,
    any.missing = FALSE,
    len = length(mu),
    lower = 0
  )
  checkmate::assert_numeric(sigrtil,
    any.missing = FALSE,
    len = length(mu),
    lower = 0
  )
  checkmate::assert_numeric(ctil, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(pd, any.missing = FALSE, len = 1)
  checkmate::assert_matrix(o_mat,
    mode = "numeric", any.missing = FALSE,
    nrows = length(mu), ncols = length(mu)
  )

  p <- length(mu)
  q <- p * (p - 1L) / 2L

  # --- Transform to math scale ---
  mu_math <- as.numeric(mu)
  sigltil_math <- log(as.numeric(sigltil))
  sigrtil_math <- log(as.numeric(sigrtil))
  ctil_math <- as.numeric(ctil)
  pd_math <- logit(as.numeric(pd))

  # o_par from orthogonal matrix
  o_par_math <- if (p == 1L) {
    numeric(0)
  } else {
    extract_orthogonal_matrix_parameters(o_mat)
  }

  # --- Build vector in canonical order ---
  out <- make_mask_names(p) # NA vector with ordered names
  out[grep("^mu\\d+$", names(out))] <- mu_math
  out[grep("^sigltil\\d+$", names(out))] <- sigltil_math
  out[grep("^sigrtil\\d+$", names(out))] <- sigrtil_math
  out["ctil"] <- ctil_math
  out["pd"] <- pd_math
  # o_par
  if (q > 0L) {
    o_inds <- grep("^o_par\\d+$", names(out))
    if (length(o_par_math) != length(o_inds)) {
      stop(
        "Length mismatch: extracted o_par has length ", length(o_par_math),
        " but expected ", length(o_inds), "."
      )
    }
    out[o_inds] <- o_par_math
  }


  if (anyNA(out)) {
    miss <- names(out)[is.na(out)]
    stop(sprintf(
      "bio_to_math produced NA values for: %s.",
      paste(miss, collapse = ", ")
    ))
  }

  out
}
