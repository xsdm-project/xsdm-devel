#' Distance in parameter space between two sets of parameters
#'
#' Computes a distance in parameter space between two parameter sets of the
#' xsdm model. The \code{o_mat}, \code{sigltil}, and \code{sigrtil} parameters
#' are only determined up to an equivalence class; this function returns the
#' minimum distance over all equivalence-class representatives of \code{p1},
#' using the Hungarian (Kuhn--Munkres) linear sum assignment algorithm to
#' avoid enumerating every permutation and sign flip. Distance is measured in
#' sum-squared-errors on the biological scale, except that \code{sigltil} and
#' \code{sigrtil} are inverted before comparison (that is the scale on which
#' distance is most meaningful for those parameters).
#'
#' The numerical backbone of this function is a pure-C++ implementation in
#' \code{src/dist_between_params.cpp} that (a) builds the cost matrix over
#' column pairings and (b) solves the assignment problem. A pure-R reference
#' implementation is kept as \code{xsdm:::dist_between_params_r} and is used
#' by the parity tests in \code{tests/testthat/test-dist_between_params_r_vs_cpp.R}.
#'
#' @param p1 First set of parameters. May be math-scale (a named numeric
#'   vector whose names complement \code{mask}) or biological-scale (a named
#'   list with entries \code{mu}, \code{sigltil}, \code{sigrtil}, \code{ctil},
#'   \code{pd}, \code{o_mat}).
#' @param p2 Second set of parameters; same format options as \code{p1}.
#' @param mask Same format as the \code{mask} argument to \code{\link{loglik_math}}
#'   and \code{\link{start_parms}}. Ignored if both \code{p1} and \code{p2} are
#'   on the biological scale. Otherwise the names of \code{mask} must exactly
#'   complement the names of whichever of \code{p1} / \code{p2} is on the
#'   math scale.
#' @param give_closest_rep If \code{TRUE}, also returns the member of the
#'   equivalence class of \code{p1} that attains the minimum distance
#'   (biological scale). Default \code{FALSE}.
#'
#' @return If \code{give_closest_rep} is \code{FALSE}, a single number: the
#'   distance. Otherwise a list with entries \code{distance} and
#'   \code{representative}.
#'
#' @references
#' H. W. Kuhn (1955). The Hungarian Method for the Assignment Problem.
#' \emph{Naval Research Logistics Quarterly} 2(1-2), 83--97.
#'
#' J. Munkres (1957). Algorithms for the Assignment and Transportation
#' Problems. \emph{Journal of the SIAM} 5(1), 32--38.
#'
#' R. Jonker and A. Volgenant (1987). A Shortest Augmenting Path Algorithm for
#' Dense and Sparse Linear Assignment Problems. \emph{Computing} 38, 325--340.
#'
#' K. Hornik (2005). A CLUE for CLUster Ensembles. \emph{Journal of
#' Statistical Software} 14(12). (See also the \pkg{clue} package for an
#' alternative R-level implementation of the same LSAP algorithm that
#' \code{xsdm:::dist_between_params_r} calls via \code{clue::solve_LSAP}.)
#'
#' @export
#'
#' @examples
#' # Using lists on the biological scale
#' par_list <- math_to_bio(examples$optim_par_vec)
#' par_list_equivalent <- math_to_bio(examples$optim_par_vec_equivalent)
#' dist_between_params(
#'   p1 = par_list,
#'   p2 = par_list_equivalent
#' )
#'
#' # Using vectors on the math scale
#' dist_between_params(
#'   p1 = examples$optim_par_vec,
#'   p2 = examples$optim_par_vec_equivalent
#' )
dist_between_params <- function(p1, p2, mask = NULL, give_closest_rep = FALSE) {
  # --- basic input checks -------------------------------------------------
  checkmate::assert_flag(give_closest_rep)
  if (!((is.null(mask)) || (is.numeric(mask) && (sum(is.na(mask)) == 0)))) {
    stop("mask must be NULL or a numeric vector with no missing values")
  }
  checkmate::assert_true(is.list(p1) || is.numeric(p1))
  checkmate::assert_true(is.list(p2) || is.numeric(p2))

  # --- math -> biological conversion (identical to legacy behaviour) ------
  if (is.numeric(p1)) {
    checkmate::assert_integerish(sqrt(9 + 8 * (length(p1) + length(mask))))
    p <- round((-5 + sqrt(9 + 8 * (length(p1) + length(mask)))) / 2)
    checkmate::assert_true(setequal(c(names(mask), names(p1)),
                                    names(make_mask_names(p))))
    checkmate::assert_numeric(p1, any.missing = FALSE, finite = TRUE)
    allnames <- names(make_mask_names(p))
    checkmate::assert_true(all(names(mask[is.infinite(mask)]) %in%
      c(allnames[grepl("^sig", allnames)], "pd")))
    p1 <- math_to_bio(create_param_vector_masked(p1, mask, p))
  }
  if (is.numeric(p2)) {
    checkmate::assert_integerish(sqrt(9 + 8 * (length(p2) + length(mask))))
    p <- round((-5 + sqrt(9 + 8 * (length(p2) + length(mask)))) / 2)
    checkmate::assert_true(setequal(c(names(mask), names(p2)),
                                    names(make_mask_names(p))))
    checkmate::assert_numeric(p2, any.missing = FALSE, finite = TRUE)
    allnames <- names(make_mask_names(p))
    checkmate::assert_true(all(names(mask[is.infinite(mask)]) %in%
      c(allnames[grepl("^sig", allnames)], "pd")))
    p2 <- math_to_bio(create_param_vector_masked(p2, mask, p))
  }

  # --- biological-scale checks -------------------------------------------
  for (obj in list(p1, p2)) {
    checkmate::assert_names(names(obj),
      must.include = c("mu", "sigltil", "sigrtil", "ctil", "pd", "o_mat"))
    checkmate::assert_numeric(obj$mu, finite = TRUE, any.missing = FALSE)
    checkmate::assert_numeric(obj$sigltil, any.missing = FALSE)
    checkmate::assert_numeric(obj$sigrtil, any.missing = FALSE)
    checkmate::assert_numeric(obj$o_mat, finite = TRUE, any.missing = FALSE)
    checkmate::assert_numeric(obj$ctil, finite = TRUE, any.missing = FALSE, len = 1)
    checkmate::assert_numeric(obj$pd, finite = TRUE, any.missing = FALSE, len = 1)
  }
  p <- length(p1$mu)
  checkmate::assert_true(length(p1$sigltil) == p)
  checkmate::assert_true(length(p1$sigrtil) == p)
  checkmate::assert_true(all(dim(p1$o_mat) == c(p, p)))
  checkmate::assert_true(length(p2$mu) == p)
  checkmate::assert_true(length(p2$sigltil) == p)
  checkmate::assert_true(length(p2$sigrtil) == p)
  checkmate::assert_true(all(dim(p2$o_mat) == c(p, p)))

  # --- C++ backend: cost matrix + LSAP + distance ------------------------
  res <- .dist_between_params_cpp(
    mu1      = as.numeric(p1$mu),
    sigltil1 = as.numeric(p1$sigltil),
    sigrtil1 = as.numeric(p1$sigrtil),
    o_mat1   = as.matrix(p1$o_mat),
    ctil1    = as.numeric(p1$ctil),
    pd1      = as.numeric(p1$pd),
    mu2      = as.numeric(p2$mu),
    sigltil2 = as.numeric(p2$sigltil),
    sigrtil2 = as.numeric(p2$sigrtil),
    o_mat2   = as.matrix(p2$o_mat),
    ctil2    = as.numeric(p2$ctil),
    pd2      = as.numeric(p2$pd)
  )

  if (!give_closest_rep) {
    return(res$distance)
  }

  # --- reconstruct the closest equivalence-class representative ----------
  perm    <- as.integer(res$perm)
  posneg  <- res$posneg
  pairing <- cbind(seq_len(nrow(posneg)), perm)
  posnegs <- posneg[pairing]
  perm_inv <- order(perm)
  posnegs <- posnegs[perm_inv]
  flip <- (posnegs == -1)

  rep_ec <- convert_equivalence_class(p1, flip = flip, perm = perm)
  rep_out <- list(
    mu      = p1$mu,
    sigltil = rep_ec$sigltil,
    sigrtil = rep_ec$sigrtil,
    ctil    = p1$ctil,
    pd      = p1$pd,
    o_mat   = rep_ec$o_mat
  )
  list(distance = res$distance, representative = rep_out)
}
