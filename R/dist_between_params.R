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
#' All numerics besides the linear sum assignment are computed in R. The
#' assignment problem itself is solved by an in-package C++ implementation of
#' the classical O(n^3) Hungarian algorithm, exposed (unexported) as
#' \code{xsdm:::.solve_lsap_cpp}. An R-level alternative is
#' \code{clue::solve_LSAP}; the two are compared in
#' \code{tests/testthat/test-solve_lsap_cpp.R}.
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
#' Statistical Software} 14(12). (See also the \pkg{clue} package on CRAN
#' for an alternative R-level LSAP implementation.)
#'
#' @export
#'
#' @examples
#' # Using lists on the biological scale
#' par_list <- math_to_bio(example_1$optim_par_vec)
#' par_list_equivalent <- math_to_bio(example_1$optim_par_vec_equivalent)
#' dist_between_params(
#'   p1 = par_list,
#'   p2 = par_list_equivalent
#' )
#'
#' # Using vectors on the math scale
#' dist_between_params(
#'   p1 = example_1$optim_par_vec,
#'   p2 = example_1$optim_par_vec_equivalent
#' )
dist_between_params <- function(p1, p2, mask = NULL, give_closest_rep = FALSE) {
  # --- basic input checks -----------------------------------------------
  checkmate::assert_flag(give_closest_rep)
  if (!((is.null(mask)) || (is.numeric(mask) && (sum(is.na(mask)) == 0)))) {
    stop("mask must be NULL or a numeric vector with no missing values")
  }
  checkmate::assert_true(is.list(p1) || is.numeric(p1))
  checkmate::assert_true(is.list(p2) || is.numeric(p2))

  # --- math -> biological conversion ------------------------------------
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

  # --- biological-scale checks -----------------------------------------
  checkmate::assert_names(names(p1),
    must.include = c("mu", "sigltil", "sigrtil", "ctil", "pd", "o_mat"))
  checkmate::assert_numeric(p1$mu, finite = TRUE, any.missing = FALSE)
  checkmate::assert_numeric(p1$sigltil, any.missing = FALSE)
  checkmate::assert_numeric(p1$sigrtil, any.missing = FALSE)
  checkmate::assert_numeric(p1$o_mat, finite = TRUE, any.missing = FALSE)
  checkmate::assert_numeric(p1$ctil, finite = TRUE, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(p1$pd, finite = TRUE, any.missing = FALSE, len = 1)
  p <- length(p1$mu)
  checkmate::assert_true(length(p1$sigltil) == p)
  checkmate::assert_true(length(p1$sigrtil) == p)
  checkmate::assert_true(all(dim(p1$o_mat) == c(p, p)))

  checkmate::assert_names(names(p2),
    must.include = c("mu", "sigltil", "sigrtil", "ctil", "pd", "o_mat"))
  checkmate::assert_numeric(p2$mu, finite = TRUE, any.missing = FALSE)
  checkmate::assert_numeric(p2$sigltil, any.missing = FALSE)
  checkmate::assert_numeric(p2$sigrtil, any.missing = FALSE)
  checkmate::assert_numeric(p2$o_mat, finite = TRUE, any.missing = FALSE)
  checkmate::assert_numeric(p2$ctil, finite = TRUE, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(p2$pd, finite = TRUE, any.missing = FALSE, len = 1)
  checkmate::assert_true(length(p2$mu) == p)
  checkmate::assert_true(length(p2$sigltil) == p)
  checkmate::assert_true(length(p2$sigrtil) == p)
  checkmate::assert_true(all(dim(p2$o_mat) == c(p, p)))

  sigdistsq <- function(x, y) (1 / x - 1 / y) ^ 2

  dd <- dim(p1$o_mat)[1]
  cost <- matrix(NA_real_, dd, dd)
  posneg <- matrix(NA_integer_, dd, dd)
  for (cc2 in seq_len(dd)) {
    for (cc1 in seq_len(dd)) {
      pos <- sum((p2$o_mat[, cc2] - p1$o_mat[, cc1]) ^ 2) +
        sigdistsq(p2$sigltil[cc2], p1$sigltil[cc1]) +
        sigdistsq(p2$sigrtil[cc2], p1$sigrtil[cc1])
      neg <- sum((p2$o_mat[, cc2] + p1$o_mat[, cc1]) ^ 2) +
        sigdistsq(p2$sigltil[cc2], p1$sigrtil[cc1]) +
        sigdistsq(p2$sigrtil[cc2], p1$sigltil[cc1])
      cost[cc2, cc1] <- min(pos, neg)
      posneg[cc2, cc1] <- if (neg < pos) -1L else 1L
    }
  }

  perm <- as.integer(.solve_lsap_cpp(cost))
  sq_dist_other_params <- sum((p1$mu - p2$mu) ^ 2) +
    (p1$ctil - p2$ctil) ^ 2 +
    (p1$pd - p2$pd) ^ 2

  pairing <- cbind(seq_len(nrow(cost)), perm)
  distance <- sqrt(sum(cost[pairing]) + sq_dist_other_params)
  if (!give_closest_rep) {
    return(distance)
  }

  perm_inv <- order(as.numeric(perm))
  posnegs <- posneg[pairing]
  posnegs <- posnegs[perm_inv]
  flip <- (posnegs == -1)
  rep_ec <- convert_equivalence_class(p1, flip = flip, perm = perm)
  rep_out <- list(
    mu = p1$mu,
    sigltil = rep_ec$sigltil,
    sigrtil = rep_ec$sigrtil,
    ctil = p1$ctil,
    pd = p1$pd,
    o_mat = rep_ec$o_mat
  )
  list(distance = distance, representative = rep_out)
}
