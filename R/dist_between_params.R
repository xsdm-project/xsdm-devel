#' Provides a measure of distance in parameter space between two sets of
#' parameters, using the Hungarian algorithm.
#'
#' Computes a distance in parameter space between two sets of parameters. The
#' o_mat, sigltil, and sigrtil parameters are only determined up to equivalence
#' class. This function runs through all representatives of the equivalence
#' class of p1 and computes the distance to p2 for each, returning the minimum.
#' Optionally also returns the equivalence class representative that gives that
#' minimum. Only it does not literally compute all that, it uses the Hungarian
#' algorithm to speed it up. Distance is reckoned using sum-squared errors, on
#' the biological scale, except for sigltil and sigrtil which are inverted
#' before measuring distance because that is the scale it makes most sense to
#' measure distance for those parameters.
#'
#' @param p1 First set of parameters. May be math- or biological-scale parameters,
#' following the associated formatting requirements in each case.
#' @param p2 Second set of parameters, same format options as the first
#' @param mask Same format as the \code{mask} argument to \code{\link{loglik_math}}
#' and \code{\link{start_parms}}. Ignored if both \code{p1} and \code{p2} are on
#' the biological scale. Otherwise the parameters contained in \code{mask} must
#' exactly complement those contained in whichever of \code{p1} and \code{p2} are on the
#' math scale.
#' @param give_closest_rep TRUE if you also want the member of the equivalence class
#' of p1 which is closest to p2. Gives the answer on the biological scale.
#' Default FALSE.
#'
#' @return If \code{give_closest_rep} is FALSE then a single number which is the
#' distance. Otherwise a list with elements for the distance and for the
#' representative which gives that distance.
#' @export
#'
#' @examples
#' # Using list in biological scale
#' par_list <- math_to_bio(examples$optim_par_vec)
#' par_list_equivalent <- math_to_bio(examples$optim_par_vec_equivalent)
#' dist_between_params(
#'   p1 = par_list,
#'   p2 = par_list_equivalent
#' )
#' # Using vectors in math scale
#' dist_between_params(
#'   p1 = examples$optim_par_vec,
#'   p2 = examples$optim_par_vec_equivalent
#' )
dist_between_params <- function(p1, p2, mask = NULL, give_closest_rep = FALSE) {
  # --- Check inputs, basic ---
  checkmate::assert_flag(give_closest_rep)
  if (!((is.null(mask)) || (is.numeric(mask) && (sum(is.na(mask)) == 0)))) {
    stop("mask must be NULL or a numeric vector with no missing values")
  }
  checkmate::assert_true(is.list(p1) || is.numeric(p1))
  checkmate::assert_true(is.list(p2) || is.numeric(p2))
  
  # --- now convert to the biological scale if needed, and also do more input checking ---
  if (is.numeric(p1)) {
    # get names in p1 and names in mask and expect them to complement each other
    checkmate::assert_integerish(sqrt(9 + 8 * (length(p1) + length(mask)))) # check that the length is allowed
    p <- round((-5 + sqrt(9 + 8 * (length(p1) + length(mask)))) / 2)
    nmask <- names(mask)
    np1 <- names(p1)
    checkmate::assert_true(setequal(c(names(mask), names(p1)), names(make_mask_names(p))))
    
    # checks on which values can be infinite
    checkmate::assert_numeric(p1, any.missing = FALSE, finite = TRUE) # p1 has to have all finite values
    allnames <- names(make_mask_names(p))
    maskinfnames <- names(mask[is.infinite(mask)])
    checkmate::assert_true(all(names(mask[is.infinite(mask)]) %in%
                                 c(allnames[grepl("^sig", allnames)], "pd"))) # mask can only have infinite values for pd and the sigs
    
    # now convert to bio scale
    p1 <- math_to_bio(create_param_vector_masked(p1, mask, p))
  }
  if (is.numeric(p2)) {
    # get names in p2 and names in mask and expect them to complement each other
    checkmate::assert_integerish(sqrt(9 + 8 * (length(p2) + length(mask)))) # check that the length is allowed
    p <- round((-5 + sqrt(9 + 8 * (length(p2) + length(mask)))) / 2)
    nmask <- names(mask)
    np2 <- names(p2)
    checkmate::assert_true(setequal(c(names(mask), names(p2)), names(make_mask_names(p))))
    
    # checks on which values can be infinite
    checkmate::assert_numeric(p2, any.missing = FALSE, finite = TRUE) # p2 has to have all finite values
    allnames <- names(make_mask_names(p))
    maskinfnames <- names(mask[is.infinite(mask)])
    checkmate::assert_true(all(names(mask[is.infinite(mask)]) %in%
                                 c(allnames[grepl("^sig", allnames)], "pd"))) # mask can only have infinite values for pd and the sigs
    
    # now convert to bio scale
    p2 <- math_to_bio(create_param_vector_masked(p2, mask, p))
  }
  
  # --- now do biological-scale checks ---
  checkmate::assert_names(names(p1), must.include = c("mu", "sigltil", "sigrtil", "ctil", "pd", "o_mat"))
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
  
  
  checkmate::assert_names(names(p2), must.include = c("mu", "sigltil", "sigrtil", "ctil", "pd", "o_mat"))
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
  
  # -- distance to be used between log-scale/math-scale entries of sigxtil,
  # currently just squared difference
  sigdistsq <- function(x, y) {
    return((1 / x - 1 / y)^2)
  }
  
  # --- get the cost matrix ---
  dd <- dim(p1$o_mat)[1]
  
  cost <- matrix(NA, dd, dd)
  posneg <- matrix(NA, dd, dd)
  for (cc2 in 1:dd)
  {
    for (cc1 in 1:dd)
    {
      pos <- sum((p2$o_mat[, cc2] - p1$o_mat[, cc1])^2) + sigdistsq(p2$sigltil[cc2], p1$sigltil[cc1]) + sigdistsq(p2$sigrtil[cc2], p1$sigrtil[cc1])
      neg <- sum((p2$o_mat[, cc2] + p1$o_mat[, cc1])^2) + sigdistsq(p2$sigltil[cc2], p1$sigrtil[cc1]) + sigdistsq(p2$sigrtil[cc2], p1$sigltil[cc1])
      cost[cc2, cc1] <- min(pos, neg)
      if (pos < neg) {
        posneg[cc2, cc1] <- 1
      } else if (neg < pos) {
        posneg[cc2, cc1] <- -1
      } else {
        posneg[cc2, cc1] <- 1 # tie: either sign works, pick positive
      }
    }
  }
  
  # --- apply the Hungarian algorithm ---
  
  # got a special object type, see help for solve_LSAP
  perm <- clue::solve_LSAP(cost)
  
  # compute square distance for all parameters except o_mat, sigltil, and sigrtil
  sq_dist_other_params <- sum((p1$mu - p2$mu)^2) +
    (p1$ctil - p2$ctil)^2 +
    (p1$pd - p2$pd)^2
  
  # --- now get the answer(s) and return
  
  pairing <- cbind(seq_len(nrow(cost)), perm)
  costs <- cost[pairing]
  distance <- sqrt(sum(costs) + sq_dist_other_params)
  if (!give_closest_rep) {
    return(distance)
  }
  
  perm_inv <- order(as.numeric(perm))
  posnegs <- posneg[pairing]
  posnegs <- posnegs[perm_inv]
  flip <- (posnegs == -1)
  rep <- convert_equivalence_class(p1, flip = flip, perm = perm)
  rep <- list(
    mu = p1$mu,
    sigltil = rep$sigltil,
    sigrtil = rep$sigrtil,
    ctil = p1$ctil,
    pd = p1$pd,
    o_mat = rep$o_mat
  )
  res <- list(distance = distance, representative = rep)
  return(res)
}
