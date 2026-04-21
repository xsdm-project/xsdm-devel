#' Provides a measure of distance in parameter space between two sets of parameters.
#'
#' Computes a distance in parameter space between two sets of parameters. The
#' o_mat, sigltil, and sigrtil parameters are only determined up to equivalence
#' class. This function runs through all representatives of the equivalence class
#' of p1 and computes the distance to p2 for each, returning the minimum.
#' Optionally also returns the representative that gives that minimum. Distance
#' is reckoned using sum-squared errors, on the biological scale, except for
#' sigltil and sigrtil which are inverted before measuring distance because that
#' is the scale it makes most sense to measure distance for those parameters.
#'
#' @param p1 First set of parameters. Named list with entries mu, sigltil, sigrtil,
#' ctil, pd, and o_mat.
#' @param p2 Second set of parameters, same format as the first
#' @param GiveClosestRep TRUE if you also want the member of the equivalence class
#' of p1 which is closest to p2. Default FALSE.
#'
#' @return If \code{GiveClosestRep} is FALSE then a single number which is the
#' distance. Otherwise a list with elements for the distance and for the
#' representative which gives that distance.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # This function is internal and mainly used for testing.
#' # Use the exported `dist_between_params()` for general purposes.
#'
#' # Create two equivalent parameter sets (p = 2)
#' p1 <- list(
#'   mu = c(10, 20),
#'   sigltil = c(1.5, 2.0),
#'   sigrtil = c(2.5, 1.0),
#'   ctil = -5,
#'   pd = 0.8,
#'   o_mat = diag(2)
#' )
#'
#' # Build an equivalent p2 by swapping columns (2,1) and flipping the sign
#' # of the first column.
#' perm <- c(2, 1)
#' flip <- 1
#'
#' o_mat_swapped <- p1$o_mat[, perm]
#' o_mat_flipped <- o_mat_swapped
#' o_mat_flipped[, flip] <- -o_mat_flipped[, flip]
#'
#' sigltil_swapped <- p1$sigltil[perm]
#' sigrtil_swapped <- p1$sigrtil[perm]
#' sigltil_final <- sigltil_swapped
#' sigrtil_final <- sigrtil_swapped
#' sigltil_final[flip] <- sigrtil_swapped[flip]
#' sigrtil_final[flip] <- sigltil_swapped[flip]
#'
#' p2 <- list(
#'   mu = p1$mu,
#'   sigltil = sigltil_final,
#'   sigrtil = sigrtil_final,
#'   ctil = p1$ctil,
#'   pd = p1$pd,
#'   o_mat = o_mat_flipped
#' )
#'
#' # Distance should be zero (or near machine precision)
#' xsdm:::distance_between_params(p1, p2)
#' 
#' # Using the package's example data.  Distance should be zero
#' # (or near machine precision)
#' bio1 <- math_to_bio(example_1_optim_param_vector)
#' bio2 <- math_to_bio(example_1_optim_param_vector_equivalent)
#' xsdm:::distance_between_params(bio1, bio2)
#' }
distance_between_params <- function(p1, p2, GiveClosestRep = FALSE) {
  # No error checking of inputs because this is not an exported function, it is
  # only used to test distance_between_params_hungarian.
  
  # Warning
  dd <- dim(p1$o_mat)[1]
  if (dd > 5) {
    warning("In distance_between_params: this function is very slow for large
            numbers of environmental variables")
  }
  
  # compute square distance for all parameters except o_mat, sigltil, and sigrtil
  sq_dist_other_params <- sum((p1$mu - p2$mu)^2) + (p1$ctil - p2$ctil)^2 +
    (p1$pd - p2$pd)^2
  
  # get lists of all the sign flip patterns and all the permutations which are
  # to be implemented in all possible combinations
  all_perms <- permutations(dd, dd, 1:dd)
  all_flips <- as.matrix(expand.grid(rep(list(c(0, 1)), dd)))
  
  # now iterate through representatives of the equivalence class of p1 and
  # compute distances to p2
  bestsqdist <- Inf
  bestparams <- NA
  for (fc in 1:(dim(all_flips)[1]))
  {
    for (pc in 1:(dim(all_perms)[1]))
    {
      p1_ec <- convert_equivalence_class(p1, all_flips[fc, ], all_perms[pc, ])
      thissqdist <- sum((p1_ec$o_mat - p2$o_mat)^2) +
        sum((1 / p1_ec$sigltil - 1 / p2$sigltil)^2) +
        sum((1 / p1_ec$sigrtil - 1 / p2$sigrtil)^2)
      if (thissqdist < bestsqdist) {
        bestsqdist <- thissqdist
        bestparams <- p1_ec
      }
    }
  }
  
  # now return
  if (!GiveClosestRep) {
    return(sqrt(bestsqdist + sq_dist_other_params))
  }
  rep <- list(
    mu = p1$mu, sigltil = bestparams$sigltil, sigrtil = bestparams$sigrtil,
    ctil = p1$ctil, pd = p1$pd, o_mat = bestparams$o_mat
  )
  res <- list(
    distance = sqrt(bestsqdist + sq_dist_other_params),
    representative = rep
  )
  return(res)
}