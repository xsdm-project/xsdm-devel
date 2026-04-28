#' Back-compatibility alias for `like_neg_ltsgr`
#'
#' Identical to \code{\link{like_neg_ltsgr}}. Kept as an unexported alias
#' so that any internal code, tests, or downstream scripts that still
#' reference \code{xsdm:::like_neg_ltsgr_cpp} continue to work after the
#' canonical rename. New code should call \code{like_neg_ltsgr} directly.
#'
#' @keywords internal
#' @noRd
like_neg_ltsgr_cpp <- function(env_dat,
                               mu,
                               sigltil,
                               sigrtil,
                               o_mat,
                               num_threads = RcppParallel::defaultNumThreads()) {
  like_neg_ltsgr(
    env_dat     = env_dat,
    mu          = mu,
    sigltil     = sigltil,
    sigrtil     = sigrtil,
    o_mat       = o_mat,
    num_threads = num_threads
  )
}
