#' Create a parameter mask aligned with `make_mask_names()`
#'
#' Constructs a named numeric vector whose names and length match the canonical
#' schema returned by `make_mask_names()`. Optionally fills selected entries
#' from a user-supplied named vector `mask`. The output is intended for use
#' within downstream functions (e.g., `loglik_math`) that need parameters in a
#' fixed order and with standard names: `mu1` \ldots `mup`, `sigltil1`, \ldots,
#' `sigltilp`, `sigrtil1`, \ldots, `sigrtilp`,`ctil`, `pd`, and \code{o_mati}
#' for \code{i} ranging from 1 to \code{(p^2-p)/2}.
#'
#' @param mask Named numeric vector (default `NULL`). Names must be a subset of
#' those produced by `make_mask_names(p)`. Values are inserted into the
#' corresponding positions; all other entries of the output are `NA_real_`.
#' @param p A positive integer representing the number of environmental
#' variables to be used in the xsdm model, i.e., \code{dim(env_dat)[3]} for the
#' \code{env_dat} argument of the \code{loglik_math} function.
#'
#' @returns A named numeric vector of length `num_par(p)` with names in
#' the canonical order given above. Entries are initialized to `NA_real_` except
#' for those provided in `mask`.
#'
#' @examples
#' # Empty mask for p = 2 (all NA values)
#' create_mask(p = 2)
#'
#' # Partially filled mask; unspecified entries remain NA
#' create_mask(mask = c(mu1 = 11, sigltil1 = Inf, pd = 1, ctil = -2), p = 2)
#'
#' # p = 1 has no o_par entries
#' create_mask(mask = c(mu1 = 7, pd = 0.5), p = 1)
#'
#' @seealso [make_mask_names()], [num_par()]
#' @export
create_mask <- function(mask = NULL, p = 1) {
  checkmate::assert_count(p, positive = TRUE)
  # canonical names + NA_real_
  out <- make_mask_names(p)

  if (is.null(mask)) {
    return(out)
  }

  checkmate::assert_numeric(mask, any.missing = FALSE, names = "named")
  bad <- setdiff(names(mask), names(out))
  if (length(bad) > 0) {
    stop(
      "Unexpected parameter name(s): ", paste(bad, collapse = ", "),
      ". Allowed names are: ", paste(names(out), collapse = ", "), "."
    )
  }

  out[names(mask)] <- mask
  out
}
