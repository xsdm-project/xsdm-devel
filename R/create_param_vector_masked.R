#' Create a complete parameter vector with canonical names (no NAs allowed)
#'
#' Builds the full parameter vector for use within `loglik_math` using canonical
#' names from `make_mask_names(p)` via `create_mask(mask = mask, p)`. The output
#' **always contains all canonical names** (length = `num_par(p)`). `mask`
#' is applied first (optional), then `param_vector` overrides (required).
#' The final vector must have **no NA values**; if any entry remains `NA`, the
#' function throws an error listing which names are missing.
#'
#' @param param_vector Named numeric vector (**required**). Names must be a
#' subset of the canonical names returned by `create_mask(mask = NULL, p)`.
#' Values in `param_vector` override those set by `mask`.
#' @param mask Named numeric vector (optional). Names must be canonical and
#'   are applied before `param_vector`.
#' @param p A positive integer representing the number of environmental
#' variables to be used in the xsdm model, i.e., \code{dim(env_dat)[3]} for the
#' \code{env_dat} argument of the \code{loglik_math} function
#'
#' @return A named numeric vector of length `num_par(p)` with **no NA** entries
#' and canonical names in the documented order.
#'
#' @details
#' Canonical names follow `loglik_math` conventions - see Details of that
#' function for the canonical ordering we use.
#'
#' This is a thin R wrapper around the internal C++ implementation
#' \code{xsdm:::.build_canonical_param_vector_cpp}. The pre-port pure-R
#' implementation is preserved internally as
#' \code{xsdm:::create_param_vector_masked_r} for parity testing.
#'
#' @seealso [make_mask_names()], [create_mask()], and `loglik_math` (Details).
#' @export
#' @examples
#' ## --- p = 1 ---
#' p1 <- 1
#' # Canonical names typically: mu1, sigltil1, sigrtil1, ctil, pd
#' pv1 <- c(sigltil1 = 1.0, sigrtil1 = 2.0, ctil = 0.2) # fills remaining slots
#' mask1 <- c(mu1 = -1, pd = 0.5)
#' out1 <- create_param_vector_masked(param_vector = pv1, mask = mask1, p = p1)
#'
#' ## --- p = 2 (includes o_par1) ---
#' p2 <- 2
#' pv2 <- c(
#'   sigltil1 = 1.0, sigltil2 = 1.1, sigrtil1 = 2.0, sigrtil2 = 2.2,
#'   ctil = 0.3, o_par1 = 0.0
#' )
#' mask2 <- c(mu1 = 0.1, mu2 = 0.2, pd = 0.05)
#' out2 <- create_param_vector_masked(param_vector = pv2, mask = mask2, p = p2)
#'
#' ## --- p = 3 (includes o_par1..3) ---
#' p3 <- 3
#' pv3 <- c(
#'   sigltil1 = 1.0, sigltil2 = 1.1, sigltil3 = 1.2,
#'   sigrtil1 = 2.0, sigrtil2 = 2.1, sigrtil3 = 2.2,
#'   ctil = 0.4, o_par1 = -0.2, o_par2 = 0.0, o_par3 = 0.15
#' )
#' mask3 <- c(mu1 = 0.1, mu2 = 0.2, mu3 = 0.3, pd = 0.01)
#' out3 <- create_param_vector_masked(param_vector = pv3, mask = mask3, p = p3)
create_param_vector_masked <- function(param_vector, mask = NULL, p) {
  # R-side validation produces the user-facing checkmate-style error
  # messages historically returned by this function. The actual merge
  # (mask-first, then param_vector overrides into the canonical skeleton)
  # is performed by the internal C++ kernel
  # `xsdm:::.build_canonical_param_vector_cpp`.
  checkmate::assert_count(p, positive = TRUE)
  if (is.null(param_vector)) {
    stop("`param_vector` must not be NULL.")
  }
  checkmate::assert_numeric(param_vector, any.missing = FALSE, names = "named")

  allowed_names <- names(create_mask(mask = NULL, p = p))
  if (is.null(allowed_names)) {
    stop("`create_mask(mask = NULL, p)` must return a named vector;
         names are canonical.")
  }

  bad_pv <- setdiff(names(param_vector), allowed_names)
  if (length(bad_pv) > 0) {
    stop(
      "Unexpected name(s) in `param_vector`: ", paste(bad_pv, collapse = ", "),
      ". Allowed: ", paste(allowed_names, collapse = ", "), "."
    )
  }

  if (!is.null(mask)) {
    checkmate::assert_numeric(mask, any.missing = FALSE, names = "named")
    bad_mask <- setdiff(names(mask), allowed_names)
    if (length(bad_mask) > 0) {
      stop(
        "Unexpected parameter name(s) in `mask`: ",
        paste(bad_mask, collapse = ", "),
        ". Allowed: ", paste(allowed_names, collapse = ", "), "."
      )
    }
    overlap <- intersect(names(param_vector), names(mask))
    if (length(overlap) > 0) {
      stop(
        "`param_vector` and `mask` must be complementary (disjoint). ",
        "Overlapping names: ", paste(overlap, collapse = ", "), "."
      )
    }
  }

  # Detect missing canonical slots up front so we can produce the
  # historical "contains NA values" error message instead of the
  # per-slot message produced by the C++ kernel on the first hit.
  provided <- union(names(param_vector), names(mask))
  missing_names <- setdiff(allowed_names, provided)
  if (length(missing_names) > 0) {
    stop(
      "The final parameter vector contains NA values for: ",
      paste(missing_names, collapse = ", "),
      ". Provide values via `mask` and/or `param_vector` so all positions are
      filled."
    )
  }

  .build_canonical_param_vector_cpp(
    param_vector = as.numeric(param_vector) |>
      stats::setNames(names(param_vector)),
    mask         = if (is.null(mask)) NULL else
      as.numeric(mask) |> stats::setNames(names(mask)),
    p            = as.integer(p)
  )
}
