#' Validate the environmental data array
#'
#' Checks that \code{env_dat} is a 3-dimensional array with dimensions
#' \code{n_loc x n_time x p} (number of locations x time-series length x number
#' of environmental variables) and contains no missing values.  Throws an
#' informative error if any condition is violated.
#'
#' @param env_dat The object to validate.
#' @param name Character scalar.  Variable name shown in error messages.
#'   Defaults to \code{"env_dat"}.
#'
#' @returns \code{env_dat} invisibly (allows use in a pipe).
#' @keywords internal
check_env_array <- function(env_dat, name = "env_dat") {
  checkmate::assert_array(
    env_dat,
    d = 3,
    any.missing = FALSE,
    .var.name = name
  )
  invisible(env_dat)
}
