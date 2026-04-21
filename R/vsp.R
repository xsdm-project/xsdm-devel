#' Generate a virtual species probability map
#'
#' Creates a virtual species probability-of-detection map based on environmental
#' time-series data and a set of species-specific parameters.
#'
#' @param env_data A named list of time-series raster objects (e.g., bioclimatic
#'  variables). Each element should be a `SpatRaster` or similar object from the
#'  `terra` package.
#' @param param_list A named list of parameters required by `log_prob_detect()`.
#'   Must include `mu`, `sigltil`, `sigrtil`, `ctil`, `pd`, and `o_mat`.
#' @param return_raster Logical. If `TRUE`, returns a `SpatRaster` object with
#' probabilities. If `FALSE`, returns a tibble with columns `x`, `y`, and
#' `probs`.
#'
#' @return Either:
#'   * A `SpatRaster` object (if `return_raster = TRUE`), or
#'   * A tibble with coordinates and probability values
#'   (if `return_raster = FALSE`).
#'
#' @details
#' Internally, the function:
#' \enumerate{
#'   \item Converts the list of rasters into an array using `env_data_array()`.
#'   \item Applies `log_prob_detect()` with the provided parameters.
#'   \item Exponentiates the log-probabilities to obtain detection
#'   probabilities.
#' }
#'
#' @examples
#' # Example using chelsa dataset preloaded in xsdm package:
#' bio1_ts <- terra::unwrap(example_1_bio01)
#' bio12_ts <- terra::unwrap(example_1_bio12)
#' bio1_ts <- bio1_ts / 100
#' bio12_ts <- bio12_ts / 100
#' env_data <- list(bio1 = bio1_ts, bio12 = bio12_ts)
#' param_list_example <- list(
#'   mu = c(14, 6.5),
#'   sigltil = c(0.46, 1.08),
#'   sigrtil = c(0.105, 0.9),
#'   ctil = -18.14,
#'   pd = 0.89,
#'   o_mat = matrix(c(-0.18, 0.983, -0.983, -0.18), nrow = 2, ncol = 2)
#' )
#' vsp(env_data, param_list_example, return_raster = TRUE)
#'
#' @seealso [env_data_array()], [log_prob_detect()], [terra::rast()]
#' @export
vsp <- function(env_data, param_list, return_raster) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required for this function. Install it with install.packages('terra').")
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required for this function. Install it with install.packages('tibble').")
  }
  # Validate inputs using checkmate
  checkmate::assert_list(env_data,
    types = "SpatRaster",
    min.len = 1,
    any.missing = FALSE
  )
  checkmate::assert_list(param_list, names = "unique", any.missing = FALSE)
  checkmate::assert_true(
    all(
      c("mu", "sigltil", "sigrtil", "ctil", "pd", "o_mat") %in% names(param_list)
    ),
    .var.name = "param_list must contain: mu, sigltil, sigrtil, ctil, pd, o_mat"
  )
  checkmate::assert_flag(return_raster)
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required for vsp(). Install it with: install.packages('terra')")
  }

  # Convert environmental data to array
  env_m <- env_data_array(env_data)

  # Generate function for probability calculation
  f <- function(env_) {
    function(mu, sigltil, sigrtil, o_mat, ctil, pd) {
      log_prob_detect(env_, mu, sigltil, sigrtil, o_mat, ctil, pd)
    }
  }

  # Apply the function to the environmental array

  f_par <- f(env_m)

  # Extract coordinates and CRS
  coords <- terra::crds(env_data[[1]])
  crs_val <- terra::crs(env_data[[1]])

  # Compute probabilities
  probs <- suppressWarnings(do.call(f_par, args = param_list))
  probs <- exp(probs)

  # Return result
  if (!return_raster) {
    data.frame(coords, probs) |> tibble::as_tibble()
  } else {
    data.frame(coords, probs) |>
      terra::rast(crs = crs_val)
  }
}
