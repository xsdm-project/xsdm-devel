#' Get an array of environmental data from presence-absence points.
#'
#' @param env_data List of environmental variables time series stacks (each a
#'   SpatRaster with multiple layers).
#' @param occ Occurrence data frame. Should contain columns "name", "x", "y",
#'   "presence". If NULL, returns data for all raster cells.
#'
#' @return A 3D array of dimensions M (points or cells) × N (time steps) × P
#'   (environmental variables). The first dimension has no dimnames; the second
#'   is named "time" with layer names from the first raster; the third is named
#'   "var" with the names of `env_data`.
#' @export
#'
#' @examples
#' bio1_ts <- terra::unwrap(examples$bio01)
#' bio12_ts <- terra::unwrap(examples$bio12)
#' env_data <- list(bio1 = bio1_ts, bio12 = bio12_ts)
#' occ <- examples$occ_df[1:5, ]
#' # Return array correspoding to each presence absence provided
#' env_data_array(env_data, occ)
#' # Return all the environmental in the rasters
#' env_data_array(env_data, occ)
env_data_array <- function(env_data, occ = NULL) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required for this function. Install it with install.packages('terra').")
  }
  checkmate::assert_list(env_data,
    any.missing = FALSE,
    null.ok = FALSE,
    min.len = 1
  )
  
  checkmate::assert_data_frame(occ, any.missing = FALSE, null.ok = TRUE)
  
  if (!is.null(occ)) {
    checkmate::assert_names(
      names(occ),
      must.include = c("name", "x", "y", "presence")
    )
  }
  
  n_vars <- length(env_data)
  time_names <- names(env_data[[1]])
  var_names <- names(env_data)
  
  # Create spatial points if occurrences are provided
  if (!is.null(occ)) {
    checkmate::assert_names(names(occ),
      must.include = c("name", "x", "y", "presence")
    )
    pts <- terra::vect(occ, geom = c("x", "y"))
  } else {
    pts <- NULL
  }
  # Helper to extract data: returns matrix M x N (or C x N)
  extract_one <- function(r, pts) {
    if (!is.null(pts)) {
      as.matrix(terra::extract(r, pts, cell = FALSE, ID = FALSE))
    } else {
      as.matrix(terra::as.data.frame(r))
    }
  }
  if (n_vars == 1) {
    mat <- extract_one(env_data[[1]], pts)
    arr <- array(mat, dim = c(nrow(mat), ncol(mat), 1))
  } else {
    extracted <- lapply(env_data, extract_one, pts = pts)
    d1 <- nrow(extracted[[1]])
    d2 <- ncol(extracted[[1]])
    # Safety check: all matrices must have identical dimensions
    if (
      !all(vapply(extracted, nrow, integer(1)) == d1) ||
        !all(vapply(extracted, ncol, integer(1)) == d2)) {
      stop("Extracted matrices have inconsistent dimensions.")
    }
    arr <- array(unlist(extracted), dim = c(d1, d2, n_vars))
  }
  dimnames(arr) <- list(NULL, time = time_names, var = var_names)
  arr
}
