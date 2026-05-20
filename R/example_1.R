#' Consolidated example data for the xsdm package
#'
#' A named list containing all example datasets used in the package's
#' documentation and examples.
#'
#' @format A list of 9 objects:
#' \describe{
#'   \item{\code{par_vec}}{Named numeric vector of length 9. Math-scale
#'     parameters for a 2-variable model (p = 2). Canonical names:
#'     \code{mu1}, \code{mu2}, \code{sigltil1}, \code{sigltil2},
#'     \code{sigrtil1}, \code{sigrtil2}, \code{ctil}, \code{pd},
#'     \code{o_par1}.}
#'   \item{\code{bio01}}{A packed \code{SpatRaster} (use
#'     \code{terra::unwrap()}) with 128 x 123 cells and 39 layers.
#'     Annual average temperature (bio1) for 1980-2018, CHELSA 2.1 data,
#'     centred on southern New Mexico, USA.}
#'   \item{\code{bio12}}{A packed \code{SpatRaster} (use
#'     \code{terra::unwrap()}) with 128 x 123 cells and 39 layers.
#'     Annual precipitation (bio12) for the same region and time period.}
#'   \item{\code{env_array}}{A 3-D numeric array with dimensions
#'     4000 (locations) x 39 (time) x 2 (variables). Contains the
#'     environmental data (bio1 and bio12, both divided by 100) extracted
#'     from the rasters at the \code{occ_df} locations.}
#'   \item{\code{occ_df}}{A data frame with 4000 rows and 4 columns:
#'     \code{name} (character), \code{lon} (longitude), \code{lat}
#'     (latitude), \code{presence} (0/1). Occurrence records for a
#'     virtual species.  Use \code{occ_df$presence} wherever a binary
#'     occurrence vector is needed.}
#'   \item{\code{par_list}}{A list of biological-scale parameters
#'     (the "true" parameter set used to generate the virtual species).
#'     Contains \code{mu}, \code{sigltil}, \code{sigrtil}, \code{ctil},
#'     \code{pd}, \code{o_mat}.}
#'   \item{\code{optim_par_vec}}{A named numeric vector of length 9.
#'     Math-scale MLE estimates.  Convert to the biological scale with
#'     \code{math_to_bio(example_1$optim_par_vec)}.}
#'   \item{\code{optim_par_vec_equivalent}}{A named numeric vector of
#'     length 9. A different math-scale representation that belongs to
#'     the same equivalence class as \code{optim_par_vec}. Used to test
#'     \code{dist_between_params()}.}
#'   \item{\code{par_table}}{A data.frame with 9 columns (one per
#'     math-scale parameter) and 100 rows of parameter combinations.}
#' }
#'
#' @details
#' All rasters (\code{bio01}, \code{bio12}) are stored as packed
#' \code{SpatRaster} objects to reduce package size. Before using them,
#' unpack with \code{terra::unwrap()}, e.g.:
#' \code{bio1 <- terra::unwrap(example_1$bio01)}.
#'
#' The environmental data are originally from CHELSA v2.1
#' (\url{https://www.chelsa-climate.org/}). The virtual species was
#' generated from the parameters in \code{par_list}.
#'
#' Three former convenience slots can be derived from the remaining
#' objects and were removed to keep the dataset minimal:
#' \itemize{
#'   \item \code{occ_vec}: use \code{example_1$occ_df$presence}
#'   \item \code{optim_par_list}: use
#'     \code{math_to_bio(example_1$optim_par_vec)}
#'   \item \code{par_vec_vsp}: was identical to \code{par_vec}
#' }
#'
#' @examples
#' # Access the list
#' names(example_1)
#'
#' # Unpack a raster
#' \donttest{
#' bio1 <- terra::unwrap(example_1$bio01)
#' }
#'
#' # Use a parameter set
#' math_to_bio(example_1$par_vec)
#'
#' @source Berti et al., 2025 (\doi{10.1101/2024.10.30.621023})
"example_1"
