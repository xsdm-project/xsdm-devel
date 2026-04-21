#' Consolidated example data for the xsdm package
#'
#' A named list containing all example datasets used in the package's
#' documentation and examples.
#'
#' @format A list of 11 objects:
#' \describe{
#'   \item{\code{par_vec}}{Named numeric vector of length 9. Math-scale
#'     parameters for a 2‑variable model (p = 2). Canonical names:
#'     `mu1`, `mu2`, `sigltil1`, `sigltil2`, `sigrtil1`, `sigrtil2`,
#'     `ctil`, `pd`, `o_par1`.}
#'   \item{\code{bio01}}{A packed `SpatRaster` (use `terra::unwrap()`)
#'     with 128 × 123 cells and 39 layers. Annual average temperature
#'     (bio1) for 1980–2018, CHELSA 2.1 data, centred on southern New
#'     Mexico, USA.}
#'   \item{\code{bio12}}{A packed `SpatRaster` (use `terra::unwrap()`)
#'     with 128 × 123 cells and 39 layers. Annual precipitation (bio12)
#'     for the same region and time period.}
#'   \item{\code{env_array}}{A 3‑D numeric array with dimensions
#'     `2000 (locations) × 39 (time) × 2 (variables)`. Contains the
#'     environmental data (bio1 and bio12) extracted from the rasters
#'     for all locations.}
#'   \item{\code{occ_df}}{A data frame with 1000 rows and 4 columns:
#'     `name` (character), `x` (longitude), `y` (latitude), `presence`
#'     (0/1). Occurrence records for the virtual species *Mus virtualis*.}
#'   \item{\code{occ_vec}}{An integer vector of length 2000. Binary
#'     presence/absence (0/1) for the same locations as `env_array`.}
#'   \item{\code{optim_par_list}}{A list of biological‑scale parameters
#'     (the MLE fit for the example). Contains `mu`, `sigltil`, `sigrtil`,
#'     `ctil`, `pd`, `o_mat`.}
#'   \item{\code{optim_par_vec}}{A named numeric vector of length 9.
#'     Math‑scale parameters corresponding to `optim_par_list`.}
#'   \item{\code{optim_par_vec_equivalent}}{A named numeric vector of
#'     length 9. A different math‑scale representation that belongs to
#'     the same equivalence class as `optim_par_vec`. Used to test
#'     `dist_between_params()`.}
#'   \item{\code{par_list}}{A list of biological‑scale parameters
#'     (a “true” parameter set, not necessarily the MLE). Used in
#'     examples of `interpret_parameters()`, `vsp()`, etc.}
#'   \item{\code{par_vec_vsp}}{A named numeric vector of length 9,
#'     math‑scale parameters. Currently an alias for `par_vec` (retained
#'     for backward compatibility).}
#'   \item{\code{par_table}}{A data.frame with 9 columns corresponding to
#'   a parameter in math scale and 100 rows corresponding to 100 combinations.}
#' }
#'
#' @details
#' All rasters (`bio01`, `bio12`) are stored as packed `SpatRaster`
#' objects to reduce package size. Before using them, unpack with
#' `terra::unwrap()`, e.g.:
#' \code{bio1 <- terra::unwrap(examples$bio01)}.
#'
#' The environmental data are originally from CHELSA v2.1
#' (<https://www.chelsa-climate.org/>). The virtual species was
#' generated from the parameters in `par_list`.
#'
#' @examples
#' # Access the list
#' names(examples)
#'
#' # Unpack a raster
#' \donttest{
#' bio1 <- terra::unwrap(examples$bio01)
#' }
#'
#' # Use a parameter set
#' math_to_bio(examples$par_vec)
#'
#' @source Berti e al,  2024 (<https://doi.org/10.1101/2024.10.30.621023>)
"examples"