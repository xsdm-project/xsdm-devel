#' Consolidated example data for the xsdm.
#' This is environmental data array  an occurence presence absence vector
#' Blarina carolinensis
#' A named list containing all example datasets used in the package's
#' documentation and examples.
#'
#' @format A list of 2 objects:
#' \describe{
#'   \item{\code{env_array}}{A 3‑D numeric array with dimensions
#'     `1156 (locations) × 39 (time) × 2 (variables)`. Contains the
#'     environmental data (bio1 and bio12) extracted from the rasters
#'     for all locations.}
#'   \item{\code{occ_vec}}{An integer vector of length 1156 Binary
#'     presence/absence (0/1) for the same locations as `env_array`.}
#' }
#'
#' @examples
#' # Access the list
#' names(example_3)
#'
#' @source Berti e al,  2024 (<https://doi.org/10.1101/2024.10.30.621023>)
"example_3"