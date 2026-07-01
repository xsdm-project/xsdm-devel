#' Generate a virtual species probability map with presence/absence sampling
#'
#' Creates a virtual species probability-of-detection map based on environmental
#' time-series data and a set of species-specific parameters, then samples
#' presence/absence points based on a user-defined probability threshold.
#'
#' @param param_list A named list of biological‑scale parameters required by
#'   `log_prob_detect()`. Must include `mu`, `sigltil`, `sigrtil`, `ctil`, `pd`,
#'   and `o_mat`. Values like `sigltil`/`sigrtil` can be `Inf`.
#' @param env_data A named list of time‑series raster objects (e.g., from the
#'   `terra` package). Each element must be a `SpatRaster` with the same
#'   geometry and number of layers.
#' @param size_presence Integer. Number of sample points to draw from cells where
#'   the detection probability **exceeds** `threshold`.
#' @param size_absence Integer. Number of sample points to draw from cells where
#'   the detection probability is **less than or equal to** `threshold`.
#' @param threshold Numeric in `[0, 1]`. Probability cutoff used to distinguish
#'   presence vs. absence sampling areas. Default `0.5`.
#'
#' @return A tibble with columns `lon`, `lat`, `presence` (0/1), where each row
#'   corresponds to a sampled point. The presence/absence is drawn from a
#'   binomial distribution using the habitat suitability value as the success
#'   probability.
#'
#' @details
#' Internally the function:
#' \enumerate{
#'   \item Computes a habitat suitability raster using `habitat_suitability()`.
#'   \item Splits the raster into two layers based on `threshold`:
#'         cells with prob > threshold (presence pool) and \eqn{\le} threshold (absence pool).
#'   \item Samples `size_presence` and `size_absence` points from each pool
#'         (without replacement), with probabilities proportional to the suitability value.
#'   \item Generates a binomial outcome for each sampled point using its suitability
#'         as the probability of success.
#' }
#'
#' @seealso [habitat_suitability()], [log_prob_detect()], [terra::spatSample()]
#' @export
#'
#' @examples
#' \donttest{
#' data("example_1", package = "xsdm")
#' bio1_ts  <- terra::unwrap(example_1$bio01) / 100
#' bio12_ts <- terra::unwrap(example_1$bio12) / 100
#' env_data <- list(bio1 = bio1_ts, bio12 = bio12_ts)
#'
#' vsp(
#'   param_list    = example_1$true_par_list,
#'   env_data      = env_data,
#'   size_presence = 100,
#'   size_absence  = 100,
#'   threshold     = 0.7
#' )
#' }
vsp <- function(param_list, env_data, size_presence, size_absence, threshold = 0.5) {
  # Input validation --------------------------------------------------------
  checkmate::assert_list(env_data,
                         types = "SpatRaster",
                         min.len = 1,
                         any.missing = FALSE
  )
  checkmate::assert_list(param_list,
                         names = "unique",
                         any.missing = FALSE
  )
  checkmate::assert_true(
    all(c("mu", "sigltil", "sigrtil", "ctil", "pd", "o_mat") %in% names(param_list)),
    .var.name = "param_list must contain: mu, sigltil, sigrtil, ctil, pd, o_mat"
  )
  checkmate::assert_count(size_presence, positive = TRUE)
  checkmate::assert_count(size_absence, positive = TRUE)
  checkmate::assert_number(threshold, lower = 0, upper = 1, finite = TRUE)
  
  # Compute habitat suitability raster --------------------------------------
  r <- habitat_suitability(
    param_list = param_list,
    env_list   = env_data,
    return_prob = TRUE
  )
  
  # Split cells by threshold using cell indices ------------------------------
  vals <- terra::values(r)
  
  # Helper to sample safely and generate presence/absence
  sample_group <- function(cell_indices, sample_size, prob_type) {
    # cell_indices: indices of cells in the target group (presence or absence)
    # prob_type: "presence" or "absence" – only used for warning messages
    n_cells <- length(cell_indices)
    if (n_cells == 0) {
      if (sample_size > 0) {
        warning(sprintf("No cells available for %s sampling (threshold = %f). Returning empty data frame.",
                        prob_type, threshold))
      }
      return(data.frame(x = numeric(0), y = numeric(0), prob = numeric(0)))
    }
    if (sample_size > n_cells) {
      warning(sprintf("Requested sample size (%d) for %s exceeds available cells (%d). Sampling all cells without replacement.",
                      sample_size, prob_type, n_cells))
      sample_size <- n_cells
    }
    # Sample cell indices uniformly without replacement. Suitability drives the
    # presence/absence outcome downstream via rbinom(), so it must not also bias
    # the spatial sampling. Index into cell_indices via sample.int() to avoid the
    # sample() footgun when cell_indices has length 1.
    chosen <- cell_indices[sample.int(n_cells, size = sample_size, replace = FALSE)]
    coords <- terra::xyFromCell(r, chosen)
    data.frame(x = coords[, "x"], y = coords[, "y"], prob = vals[chosen])
  }
  
  # Sample presence and absence groups
  presence_pts <- sample_group(which(vals > threshold),  size_presence, "presence")
  absence_pts  <- sample_group(which(vals <= threshold), size_absence,  "absence")
  
  # Generate binomial outcomes
  generate_binom <- function(pts) {
    if (nrow(pts) == 0) return(pts)
    pts$occurrence <- stats::rbinom(nrow(pts), size = 1, prob = pts$prob)
    pts$prob <- NULL  # remove the probability column
    return(pts)
  }
  
  presence_pts <- generate_binom(presence_pts)
  absence_pts  <- generate_binom(absence_pts)
  
  # Combine and rename columns
  occ_df <- rbind(presence_pts, absence_pts)
  if (nrow(occ_df) > 0) {
    names(occ_df) <- c("lon", "lat", "presence")
  } else {
    # Create an empty data frame with correct column names
    occ_df <- data.frame(lon = numeric(0), lat = numeric(0), presence = integer(0))
  }
  
  # Return as tibble
  tibble::as_tibble(occ_df)
}