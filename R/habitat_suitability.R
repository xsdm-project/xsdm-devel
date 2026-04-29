#' Tiled habitat-suitability map from environmental raster stacks
#'
#' Evaluates the log detection probability (or its exponential, the
#' probability of detection) for every cell of a list of multi-layer
#' \code{terra::SpatRaster} objects, processing the inputs in
#' memory-bounded blocks so that arbitrarily large grids can be handled
#' without loading the entire dataset into R memory.  Each block is
#' forwarded to \code{\link{log_prob_detect_cpp}}, the xtensor-backed C++
#' kernel that consolidates the
#' \code{like_neg_ltsgr() -> like_ltsg()} call chain.
#'
#' @param param_list A named list of biological-scale parameters.  Must
#'   contain \code{mu}, \code{sigltil}, \code{sigrtil}, \code{o_mat},
#'   \code{ctil} and \code{pd}.  See \code{\link{log_prob_detect}} for
#'   details of each element.
#' @param env_list A list of \code{\link[terra]{SpatRaster}} objects, one
#'   per environmental variable.  Each raster must have the same number
#'   of layers (time steps) and identical spatial geometry (extent,
#'   resolution, CRS).  Minimum length 1.
#' @param output Character scalar.  File path for the output GeoTIFF.
#'   The empty string \code{""} (default) creates an in-memory
#'   \code{SpatRaster}.
#' @param overwrite Logical scalar.  If \code{TRUE}, an existing file at
#'   \code{output} is overwritten.  Default \code{FALSE}.
#' @param return_prob Logical scalar.  If \code{TRUE} (default), the
#'   output cell values are probabilities of detection (range
#'   \eqn{[0, 1]}).  If \code{FALSE}, the cell values are
#'   log-probabilities (range \eqn{(-\infty, 0]}).
#' @param threads Integer scalar.  Number of parallel threads forwarded
#'   to \code{\link{log_prob_detect_cpp}}.  Use \code{0} (default) to let
#'   \pkg{RcppParallel} pick the number of threads automatically.
#' @param wopt List.  Additional write options forwarded to
#'   \code{\link[terra]{writeStart}}.  Default \code{list()}.
#'
#' @return A \code{SpatRaster} with one layer named either
#'   \code{"habitat_suitability"} (when \code{return_prob = TRUE}) or
#'   \code{"log_prob_detect"} (when \code{return_prob = FALSE}).  The
#'   raster is returned invisibly when \code{output != ""}.
#'
#' @details
#' Internally the function uses \pkg{terra}'s streaming block-loop API:
#' \enumerate{
#'   \item \code{\link[terra]{readStart}} is called on every raster in
#'     \code{env_list}.
#'   \item \code{\link[terra]{writeStart}} is called on the output
#'     raster, which returns a block schedule chosen by terra's memory
#'     manager.
#'   \item For each block, \code{\link[terra]{readValues}} reads a
#'     horizontal strip from every input raster into a matrix; the
#'     strips are packed into a flat column-major vector and passed to
#'     \code{\link{log_prob_detect_cpp}}.  Cells that are NA in any
#'     variable or time step are masked out and re-inserted as NA in the
#'     output.
#'   \item \code{\link[terra]{writeValues}} writes the per-cell results.
#'   \item \code{\link[terra]{readStop}} and
#'     \code{\link[terra]{writeStop}} are called via \code{\link{on.exit}}
#'     to ensure file handles are released even if an error occurs.
#' }
#' At most one block of pixels is held in R memory at any time, making
#' the function suitable for continental or global rasters.
#'
#' @seealso \code{\link{log_prob_detect_cpp}}, \code{\link{log_prob_detect}},
#'   \code{\link{vsp}}, \code{\link[terra]{writeStart}}
#'
#' @examples
#' \donttest{
#' data("example_1", package = "xsdm")
#' bio01 <- terra::unwrap(example_1$bio01) / 100
#' bio12 <- terra::unwrap(example_1$bio12) / 100
#' env_list <- list(bio01 = bio01, bio12 = bio12)
#' suit <- habitat_suitability(
#'   param_list  = example_1$par_list,
#'   env_list    = env_list,
#'   return_prob = TRUE
#' )
#' suit
#' }
#' @export
habitat_suitability <- function(
    param_list,
    env_list,
    output      = "",
    overwrite   = FALSE,
    return_prob = TRUE,
    threads     = 0L,
    wopt        = list()
) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop(
      "Package 'terra' is required for habitat_suitability(). ",
      "Install it with: install.packages('terra')"
    )
  }

  # ---- Input validation ----
  checkmate::assert_list(param_list, names = "unique", any.missing = FALSE)
  checkmate::assert_true(
    all(c("mu", "sigltil", "sigrtil", "o_mat", "ctil", "pd") %in%
          names(param_list)),
    .var.name = "param_list must contain: mu, sigltil, sigrtil, o_mat, ctil, pd"
  )
  checkmate::assert_list(
    env_list,
    types       = "SpatRaster",
    min.len     = 1L,
    any.missing = FALSE
  )
  checkmate::assert_string(output)
  checkmate::assert_flag(overwrite)
  checkmate::assert_flag(return_prob)
  checkmate::assert_count(threads, positive = FALSE)
  checkmate::assert_list(wopt)

  # ---- Geometry consistency ----
  ref <- env_list[[1L]]
  if (length(env_list) > 1L) {
    for (k in seq.int(2L, length(env_list))) {
      if (!terra::compareGeom(ref, env_list[[k]], stopOnError = FALSE)) {
        stop(
          "All rasters in env_list must have the same geometry ",
          "(extent, resolution, CRS)."
        )
      }
    }
  }

  p  <- length(env_list)
  ts <- terra::nlyr(ref)
  nc <- terra::ncol(ref)

  # ---- Build output SpatRaster ----
  out <- terra::rast(ref, nlyr = 1L)
  names(out) <- if (return_prob) "habitat_suitability" else "log_prob_detect"

  # ---- Open input rasters for streaming reads ----
  for (k in seq_len(p)) terra::readStart(env_list[[k]])
  on.exit({
    for (k in seq_len(p)) {
      tryCatch(terra::readStop(env_list[[k]]), error = function(e) NULL)
    }
  }, add = TRUE)

  # ---- Open output for streaming writes ----
  b <- terra::writeStart(
    out,
    filename  = output,
    overwrite = overwrite,
    wopt      = wopt
  )
  on.exit(tryCatch(terra::writeStop(out), error = function(e) NULL), add = TRUE)

  # ---- Block loop ----
  for (i in seq_len(b$n)) {
    n_tile <- b$nrows[i] * nc

    # Pack each variable's tile into a flat (n_tile x ts x p) column-major vector.
    # Within variable k, the (n_tile x ts) matrix returned by readValues(mat=TRUE)
    # is already in the column-major layout that log_prob_detect_cpp expects:
    # env_dat[l, t, k] = env_dat_vec[l + n_tile*t + n_tile*ts*k].
    env_vec <- numeric(n_tile * ts * p)
    valid   <- rep(TRUE, n_tile)

    for (k in seq_len(p)) {
      tile_k <- terra::readValues(
        env_list[[k]],
        row   = b$row[i],
        nrows = b$nrows[i],
        col   = 1L,
        ncols = nc,
        mat   = TRUE
      )
      offset <- (k - 1L) * n_tile * ts
      env_vec[offset + seq_len(n_tile * ts)] <- as.vector(tile_k)
      # A cell is invalid if it has any NA across time for any variable.
      valid <- valid & !apply(tile_k, 1L, anyNA)
    }

    block_result <- rep(NA_real_, n_tile)
    n_valid <- sum(valid)

    if (n_valid > 0L) {
      if (n_valid < n_tile) {
        # Compact env_vec to only valid rows, preserving (n_valid x ts x p) layout.
        valid_idx <- which(valid)
        compact <- numeric(n_valid * ts * p)
        for (k in seq_len(p)) {
          src_offset <- (k - 1L) * n_tile * ts
          dst_offset <- (k - 1L) * n_valid * ts
          for (t in seq_len(ts)) {
            compact[dst_offset + (t - 1L) * n_valid + seq_len(n_valid)] <-
              env_vec[src_offset + (t - 1L) * n_tile + valid_idx]
          }
        }
        result_valid <- log_prob_detect_cpp(
          env_dat_vec  = compact,
          env_dat_dims = as.integer(c(n_valid, ts, p)),
          mu           = param_list$mu,
          sigltil      = param_list$sigltil,
          sigrtil      = param_list$sigrtil,
          o_mat        = param_list$o_mat,
          ctil         = param_list$ctil,
          pd           = param_list$pd,
          return_prob  = return_prob,
          num_threads  = as.integer(threads)
        )
        block_result[valid] <- result_valid
      } else {
        block_result <- log_prob_detect_cpp(
          env_dat_vec  = env_vec,
          env_dat_dims = as.integer(c(n_tile, ts, p)),
          mu           = param_list$mu,
          sigltil      = param_list$sigltil,
          sigrtil      = param_list$sigrtil,
          o_mat        = param_list$o_mat,
          ctil         = param_list$ctil,
          pd           = param_list$pd,
          return_prob  = return_prob,
          num_threads  = as.integer(threads)
        )
      }
    }

    terra::writeValues(out, block_result, b$row[i], b$nrows[i])
  }

  if (nzchar(output)) invisible(out) else out
}
