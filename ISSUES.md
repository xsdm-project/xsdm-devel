# Proposed Feature Issues

---

## Issue 1

**Title:** feat: add log_prob_detect_cpp — standalone C++ Rcpp function for log detection probabilities

**Body:**

## Summary

Add a new exported C++ Rcpp function `log_prob_detect_cpp` that computes per-location **log detection probabilities** directly from a 3-D environmental array, consolidating the logic currently spread across the R wrappers `like_neg_ltsgr_cpp()` and `log_prob_detect()` with the underlying `like_ltsg()` C++ kernel.

## Motivation

The current pipeline passes through four R layers before reaching C++:

```
loglik_math()             R/loglik_math.R
  → loglik_bio()          R/loglik_bio.R
    → log_prob_detect()   R/log_prob_detect.R
      → like_neg_ltsgr_cpp()  R/like_neg_ltsgr_cpp.R  (R wrapper)
        → like_ltsg()    src/like_ltsg.cpp             (C++ kernel)
```

Each layer adds:
- Redundant `checkmate` validation repeated at every level
- `RcppParallel::setThreadOptions` calls with `on.exit` restores at each level
- An R-side `aperm(env_dat, c(3,2,1))` + `matrix()` reshape before handing off to C++
- The formula `log(pd) - log1pexp(ctil + h)` lives in R even though `h` is fully computed in C++

For the planned `habitat_suitability()` tiled raster function (see companion issue), `log_prob_detect` will be called block-by-block in a tight loop. Each call currently pays the full R-wrapper overhead. A single C++ entry point eliminates that cost.

## Proposed Implementation

### New file: `src/log_prob_detect.cpp`

Proposed signature:

```cpp
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::export]]
Rcpp::NumericVector log_prob_detect_cpp(
    Rcpp::NumericVector env_dat,  // 3-D array (n_loc x n_time x p), R column-major
    Rcpp::NumericVector mu,       // length p
    Rcpp::NumericVector sigltil,  // length p, positive or +Inf
    Rcpp::NumericVector sigrtil,  // length p, positive or +Inf
    Rcpp::NumericMatrix o_mat,    // p x p orthogonal
    double ctil,
    double pd
);
```

Inline an extended `ColumnWorker` (RcppParallel) from `like_ltsg.cpp` that reads dims directly from `env_dat`'s `dim` attribute, eliminating the R-side `aperm + matrix()`. Aggregate over the time dimension and apply `log(pd) - log1pexp_cpp(ctil + h)` per location in C++. `Inf` entries in `sigltil`/`sigrtil` are handled automatically (`1.0 / +Inf == 0.0` under IEEE 754).

## Acceptance Criteria

- [ ] New `log_prob_detect_cpp` implemented in `src/log_prob_detect.cpp`
- [ ] `Rcpp::compileAttributes()` run; `R/RcppExports.R` and `src/RcppExports.cpp` regenerated
- [ ] `export(log_prob_detect_cpp)` added to `NAMESPACE`
- [ ] Test file `tests/testthat/test-log_prob_detect_cpp.R`:
  - Numerical parity with `log_prob_detect()` at tolerance <= 1e-12
  - Correct handling of `sigltil = Inf` and `sigrtil = Inf`
  - Return vector length equals `dim(env_dat)[1]`
- [ ] All existing tests pass unchanged

## Affected Files

| File | Change |
|---|---|
| `src/log_prob_detect.cpp` | **new** — C++ implementation |
| `src/RcppExports.cpp` | regenerated |
| `R/RcppExports.R` | regenerated |
| `NAMESPACE` | add `export(log_prob_detect_cpp)` |
| `tests/testthat/test-log_prob_detect_cpp.R` | **new** — parity + edge-case tests |

## Notes

The existing `like_neg_ltsgr_cpp()` / `like_ltsg()` stack is **kept unchanged** for backward compatibility. `log_prob_detect_cpp` is a new, higher-level C++ entry point.

---

## Issue 2

**Title:** feat: add habitat_suitability() — tiled raster evaluation using log_prob_detect_cpp

**Body:**

## Summary

Add a new exported R function `habitat_suitability()` that evaluates log detection probabilities block-by-block over a list of multi-layer `terra::SpatRaster` objects using `log_prob_detect_cpp` (see companion issue) and terra's streaming block-loop API.

## Motivation

### Current `vsp()` bottleneck

`vsp()` calls `env_data_array(env_data)` which converts the full raster to a dense R array via `terra::as.data.frame(r)`:

```r
# env_data_array.R (occ = NULL path)
as.matrix(terra::as.data.frame(r))
```

For a 10000 x 10000 raster with 39 time layers and 2 variables:
10000 x 10000 x 39 x 2 x 8 bytes = ~62 GB — unusable on any workstation.
Even the shipped 128 x 123 example raster only works because it is tiny.
The function cannot scale to real SDM use-cases.

### Terra's block-loop pattern

`terra` exposes a low-allocation streaming API:

```r
bs <- terra::blockSize(r)
terra::readStart(r)
for (b in seq_len(bs$n)) {
  chunk <- terra::readValues(r, row = bs$row[b], nrows = bs$nrows[b], mat = TRUE)
  terra::writeValues(out, result, bs$row[b], bs$nrows[b])
}
terra::readStop(r)
```

Each block fits in available RAM; calling `log_prob_detect_cpp` per block keeps memory proportional to block size, not raster size.

## Proposed Implementation

### New file: `R/habitat_suitability.R`

```r
habitat_suitability <- function(env_data, param_list, filename = "", ...) {
  # 1. Validate: env_data is a named list of SpatRaster (each with n_time layers)
  #    param_list has mu, sigltil, sigrtil, ctil, pd, o_mat
  # 2. Create output SpatRaster (1 layer, same extent/CRS as env_data[[1]])
  # 3. terra::readStart() for every raster in env_data
  # 4. terra::writeStart() for output
  # 5. bs <- terra::blockSize(env_data[[1]])
  # 6. For b in seq_len(bs$n):
  #      a. readValues() each raster -> matrix (n_cells x n_time)
  #      b. Assemble 3-D array: dim = c(n_cells, n_time, p)
  #      c. log_p <- log_prob_detect_cpp(block_arr, mu, sigltil, sigrtil,
  #                                      o_mat, ctil, pd)
  #      d. writeValues(out, exp(log_p), bs$row[b], bs$nrows[b])
  # 7. terra::writeStop(out); terra::readStop() for every input raster
  # 8. Return output SpatRaster
}
```

**NA handling:** Cells with NA in any layer/variable propagate NA in the output.
**`filename` argument:** If non-empty, write directly to disk; otherwise return an in-memory `SpatRaster`.

## Acceptance Criteria

- [ ] New `habitat_suitability()` in `R/habitat_suitability.R`
- [ ] Uses terra `readStart` / `blockSize` / `readValues` / `writeValues` / `writeStop` — no `terra::as.data.frame()` or full-raster loads
- [ ] Calls `log_prob_detect_cpp` per block
- [ ] Exported in `NAMESPACE`
- [ ] Test file `tests/testthat/test-habitat_suitability.R`:
  - Numerically identical to `exp(log_prob_detect(...))` on small in-memory rasters (same fixtures as `test-vsp.R`)
  - NA propagation: a cell with NA in any layer -> NA in output
  - Output extent, resolution, and CRS match the input
- [ ] `vsp()` kept as-is for backward compatibility
- [ ] All existing tests pass unchanged

## Affected Files

| File | Change |
|---|---|
| `R/habitat_suitability.R` | **new** — tiled raster function |
| `NAMESPACE` | add `export(habitat_suitability)` |
| `tests/testthat/test-habitat_suitability.R` | **new** |

## Dependencies

- **`terra`** already in `Suggests`; `habitat_suitability()` can keep the `requireNamespace("terra", quietly=TRUE)` guard matching `vsp()`.
- **`log_prob_detect_cpp`** (companion issue) must be merged first, or both can be implemented in the same PR.

## Depends on

Companion issue: feat: add log_prob_detect_cpp
