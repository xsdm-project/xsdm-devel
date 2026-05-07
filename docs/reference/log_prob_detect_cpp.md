# Compute log detection probabilities from a flat environmental data vector

C++ implementation of
[`log_prob_detect()`](https://xsdm-project.github.io/xsdm-devel/reference/log_prob_detect.md)
that accepts environmental data as a flat numeric vector with explicit
dimension metadata. This signature is designed for block-by-block raster
evaluation where each block is passed as a contiguous vector rather than
a 3-D R array.

## Usage

``` r
log_prob_detect_cpp(
  env_dat_vec,
  env_dat_dims,
  mu,
  sigltil,
  sigrtil,
  o_mat,
  ctil,
  pd,
  return_prob = FALSE,
  num_threads = 0L
)
```

## Arguments

- env_dat_vec:

  Numeric vector. Column-major flat representation of a 3-D array with
  logical dimensions `c(n_loc, ts_length, p)`: variable `k` (1-indexed)
  occupies positions `(k-1)*n_loc*ts_length + 1` to `k*n_loc*ts_length`,
  and within that block pixels (locations) vary fastest.

- env_dat_dims:

  Integer vector of length 3: `c(n_loc, ts_length, p)`.

- mu:

  Numeric vector of length `p`. Optimal environmental values.

- sigltil:

  Numeric vector of length `p`. Positive; `Inf` entries are allowed
  (treated as zero inverse-scale).

- sigrtil:

  Numeric vector of length `p`. Positive; `Inf` entries are allowed.

- o_mat:

  Numeric matrix, `p x p` orthogonal.

- ctil:

  Scalar. Center of the detection-link function.

- pd:

  Scalar in `(0, 1]`. Maximum probability of detection.

- return_prob:

  Logical. If `TRUE`, return probabilities; if `FALSE` (default) return
  log-probabilities.

- num_threads:

  Integer. Number of parallel threads. `0` (default) uses
  [`RcppParallel::defaultNumThreads()`](https://rdrr.io/pkg/RcppParallel/man/setThreadOptions.html).

## Value

Numeric vector of length `n_loc`.

## Details

Collapses the R call chain `like_neg_ltsgr_cpp() -> like_ltsg()` into a
single xtensor-accelerated C++ function.
