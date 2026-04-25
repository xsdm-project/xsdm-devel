# 03 — `terra` block-wise I/O and raw pointers into `log_prob_detect_cpp`

## Why

A continental-scale environmental raster stack is far larger than RAM:

```
terra::SpatRaster
├─ bio01:  9000 x 9000 pixels x 30 layers x 8 bytes = 19.5 GB
├─ bio12:  9000 x 9000 pixels x 30 layers x 8 bytes = 19.5 GB
└─ ...
```

We cannot materialise the whole `(n_loc, ts_length, p)` array in R memory
to call `log_prob_detect_cpp` once. Instead we use `terra`'s streaming
block-loop API — read a horizontal strip, evaluate the kernel on it, write
the strip out, drop it, and move on. Memory stays bounded by the largest
block, regardless of total raster size.

Reference implementation: `R/habitat_suitability.R` (lines 79–227).

## Prerequisites

- `terra >= 1.7-46` (block-loop API stable since 2023).
- The C++ kernel from [recipe 2](02-log_prob_detect_cpp.md) accepts a
  flat column-major buffer and explicit dim triple.
- `Imports: terra, checkmate` in `DESCRIPTION`.

## Step-by-step recipe

### 1. Validate inputs and geometry

Before opening any file handles, check that all rasters share the same
extent / resolution / CRS:

```r
ref <- env_list[[1L]]
for (k in seq.int(2L, length(env_list))) {
  if (!terra::compareGeom(ref, env_list[[k]], stopOnError = FALSE)) {
    stop("All rasters in env_list must have the same geometry.")
  }
}
```

Pull `p`, `ts`, `nc` (cols per row) from the reference raster:

```r
p  <- length(env_list)
ts <- terra::nlyr(ref)
nc <- terra::ncol(ref)
```

(See `R/habitat_suitability.R` lines 114–129.)

### 2. Open every input raster for streaming reads

```r
for (k in seq_len(p)) terra::readStart(env_list[[k]])
on.exit({
  for (k in seq_len(p)) {
    tryCatch(terra::readStop(env_list[[k]]), error = function(e) NULL)
  }
}, add = TRUE)
```

`readStart` opens the file (or marks the in-memory object) for cursor-style
reads. Pair it with an `on.exit` that calls `readStop` so a thrown error
during the loop still releases the file handle.

### 3. Open the output for streaming writes

```r
out <- terra::rast(ref, nlyr = 1L)
b <- terra::writeStart(out, filename = output, overwrite = overwrite, wopt = wopt)
on.exit(tryCatch(terra::writeStop(out), error = function(e) NULL), add = TRUE)
```

`writeStart` returns a list with the block schedule chosen by terra's
memory manager:

```
b$n      : number of blocks
b$row    : starting row of block i (1-indexed)
b$nrows  : number of rows in block i
```

The schedule is automatic — terra picks `b$nrows[i]` so each block fits
within `terra::terraOptions("memmin")`.

### 4. The block loop

```r
for (i in seq_len(b$n)) {
  n_tile <- b$nrows[i] * nc
  env_vec <- numeric(n_tile * ts * p)
  valid   <- rep(TRUE, n_tile)
  ...
}
```

`n_tile` = pixels in this block. We allocate ONE flat numeric vector of
length `n_tile * ts * p` for the entire block and fill it in column-major
order to match the kernel's contract.

### 5. Read each variable's tile and pack into the flat buffer

```r
for (k in seq_len(p)) {
  tile_k <- terra::readValues(
    env_list[[k]],
    row = b$row[i], nrows = b$nrows[i],
    col = 1L,      ncols = nc,
    mat = TRUE                        # returns n_tile x ts matrix
  )
  offset <- (k - 1L) * n_tile * ts
  env_vec[offset + seq_len(n_tile * ts)] <- as.vector(tile_k)
  valid <- valid & !apply(tile_k, 1L, anyNA)
}
```

(See `R/habitat_suitability.R` lines 163–176.)

`mat = TRUE` makes `readValues` return a matrix with one row per pixel and
one column per time step. R stores it column-major, which is exactly what
the kernel expects: variable `k` occupies the slab
`[k * n_tile * ts, (k+1) * n_tile * ts)` of `env_vec`, and within that slab
pixels vary fastest within a time step.

### 6. NA masking via compaction

Cells that are NA in any variable / any time step must not be passed to
the kernel — they would propagate `NaN` and pollute neighbouring numerical
results. **Compact the buffer before calling the kernel** and re-insert
`NA` afterward:

```r
n_valid <- sum(valid)
block_result <- rep(NA_real_, n_tile)

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
      mu = ..., sigltil = ..., sigrtil = ..., o_mat = ...,
      ctil = ..., pd = ...,
      return_prob = return_prob,
      num_threads = as.integer(threads)
    )
    block_result[valid] <- result_valid
  } else {
    block_result <- log_prob_detect_cpp(env_dat_vec = env_vec, ...)
  }
}
```

(See `R/habitat_suitability.R` lines 178–221.)

### 7. Write the block

```r
terra::writeValues(out, block_result, b$row[i], b$nrows[i])
```

`writeValues` accepts a length-`n_tile` numeric vector for a single layer.
The cell order matches `readValues` (row-major within the block).

### 8. Return appropriately

If the user supplied a filename, return invisibly so they don't get the
whole `SpatRaster` printed at the REPL:

```r
if (nzchar(output)) invisible(out) else out
```

## How raw pointers reach the C++ kernel

The R-to-C++ pointer hand-off happens entirely inside the
`[[Rcpp::export]]` wrapper of `log_prob_detect_cpp`:

```
R                              C++
─                              ───
env_vec (numeric)         ──►  Rcpp::NumericVector env_dat_vec
                               REAL(env_dat_vec)              ──►  const double* env_dat_ptr
env_dat_dims (integer)    ──►  Rcpp::IntegerVector env_dat_dims
                               env_dat_dims[0..2]             ──►  int n_loc, ts_length, p
```

`REAL(env_dat_vec)` is `&env_dat_vec[0]` — a raw `double*` into R's heap.
The kernel reads this pointer with no copy, runs the parallel loop and
xtensor reduction, and returns a `std::vector<double>` that Rcpp wraps
back into a `NumericVector`. (See `src/log_prob_detect.cpp` lines 254–262.)

For the optimization path, the wrapping is even tighter — see
[recipe 4](04-ucminfcpp-xptr.md): the env_dat buffer is **snapshotted
once** into an owned `std::vector<double>` so the optimizer's hot loop
doesn't even cross the Rcpp boundary.

## Tests / parity checks

`tests/testthat/test-habitat_suitability.R` builds a tiny in-memory
`SpatRaster` stack, runs `habitat_suitability(...)`, and asserts:

1. Output has correct geometry, layer count, and layer name.
2. Per-pixel results equal `log_prob_detect_cpp` called on the unblocked
   array (parity at `tolerance = 1e-12`, since this is just a permutation
   check).
3. NA cells in input become NA in output.
4. `return_prob = TRUE` produces values in `[0, 1]`.

## Pitfalls

- **Forgetting `on.exit(readStop)`**. If an error fires mid-loop, the file
  handle stays open until R exits. On Windows this leaves the file locked.
- **Reading variables in the wrong order.** The kernel's column-major
  contract puts variable `k` at offset `(k - 1L) * n_tile * ts`. Flipping
  variable order between read and pack is a silent correctness bug — the
  parity test in step 2 of recipe 2 will catch it, but only if you
  actually wrote one.
- **Calling `readValues(..., mat = FALSE)`** returns a flattened vector
  with a different layout than `mat = TRUE`. Always use `mat = TRUE` and
  `as.vector(tile_k)` to get the column-major layout the kernel expects.
- **NA propagation.** Even *one* NA in a single time step of a single
  variable for a single pixel will produce `NaN` in the kernel output for
  that pixel. Always do the `valid` mask + compaction.
- **Block schedule depends on memory settings.** `b$n` and `b$nrows[i]`
  vary with `terra::terraOptions("memmin")` and the OS's free RAM. Tests
  that hard-code `b$n` will be flaky on CI; assert per-pixel output
  instead of per-block bookkeeping.
- **`writeValues` cell order.** It expects values in *row-major within the
  block* — same order as `readValues` returns. The compaction pattern in
  step 6 preserves this because `valid` is a length-`n_tile` logical mask
  in the original cell order.

## Reusable checklist

- [ ] Geometry validation across `env_list` before any file ops
- [ ] `readStart` for every input + `on.exit(readStop)`
- [ ] `writeStart` for output + `on.exit(writeStop)`
- [ ] Block loop using the schedule from `writeStart`
- [ ] Flat column-major buffer assembled from `readValues(mat = TRUE)`
- [ ] NA mask built per-pixel; compaction before the C++ kernel
- [ ] Re-insertion of NA in the output buffer before `writeValues`
- [ ] Return invisibly when writing to disk
