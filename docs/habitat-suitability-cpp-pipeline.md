# C++ Detection-Probability Pipeline & `habitat_suitability()` — Implementation Report

**Date:** 2026-04-25  
**Package:** `xsdm` (v1.0.0)  
**Repository:** [alrobles/xsdm-devel](https://github.com/alrobles/xsdm-devel)

---

## 1  Motivation

The original raster-evaluation path in `vsp()` followed this data flow:

```
Disk (GeoTIFF)
  → terra::extract() / env_data_array()   # full raster → R 3-D array  ← BOTTLENECK
    → log_prob_detect()                   # R thread management
      → like_neg_ltsgr_cpp()              # R: aperm + matrix reshape + diag math
        → like_ltsg()                     # C++ RcppParallel ColumnWorker
    → exp()                               # back to R
  → terra::rast()                         # R matrix → SpatRaster
```

**Problems:**

| # | Issue | Impact |
|---|-------|--------|
| 1 | `env_data_array()` loads all pixels × all time steps × all variables into R heap at once | ~96 MB for the example 128×123×39×2 grid; GBs for production grids |
| 2 | R-level `aperm` + `matrix` reshape inside `like_neg_ltsgr_cpp()` on every call | Redundant allocation and copy per optimisation iteration |
| 3 | Four-layer R→C++ nesting (`log_prob_detect` → `like_neg_ltsgr_cpp` → `like_ltsg`) | S4/function-call overhead; repeated thread-option round-trips |
| 4 | `parallel::parLapply`-style parallelism in `vsp()` not applicable for large rasters | Serialisation cost dominates for small tiles |

---

## 2  What Was Built

### 2.1  `log_prob_detect_cpp` — C++ kernel (PR #1, PR #3)

**Files:**
- `src/log_prob_detect.hpp`
- `src/log_prob_detect.cpp`

Collapses the entire four-layer R call chain into a single `// [[Rcpp::export]]` function.

**Signature:**

```r
log_prob_detect_cpp(
  env_dat_vec,      # numeric vector — column-major flat array (n_loc × ts × p)
  env_dat_dims,     # integer vector c(n_loc, ts, p)
  mu,               # length p
  sigltil,          # length p, positive (Inf allowed)
  sigrtil,          # length p, positive (Inf allowed)
  o_mat,            # p × p orthogonal matrix
  ctil,             # scalar
  pd,               # scalar ∈ (0, 1]
  return_prob = FALSE,
  num_threads = 0L  # 0 → RcppParallel::defaultNumThreads()
)
```

**Internal pipeline (C++ only):**

```
env_dat_vec (flat, n_loc × ts × p)
  │
  ▼  LPDColumnWorker (RcppParallel::parallelFor over locations)
  │    for each location l, time t:
  │      index directly into env_dat_vec[l + n_loc*t + n_loc*ts*k]  ← zero-copy, no aperm
  │      compute dot product with o_mat_t row
  │      apply asymmetric scaling: usym = inv_sigl*dot + inv_drl*max(0,dot)
  │      accumulate usym^2 over time steps
  │    h[l] = sum_t(usym^2) / (2 * ts)
  │
  ▼  log_p[l] = log(pd) - log1pexp_cpp(ctil + h[l])
  │
  ▼  if return_prob: exp(log_p)
  │
  ▼  NumericVector of length n_loc
```

**Key implementation details:**

- **Direct index mapping** into the raw `env_dat_vec` pointer replaces the R-side `aperm` + `matrix` reshape — zero intermediate allocation
- **`log1pexp_cpp`** implements the Mächler (2012) four-region numerically stable formula, matching `R/log1pexp.R` exactly:
  - `x ≤ -37` → `exp(x)`
  - `x ≤ 18` → `log1p(exp(x))`
  - `x ≤ 33.3` → `x + exp(-x)`
  - `x > 33.3` → `x`
- **`Inf` entries** in `sigltil`/`sigrtil` handled correctly: `1/Inf = 0`
- **Thread count**: `num_threads = 0` leaves RcppParallel's global thread state unchanged; positive values set and restore on exit

### 2.2  xtensor vendoring (prerequisite)

**Files added:**
- `src/vendor/xtensor/` — xtensor 0.25.0 headers
- `src/vendor/xtl/` — xtl 0.7.7 headers
- `src/vendor/xsimd/` — xsimd 13.0.0 headers
- `src/Makevars` — `-Ivendor/xtensor -Ivendor/xtl -Ivendor/xsimd`
- `src/Makevars.win` — same flags for Windows

xtensor is used in `log_prob_detect_tile()` via `xt::adapt` for zero-copy views over the time-aggregation step.

### 2.3  `habitat_suitability()` — tiled raster driver (PR #3)

**File:** `R/habitat_suitability.R`

A memory-bounded, scalable replacement for `vsp(..., return_raster = TRUE)`.

**Signature:**

```r
habitat_suitability(
  param_list,           # list: mu, sigltil, sigrtil, o_mat, ctil, pd
  env_list,             # list of SpatRaster (one per env variable, n_time layers each)
  output      = "",     # "" = in-memory SpatRaster; non-empty = writes GeoTIFF
  overwrite   = FALSE,
  return_prob = TRUE,   # TRUE = probability; FALSE = log-probability
  threads     = RcppParallel::defaultNumThreads(),
  wopt        = list()
)
```

**Data flow — zero full-raster materialisation in R:**

```
env_list (list of SpatRaster stacks)
  │  terra::readStart() on each stack
  │
  ▼  terra::writeStart() → block schedule b (b$n blocks, b$row, b$nrows)
  │
  │  for i in 1:b$n:
  │    for k in 1:p:
  │      tile_k ← terra::readValues(env_list[[k]], row=b$row[i], nrows=b$nrows[i])
  │      # (n_tile × ts) matrix — one block of one variable
  │    assemble env_vec (n_tile × ts × p flat vector)
  │    result ← log_prob_detect_cpp(env_vec, c(n_tile, ts, p), ...)
  │    terra::writeValues(out, result, b$row[i], b$nrows[i])
  │
  ▼  terra::readStop() / terra::writeStop()
output: SpatRaster (1 band, same CRS/extent as input)
```

**Peak R memory:** proportional to one block × `p` variables — independent of full raster size.

**Validation:**
- `checkmate` guards on all parameters (param keys, `SpatRaster` types, geometry consistency across all env stacks)
- `on.exit(..., add = TRUE)` guarantees `readStop`/`writeStop` even on error

---

## 3  Test Coverage (PR #4)

### `tests/testthat/test-log_prob_detect_cpp.R` — 6 tests

| Test | What it checks |
|------|---------------|
| Cross-language parity | `log_prob_detect_cpp` ≡ `log_prob_detect()` within 1e-10 on `examples$env_array` |
| `return_prob` round-trip | `prob_result == exp(log_result)` within 1e-12 |
| Output length | returns vector of length `n_loc` |
| Log-prob range | all values ≤ 0 |
| Probability range | all values ∈ [0, 1] |
| Single-location edge case | `n_loc = 1` works and matches R reference |

### `tests/testthat/test-habitat_suitability.R` — 6 tests

| Test | What it checks |
|------|---------------|
| Parity with `vsp()` | pixel-level agreement within 1e-8 on `examples$bio01` / `examples$bio12` |
| SpatRaster dims | 1 layer, same nrow/ncol as input |
| Log-prob range | non-NA values ≤ 0 |
| Probability range | non-NA values ∈ [0, 1] |
| File output | writes valid GeoTIFF when `output` is a path |
| Geometry mismatch | error with message matching `"geometry"` |

---

## 4  Merged Pull Requests

| PR | Title | Merged |
|----|-------|--------|
| [#1](https://github.com/alrobles/xsdm-devel/pull/1) | Add `log_prob_detect_cpp`: collapse R LTSG call chain into single xtensor-accelerated C++ export | 2026-04-25 |
| [#2](https://github.com/alrobles/xsdm-devel/pull/2) | Add ISSUES.md with proposed feature specs | 2026-04-25 |
| [#3](https://github.com/alrobles/xsdm-devel/pull/3) | feat: add `log_prob_detect_cpp` + `habitat_suitability()` tiled raster evaluator | 2026-04-25 |
| [#4](https://github.com/alrobles/xsdm-devel/pull/4) | Add `log_prob_detect_cpp` C++ function, `habitat_suitability()` tiled driver, and test suites | 2026-04-25 |

---

## 5  Files Changed Summary

```
src/
  log_prob_detect.hpp          # new — xsdm namespace, log1pexp_cpp, tile function declaration
  log_prob_detect.cpp          # new — LPDColumnWorker, tile implementation, Rcpp::export
  RcppExports.cpp              # regenerated — registers log_prob_detect_cpp
  Makevars                     # updated — -Ivendor/xtensor -Ivendor/xtl -Ivendor/xsimd
  Makevars.win                 # updated — same flags for Windows
  vendor/
    xtensor/                   # new — xtensor 0.25.0 headers
    xtl/                       # new — xtl 0.7.7 headers
    xsimd/                     # new — xsimd 13.0.0 headers

R/
  habitat_suitability.R        # new — tiled terra block-loop driver
  RcppExports.R                # regenerated — exposes log_prob_detect_cpp

NAMESPACE                      # updated — export(log_prob_detect_cpp)
                               #            export(habitat_suitability)

tests/testthat/
  test-log_prob_detect_cpp.R   # new — 6 cross-language / range tests
  test-habitat_suitability.R   # new — 6 parity / output / validation tests
```

---

## 6  Usage Examples

### Evaluate log detection probability for extracted locations

```r
library(xsdm)

env  <- examples$env_array    # 2000 × 39 × 2 array
pl   <- examples$par_list

# C++ path (replaces log_prob_detect())
log_p <- log_prob_detect_cpp(
  env_dat_vec  = as.vector(env),
  env_dat_dims = dim(env),
  mu      = pl$mu,
  sigltil = pl$sigltil,
  sigrtil = pl$sigrtil,
  o_mat   = pl$o_mat,
  ctil    = pl$ctil,
  pd      = pl$pd
)
```

### Generate a habitat suitability map over large rasters

```r
library(xsdm)
library(terra)

bio1  <- unwrap(examples$bio01) / 100   # unpack + rescale
bio12 <- unwrap(examples$bio12) / 100

# Tiled — never loads full raster into R memory
suit <- habitat_suitability(
  param_list  = examples$par_list,
  env_list    = list(bio1, bio12),
  return_prob = TRUE,
  output      = "habitat_suitability.tif",
  overwrite   = TRUE
)

plot(suit)
```

---

## 7  Architecture Comparison

| Aspect | Before (vsp) | After (habitat_suitability) |
|--------|-------------|----------------------------|
| R memory (128×123 grid) | ~96 MB full array | ~1 block (~1–4 MB) |
| R memory (10k×10k grid) | ~62 GB | ~1 block (~100 KB–1 MB) |
| R→C++ call chain depth | 4 layers | 1 layer |
| Parallelism | none / PSOCK serialisation | RcppParallel `parallelFor` over locations |
| `aperm` + `matrix` reshape | R-level, every call | C++ index mapping, zero allocation |
| `Inf` sigma values | handled in R | handled in C++ (1/Inf = 0) |
| Thread-state side effects | global `setThreadOptions` reset | local, restored on exit |
| File output | `terra::rast()` from data.frame | direct `terra::writeValues` streaming |

---

## 8  Next Steps

- [ ] Benchmark `habitat_suitability()` vs `vsp()` on the full 128×123×39 example grid with `bench::mark()`
- [ ] Validate on a production-scale grid (0.5° global, ~259,200 cells)
- [ ] Deprecation notice on `vsp()` pointing users to `habitat_suitability()`
- [ ] Consider moving `log_prob_detect()` R chain to an internal-only function once cross-language tests are stable in CI
- [ ] GDAL direct I/O path (Disk → GDAL → xtensor → GDAL → Disk) to eliminate the remaining terra R wrapper, following the `rxbioclim` `GdalReader`/`GdalWriter` pattern
