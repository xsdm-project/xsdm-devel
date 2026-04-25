# 01 — Incorporate xtensor for C++ array computation in an Rcpp package

## Why

R's vectorised operations on multi-dimensional arrays cross the interpreter
on every primitive (`+`, `*`, `apply`, `aperm`, ...) and allocate a fresh
SEXP for each intermediate. For per-cell computations on a 3-D array of
shape `(n_loc, ts_length, p)`, this becomes the dominant cost long before
the actual maths does.

`xtensor` is a header-only C++14 library that lets us:

- adopt a contiguous R buffer as a typed N-dimensional array **without
  copying** (`xt::adapt`);
- compose lazy, fused element-wise expressions and reductions;
- return `std::vector<double>` results that Rcpp wraps back into a
  `NumericVector` cheaply.

Reference kernels in this repo:

- `src/log_prob_detect.hpp` (xtensor headers + `log_prob_detect_tile`)
- `src/log_prob_detect.cpp` (implementation, lines 99–168)
- `src/loglik_bio.hpp` (composed kernel reusing `log_prob_detect_tile`)

## Prerequisites

### `DESCRIPTION`

```
LinkingTo:
    Rcpp,
    RcppParallel
SystemRequirements: GNU make, C++17
```

`xtensor` itself is **vendored** under `inst/include/xtensor*` (along with
`xtl`, `xsimd` if you want SIMD). Vendoring keeps the package self-contained
and avoids forcing `xtensor` as a `LinkingTo:` dependency.

### `src/Makevars`

```
CXX_STD = CXX17

PKG_CPPFLAGS = -I../inst/include \
               -I../inst/include/xtensor \
               -I../inst/include/xtl
```

The `-I../inst/include` paths resolve at `R CMD INSTALL` time because
`Makevars` is run from inside `src/`.

### `NAMESPACE`

The xtensor side is invisible to R. You only export the Rcpp-level
functions that wrap your C++ kernels (see step 5).

## Step-by-step recipe

### 1. Vendor xtensor (one-time)

Drop the `xtensor` and `xtl` source trees into `inst/include/`. Track only
the headers; the package will not link to a separate xtensor library.

```
inst/include/
├── xtensor/
│   ├── xarray.hpp
│   ├── xadapt.hpp
│   ├── xmath.hpp
│   └── xreducer.hpp
└── xtl/
    └── xtype_traits.hpp
```

### 2. Pick a canonical memory layout

Decide *once* what column- vs row-major layout your kernels expect, and
state it explicitly in the header. xsdm uses column-major
`(n_loc, ts_length, p)` so that the flat index is

```
env_dat_ptr[l + n_loc*t + n_loc*ts_length*k] = env_dat[l, t, k]
```

This matches R's `array()` default (column-major) so that
`as.numeric(env_dat)` from R is already the correct buffer.

### 3. Header-only compute kernel

Write the math against raw `const double*` pointers, not against
`Rcpp::NumericVector` — this lets you call the same kernel from another
C++ kernel later. Put it in `src/<name>.hpp`:

```cpp
#pragma once
#include <vector>
#include <xtensor/xadapt.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xreducer.hpp>

namespace mypkg {

std::vector<double> my_tile_kernel(
    const double* env_dat_ptr,
    int n_loc, int ts_length, int p,
    /* params... */
);

} // namespace mypkg
```

Use `xt::adapt` to view the raw buffer as an xtensor array **with no
copy**:

```cpp
std::vector<std::size_t> shape = {
    static_cast<std::size_t>(n_loc),
    static_cast<std::size_t>(ts_length)
};
auto sums_view = xt::adapt(
    col_sums.data(),
    static_cast<std::size_t>(n_loc * ts_length),
    xt::no_ownership(),
    shape
);
auto h_xt = xt::sum(sums_view, {1}) / (2.0 * ts_length);
```

(Pattern from `src/log_prob_detect.cpp`, lines 145–155.) `xt::no_ownership`
is critical — it tells xtensor not to free the buffer on destruction.

### 4. Implement the kernel in `src/<name>.cpp`

Hot loops should run over flat indices and use plain C++ arithmetic; reach
for xtensor only where its expression-template fusion actually wins
(reductions, broadcasts, fused element-wise pipelines). Mixing the two
freely is fine.

xsdm pattern: tight inner loops in a `RcppParallel::Worker`, then
`xt::sum(...)` over the time axis.

### 5. Rcpp-level wrapper

Add a thin `[[Rcpp::export]]` wrapper that does input validation, threads,
and `std::vector → NumericVector` conversion. Do **not** put math in here;
the kernel must remain callable from other C++ TUs without going through R.

```cpp
// [[Rcpp::export]]
Rcpp::NumericVector my_kernel_cpp(
    Rcpp::NumericVector env_dat_vec,
    Rcpp::IntegerVector env_dat_dims,
    /* params... */,
    int num_threads = 0
) {
    /* validate dims, set threads, call mypkg::my_tile_kernel,
       wrap the std::vector return value */
}
```

### 6. Regenerate Rcpp glue

```r
Rcpp::compileAttributes()
```

This rewrites `src/RcppExports.cpp` and `R/RcppExports.R`. **Never
hand-edit those files.** If you need stable export names, use the export
attribute argument:

```cpp
// [[Rcpp::export(.private_kernel)]]
```

The leading dot makes the R-side wrapper unexported.

### 7. Document the layout in roxygen

The R wrapper documentation must spell out the column-major flat-vector
contract so that callers (e.g. the `terra` block loop in recipe 3) build
the right buffer:

```
@param env_dat_vec Numeric vector. Column-major flat representation of a
  3-D array with logical dimensions c(n_loc, ts_length, p): variable k
  (1-indexed) occupies positions (k-1)*n_loc*ts_length + 1 to
  k*n_loc*ts_length, and within that block pixels (locations) vary fastest.
```

(See `R/log_prob_detect.R` and the corresponding `[[Rcpp::export]]`
wrapper in `src/log_prob_detect.cpp` for the canonical example.)

## Tests / parity checks

For every C++ kernel, write a test that compares it to the R reference at
numerical tolerance:

```r
test_that("my_kernel_cpp matches my_kernel_r within 1e-10", {
  set.seed(1); env <- array(rnorm(20*5*2), c(20, 5, 2))
  expect_equal(
    my_kernel_cpp(as.numeric(env), c(20L, 5L, 2L), ...),
    my_kernel_r(env, ...),
    tolerance = 1e-10
  )
})
```

Tolerance of `1e-10` is the right default for a kernel that does
`exp / log / log1p / log1mexp`. Tighten to `1e-12` for pure linear algebra.
This repo follows that convention in
`tests/testthat/test-log_prob_detect_cpp.R` and
`tests/testthat/test-loglik_bio_cpp.R`.

## Pitfalls

- **Forgetting `xt::no_ownership`** in `xt::adapt` causes a double-free at
  scope exit — the xtensor view tries to free a buffer it doesn't own.
- **Mixing layouts.** xtensor's default container layout is row-major;
  `xt::adapt(..., shape)` adopts row-major by default. If your R buffer is
  column-major (the R default), either pass an explicit row-major shape
  that swaps the dimensions, or use `xt::layout_type::column_major`. xsdm
  works around this by using a row-major xtensor view of an intermediate
  buffer that is itself filled in row-major order from the column-major
  R buffer (see `src/log_prob_detect.cpp` worker loop).
- **Don't keep R SEXPs in long-lived state.** A captured `NumericVector`
  inside a `std::function` that outlives the calling R frame can become
  invalid if the R-side variable is replaced. Snapshot to
  `std::vector<double>` at construction time. See recipe 4 for the worked
  example in the `ucminf` XPtr closure.
- **Input parameter typing.** Don't accept `Rcpp::NumericVector` if you
  need to pass to other Rcpp wrappers — Rcpp's auto-coercion semantics
  have surprising rules around integer/logical inputs. Coerce explicitly
  in the R wrapper (`as.integer`, `as.numeric`) before calling.

## Reusable checklist

- [ ] `inst/include/xtensor*` vendored
- [ ] `src/Makevars` adds `-I../inst/include/xtensor -I../inst/include/xtl`
- [ ] `CXX_STD = CXX17` in `Makevars`
- [ ] Header file declares the kernel with `const double*` arguments
- [ ] `.cpp` file implements the kernel (header-only also fine for small
      kernels)
- [ ] Rcpp-level wrapper does input validation, thread setup, and result
      wrapping only
- [ ] Roxygen for the R wrapper specifies the buffer layout
- [ ] Parity test against the R reference at `1e-10` or tighter
