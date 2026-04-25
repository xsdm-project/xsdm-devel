# 02 — Develop `log_prob_detect_cpp` (collapse an R call chain into a single C++ kernel)

## Why

The R reference implementation of `log_prob_detect` was the composition of
several specialised functions:

```
log_prob_detect(env, params)
  └─ like_neg_ltsgr_cpp(env, params)
       ├─ matrix permutation (aperm)
       ├─ centering (env - mu)
       ├─ rotation (t(o_mat) %*% .)
       ├─ asymmetric scaling (sigltil / sigrtil split)
       ├─ squared sum across variables (like_ltsg)
       └─ time aggregation
  └─ log1pexp + log(pd) and exp() if return_prob
```

Each step allocated a fresh R array of shape `(n_loc, ts_length, p)`. For a
`60_000 × 30 × 2` block this is hundreds of MB of temporary R objects per
evaluation, plus interpreter overhead on every step.

`log_prob_detect_cpp` collapses the entire chain into a **single
C++ kernel** that runs in one pass, parallelised across `(location, time)`
columns and time-aggregated with xtensor.

Reference files in this repo:

- `src/log_prob_detect.cpp` (the Rcpp-level wrapper + kernel implementation)
- `src/log_prob_detect.hpp` (kernel signature + numerical helpers)
- `R/log_prob_detect.R` (R reference, unchanged, kept as parity oracle)
- `tests/testthat/test-log_prob_detect_cpp.R` (parity tests)

## Prerequisites

- xtensor wired in per [recipe 1](01-xtensor-cpp-arrays.md).
- `LinkingTo: RcppParallel` in `DESCRIPTION`.
- An R reference implementation that you can keep around as the
  ground-truth oracle for parity tests. **Do not delete the R version** —
  it is your most valuable debugging tool when the C++ path drifts.

## Step-by-step recipe

### 1. Identify the call chain to collapse

List every R function that runs between user input and the final scalar /
vector. For each, note: input shape, output shape, allocations.

For xsdm: the chain produced an asymmetric Mahalanobis-like distance per
`(location, time)`, summed over variables; the time aggregation came last.

### 2. Define the boundary of the C++ kernel

The kernel should take **only what survives compaction**:

- one flat numeric buffer for the environmental data
- the integer dim triple `(n_loc, ts_length, p)`
- raw pointers / scalars for parameters
- output shape (here: `n_loc`)

It should not take an `Rcpp::NumericVector` that requires the R interpreter
to be alive. (See `src/log_prob_detect.hpp` lines 52–64 for the canonical
signature.)

### 3. Pure-C++ helpers

Numerically-stable elementary functions get inline implementations in the
header. xsdm provides:

```cpp
inline double log1pexp(double x) {
    if (x <= -37.0) return std::exp(x);
    if (x <= 18.0)  return std::log1p(std::exp(x));
    if (x <= 33.3)  return x + std::exp(-x);
    return x;
}
```

Match the R reference's defaults exactly (here: `c0=-37, c1=18, c2=33.3`)
so that parity tests pass at 1e-10.

### 4. Parallel inner loop with `RcppParallel::Worker`

Define a `Worker` struct that operates on a single `(location, time)`
column at a time. The body of `operator()(begin, end)` is a flat loop
`for (j = begin; j < end; ++j)`, decoded back to `(l, t) = (j / ts, j % ts)`.

Pattern from `src/log_prob_detect.cpp` lines 30–94:

```cpp
struct LPDColumnWorker : public RcppParallel::Worker {
    const double* env_dat_ptr;
    /* ... pointers to inv_sigl, inv_drl, mu, o_mat ... */
    double*       output;

    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t j = begin; j < end; ++j) {
            int l = j / ts_length, t = j % ts_length;
            /* asymmetric Mahalanobis squared sum: write to output[j] */
        }
    }
};
```

Then dispatch with `RcppParallel::parallelFor(0, total, worker)`.

### 5. Time aggregation with xtensor

After the parallel pass produces a `n_loc × ts_length` flat buffer, fold
the time axis with xtensor. This is where xtensor pays off — fused
expressions, no extra allocation:

```cpp
auto sums_view = xt::adapt(col_sums.data(), total, xt::no_ownership(), shape);
auto h_xt      = xt::sum(sums_view, {1}) / (2.0 * ts_length);
```

(See `src/log_prob_detect.cpp` lines 142–155.)

### 6. Final scalar transform

Compute the final per-location output in a plain loop, applying
`log1pexp`, the `log(pd)` correction, and the `exp` if the user asked for
probabilities:

```cpp
const double log_pd = std::log(pd);
for (int l = 0; l < n_loc; ++l) {
    double val = log_pd - log1pexp(ctil + h_xt(l));
    result[l] = return_prob ? std::exp(val) : val;
}
```

### 7. Rcpp wrapper layer

The `[[Rcpp::export]]` function does:

1. Input validation (`Rcpp::stop` with clear messages on mismatched dims).
2. Optional thread count override via the `RcppParallel` namespace
   (`setThreadOptions`).
3. Delegate to `xsdm::log_prob_detect_tile`.
4. Wrap the `std::vector<double>` return value in a `NumericVector`.

(See `src/log_prob_detect.cpp` lines 207–268.)

### 8. Inverse-scaling vectors

Pre-compute `1/sigltil` and `1/sigrtil - 1/sigltil` once before the
parallel loop. **Handle `Inf` via IEEE 754** (i.e. `1.0 / Inf == 0.0`)
rather than special-casing — this matches the R reference behaviour:

```cpp
double isl = std::isinf(sigltil_ptr[i]) ? 0.0 : 1.0 / sigltil_ptr[i];
double isr = std::isinf(sigrtil_ptr[i]) ? 0.0 : 1.0 / sigrtil_ptr[i];
inv_sigl[i] = isl;
inv_drl[i]  = isr - isl;
```

(See `src/log_prob_detect.cpp` lines 119–126.)

## Tests / parity checks

`tests/testthat/test-log_prob_detect_cpp.R` runs the same `(env, params)`
through the R reference and the new C++ kernel and asserts numerical
equality:

```r
expect_equal(
  log_prob_detect_cpp(as.numeric(env), c(n, ts, p), mu, sigltil, sigrtil,
                      o_mat, ctil, pd, return_prob = FALSE, num_threads = 1L),
  log_prob_detect(env, list(mu = mu, sigltil = sigltil, ...)),
  tolerance = 1e-10
)
```

For multivariate cases, run with `p = 1, 2, 3, 5` to exercise the
`o_mat` rotation in non-trivial dimensions. Also add a test for `Inf`
entries in `sigltil` / `sigrtil` since those exercise the IEEE 754 path.

## Pitfalls

- **`aperm` confusion.** R's `aperm(env_dat, c(3, 2, 1))` is a logical
  permutation of dims — but the *underlying buffer* stays the same. When
  porting to C++, you don't actually need to permute; just rewrite the
  index expression. xsdm exploits this in
  `LPDColumnWorker::operator()` (lines 70–84): the original R code did an
  `aperm` and then accessed `mat[k, j]`, but the C++ kernel reads
  directly from `env_dat_ptr[l + n_loc*t + n_loc*ts_length*k]` — no
  permutation needed.
- **Column-major vs row-major surprises with `o_mat`.** R matrices are
  column-major. `t(o_mat)[i, k] == o_mat[k, i]` flat-indexes to
  `o_mat_ptr[k + p*i]` (column-major). Get this wrong and your parity
  test breaks at the third decimal — close enough that `tolerance = 1e-3`
  passes but actual results are wrong.
- **`Inf` handling.** Don't write `if (std::isinf(...)) continue;` in the
  inner loop; rely on IEEE 754 (`1.0 / Inf == 0.0`) so the asymmetric
  scaling path produces zero contribution naturally. This both matches R
  semantics and avoids branching in the hot loop.
- **`std::vector<double>` return + `Rcpp::wrap`.** This copies. For very
  large outputs use `Rcpp::NumericVector` allocated by the wrapper and
  filled by the kernel via raw pointer. For `n_loc ≈ 1e6` the copy is
  ~8 MB and rarely matters; profile before optimising.
- **Thread set-up across calls.** `RcppParallel::setThreadOptions` is
  process-global. If you change it in one wrapper, save and restore it
  around your call so the user's setting isn't clobbered. xsdm does this
  in `src/loglik_math.cpp` lines 76–95.

## Reusable checklist

- [ ] R reference function preserved unchanged
- [ ] C++ header declares the tile kernel with `const double*` args
- [ ] Worker struct decoding `j` to `(l, t)` flat indices
- [ ] xtensor used **only** for the reduction step (or wherever fusion
      genuinely helps)
- [ ] Numerically-stable `log1pexp` / `log1mexp` inline helpers
- [ ] Rcpp wrapper validates input dims and thread count
- [ ] Parity test at `1e-10` against the R reference
- [ ] Stress test with `Inf` entries in scale parameters
