# 04 — Use `ucminfcpp::ucminf_xptr` to remove R-callback overhead from optimization

## Why

The default optimizer interface accepts an R function as the objective:

```r
ucminfcpp::ucminf(par, fn, gr = NULL, control = list(...))
```

Every `(f, g)` evaluation pays an R-interpreter round-trip:

```
ucminfcpp (C++) → R interpreter → user fn (R) → R interpreter → ucminfcpp
```

For `n_free` parameters and central differences, that is `1 + 2 * n_free`
round-trips per gradient step. With `n_free = 9` and 200 BFGS iterations,
that's ~3800 R round-trips just for one optimization run — and each
round-trip allocates a fresh SEXP, blocks GC, etc.

`ucminfcpp::ucminf_xptr(par, xptr, control)` accepts an
`Rcpp::XPtr<ucminf::ObjFun>` instead. The optimizer calls the wrapped
C++ `std::function` directly. Properly built, this eliminates 100% of the
interpreter overhead inside the optimizer's inner loop.

The tricky parts are:

1. capturing the data buffers (env_dat, occ) **once** at construction so
   the closure never touches the R heap;
2. mapping the optimizer's flat `x` vector back to the named-parameter
   schema your kernel expects;
3. handling masks (parameters held fixed during optimization);
4. isolating the closure from exceptions so a bad evaluation can't crash
   R.

Reference implementation:

- `src/loglik_math_xptr_cpp.cpp` (the new pure-C++ XPtr builder)
- `src/loglik_math_xptr.cpp` (legacy R-callback version, kept for
  reference)
- `R/optimize_likelihood.R`, `optimize_loglik_math_()` (caller)
- `tests/testthat/test-loglik_math_xptr_parity.R`
- `inst/benchmarks/loglik_cpp_vs_r.R` (10–145× kernel speedup)

## Prerequisites

- A C++ implementation of the objective callable from a `const double*`
  parameter pointer (i.e. a kernel from [recipe 2](02-log_prob_detect_cpp.md)).
- `LinkingTo: Rcpp, ucminfcpp` in `DESCRIPTION`.
- Access to `inst/include/ucminf_core.hpp` from the `ucminfcpp` package
  (header-only API).

## The `ucminf::ObjFun` contract

```cpp
using ObjFun = std::function<void(
    const std::vector<double>& x,   // current point, length n
    std::vector<double>&       g,   // gradient out, length n (pre-sized)
    double&                    f    // objective value out
)>;
```

The optimizer pre-sizes `g` to `x.size()`. The callback must write `f`
**and** all `n` gradient components.

## Step-by-step recipe

### 1. Decide where to snapshot the data

Anything passed to the closure that originates as an R object must be
copied to an owned C++ container at construction time. Storing
`Rcpp::NumericVector` inside a `std::function` that survives across R
top-level calls is unsafe.

xsdm pattern: a `XptrClosureState` struct captured by `shared_ptr`:

```cpp
struct XptrClosureState {
    std::vector<double> env_dat;          // snapshot
    std::vector<int>    occ;              // snapshot
    int n_loc, ts, p;

    std::vector<std::string> free_names;  // closure parameter schema
    bool                     has_mask = false;
    std::vector<std::string> mask_names;
    std::vector<double>      mask_values;

    std::vector<double> full;             // canonical param vector,
                                          // mask values pre-merged
    std::vector<int>    free_slot_idx;    // x[i] writes to full[free_slot_idx[i]]
    bool                ready = false;    // set on first call

    BioParams           bp;               // scratch (resized once)

    double gradstep_rel = 1e-6, gradstep_abs = 1e-8;
    bool   use_central  = true;
};
```

(See `src/loglik_math_xptr_cpp.cpp` lines 57–83.)

### 2. Build the XPtr-returning Rcpp wrapper

```cpp
// [[Rcpp::export]]
SEXP make_loglik_math_xptr_cpp(
    Rcpp::NumericVector env_dat,
    Rcpp::IntegerVector occ,
    Rcpp::Nullable<Rcpp::NumericVector> mask,
    Rcpp::CharacterVector free_names,
    int num_threads = 1,
    std::string grad = "central",
    Rcpp::NumericVector gradstep = Rcpp::NumericVector::create(1e-6, 1e-8))
{
    /* validate args, fill state, build closure, return XPtr */
}
```

(See `src/loglik_math_xptr_cpp.cpp` lines 116+.)

### 3. Validate inputs at construction (clear errors > deferred crashes)

Defensive checks **at the build site** translate to R errors with stack
traces, instead of segfaults inside the optimizer:

```cpp
if (gradstep.size() != 2)       Rcpp::stop("`gradstep` must have length 2.");
if (!(gradstep[0] > 0.0) || !(gradstep[1] > 0.0))
    Rcpp::stop("`gradstep` entries must be strictly positive.");
if (free_names.size() == 0)     Rcpp::stop("`free_names` must be non-empty.");

// Critical: unnamed mask must be rejected explicitly.
if (mask.isNotNull()) {
    Rcpp::NumericVector mv(mask);
    if (mv.size() > 0) {
        SEXP names_sexp = Rf_getAttrib(mv, R_NamesSymbol);
        if (names_sexp == R_NilValue) {
            Rcpp::stop("`mask` must be a named numeric vector.");
        }
        Rcpp::CharacterVector mn(names_sexp);
        if (mn.size() != mv.size()) {
            Rcpp::stop("`mask` names length does not match values length.");
        }
        /* copy mn[i] / mv[i] into state */
    }
}
```

Why the unnamed-mask check is non-optional: `Rcpp::NumericVector::names()`
on an unnamed vector returns `R_NilValue`, which becomes a zero-length
`CharacterVector`. Iterating `mn[i]` past `mn.size()` reads out of bounds.
On Windows this can produce **a fatal R session abort** rather than an R
error.

### 4. Snapshot data into owned containers

```cpp
auto state = std::make_shared<XptrClosureState>();
state->env_dat.assign(env_dat.begin(), env_dat.end());
state->occ.assign(occ.begin(), occ.end());
state->n_loc = ...; state->ts = ...; state->p = ...;
state->bp.resize(p);
```

After this, the closure does not need any R object.

### 5. Defer canonical-vector setup to the first call (lazy init)

If you build the canonical-vector mapping eagerly inside the wrapper, you
fail fast on bogus names — but you also break tests that mock the
optimizer and never invoke the closure (e.g. `test-optimize_helpers.R`).

Lazy init resolves this: store free_names / mask_names as plain strings,
and on the first ObjFun call:

```cpp
const auto canon = canonical_names(state->p);  // e.g. "mu1"..."o_par1"
const int N = static_cast<int>(canon.size());
state->full.assign(N, std::numeric_limits<double>::quiet_NaN());

std::vector<bool> slot_filled(N, false);
if (state->has_mask) {
    for (size_t i = 0; i < state->mask_names.size(); ++i) {
        const int k = find_canonical_index(canon, state->mask_names[i]);
        if (k < 0) throw std::runtime_error(
            "ObjFun: unknown canonical name in mask: " + state->mask_names[i]);
        state->full[k] = state->mask_values[i];
        slot_filled[k] = true;
    }
}
state->free_slot_idx.clear();
state->free_slot_idx.reserve(state->free_names.size());
for (const auto& nm : state->free_names) {
    const int k = find_canonical_index(canon, nm);
    if (k < 0) throw std::runtime_error(
        "ObjFun: unknown canonical name in free_names: " + nm);
    if (slot_filled[k]) throw std::runtime_error(
        "ObjFun: free_names and mask overlap on: " + nm);
    state->free_slot_idx.push_back(k);
    slot_filled[k] = true;
    state->full[k] = 0.0;
}
for (int i = 0; i < N; ++i)
    if (!slot_filled[i]) throw std::runtime_error("ObjFun: slot uncovered");
state->ready = true;
```

(See `src/loglik_math_xptr_cpp.cpp` lines 234–272.)

`canonical_names(p)` returns the kernel's parameter schema in fixed order;
keep it pure C++ (no SEXP) so the closure has zero R dependency after
construction.

### 6. Splice and evaluate on each call

Once the canonical vector is set up, the per-call work is just an index
splice and a kernel call:

```cpp
const int n_free = static_cast<int>(x.size());
for (int i = 0; i < n_free; ++i)
    state->full[state->free_slot_idx[i]] = x[i];

const double ll = xsdm::loglik_math_eval(
    state->full.data(), state->p,
    state->env_dat.data(), state->occ.data(),
    state->n_loc, state->ts, state->bp);
f = -ll;       // ucminfcpp minimises
```

### 7. Finite-difference gradient

For each free slot, perturb `state->full[k]`, re-evaluate, and restore:

```cpp
for (int i = 0; i < n_free; ++i) {
    const int    k  = state->free_slot_idx[i];
    const double xi = state->full[k];
    const double dx = std::abs(xi) * state->gradstep_rel + state->gradstep_abs;

    state->full[k] = xi + dx;
    const double f_plus = -xsdm::loglik_math_eval(
        state->full.data(), state->p, state->env_dat.data(),
        state->occ.data(), state->n_loc, state->ts, state->bp);

    if (state->use_central) {
        state->full[k] = xi - dx;
        const double f_minus = -xsdm::loglik_math_eval(...);
        g[i] = (f_plus - f_minus) / (2.0 * dx);
    } else {
        g[i] = (f_plus - f) / dx;
    }
    state->full[k] = xi;       // restore
}
```

Because the canonical vector is mutated in place, no per-step allocation
happens. The `state->bp` BioParams is also reused. (See
`src/loglik_math_xptr_cpp.cpp` lines 314–342.)

### 8. **Catch all exceptions at the closure boundary**

The optimizer's template stack is not exception-safe: a thrown
`std::runtime_error` unwinding through `ucminf::minimize_direct<F>` leaves
the trust-region / line-search workspace partially updated. While Rcpp's
outer `BEGIN_RCPP` / `END_RCPP` wrapper does ultimately catch
`std::exception`, you should catch sooner — at the closure boundary —
and surface failure as a sentinel value:

```cpp
try {
    /* lazy init + evaluate + gradient */
} catch (const std::exception& e) {
    f = std::numeric_limits<double>::infinity();
    std::fill(g.begin(), g.end(), 0.0);
    REprintf("xsdm XPtr ObjFun error: %s\n", e.what());
} catch (...) {
    f = std::numeric_limits<double>::infinity();
    std::fill(g.begin(), g.end(), 0.0);
    REprintf("xsdm XPtr ObjFun error: unknown C++ exception\n");
}
```

A finite gradient set to zero plus `f = +Inf` is a valid "abandon this
direction" signal that the optimizer can act on without corrupting state.

### 9. Construct and return the XPtr

```cpp
auto fn = std::make_unique<ucminf::ObjFun>(
    [state](const std::vector<double>& x,
            std::vector<double>&       g,
            double&                    f) { /* body from steps 5–8 */ });

return Rcpp::XPtr<ucminf::ObjFun>(fn.release(), true);  // gc-finalize
```

The lambda captures `state` by value (i.e. by `shared_ptr` copy). When R
GCs the XPtr, the unique_ptr deleter runs, the std::function dies, the
lambda destructs, the shared_ptr's refcount goes to zero, and the state
buffer is freed.

### 10. Wire it into the R caller

```r
loglik_xptr <- make_loglik_math_xptr_cpp(
  env_dat     = env_dat,
  occ         = occ_i,
  mask        = mask,
  free_names  = names(param_vector),
  num_threads = num_threads,
  grad        = grad_ctrl$grad,
  gradstep    = grad_ctrl$gradstep
)

res <- ucminfcpp::ucminf_xptr(
  par     = param_vector,
  xptr    = loglik_xptr,
  control = ctrl,
  hessian = FALSE
)
```

(See `R/optimize_likelihood.R` lines 295–321.)

## Tests / parity checks

`tests/testthat/test-loglik_math_xptr_parity.R` covers:

1. **Optimizer parity** — running ucminf for `maxeval = 50` produces an
   objective at convergence that matches `loglik_math(par_final, negative = TRUE)`
   to `tolerance = 1e-8`.
2. **Single-step parity** — `maxeval = 1` matches `loglik_math_cpp` at
   `par_final` to the same tolerance.
3. **Mask parity** — same checks with a non-empty mask.
4. **Input validation** — unnamed mask and non-positive `gradstep` are
   rejected with clear R errors before any closure call.
5. **Crash resilience** — an XPtr built with incomplete free_names
   coverage (would throw inside the closure) surfaces `f = +Inf` and lets
   `ucminf_xptr` return gracefully instead of aborting R.

## Benchmarks

`inst/benchmarks/loglik_cpp_vs_r.R` runs the same kernel under both
implementations:

| n     | Tt | p | reps | R time | C++ time | Speedup |
|-------|----|---|------|--------|----------|---------|
| 10    | 5  | 2 | 200  | 0.724s | 0.005s   | 145×    |
| 50    | 10 | 2 | 100  | 0.090s | 0.004s   | 23×     |
| 200   | 20 | 2 | 50   | 0.053s | 0.003s   | 18×     |
| 1000  | 20 | 2 | 10   | 0.019s | 0.002s   | 10×     |

The speedup is largest for small problems where R-callback overhead
dominates real compute time.

## Pitfalls

- **Storing `Rcpp::NumericVector` inside the closure** is the single most
  common bug. The captured SEXP can become invalidated when the calling R
  frame goes out of scope or when the user reassigns the variable. Always
  snapshot to `std::vector<double>`.
- **Eager validation breaks mocked-optimizer tests.** Lazy init (step 5)
  is the cleanest fix.
- **Forgetting `Rcpp::Nullable<>`** for the optional mask. Plain
  `Rcpp::NumericVector` cannot be NULL; `Rcpp::Nullable<Rcpp::NumericVector>`
  can.
- **Default value for `Rcpp::Nullable<>` parameters.** `Rcpp::compileAttributes`
  doesn't always handle defaults gracefully; pass NULL explicitly from R
  rather than relying on a default.
- **Unnamed mask = OOB read.** See step 3.
- **Letting exceptions unwind through `minimize_direct<F>`.** See step 8.
  Catch at the closure boundary.
- **Failing to set `gc_finalize = true` on the XPtr** leaks the
  `std::function` on every call.
- **Threads.** Set thread count once in the wrapper before constructing
  the closure; do not toggle it inside the closure (it requires an R
  callback). xsdm calls `RcppParallel::setThreadOptions` once in
  `make_loglik_math_xptr_cpp` and never inside the closure.

## Reusable checklist

- [ ] `XptrClosureState` POD struct holding owned vectors + scratch
- [ ] `[[Rcpp::export]]` builder accepts `Rcpp::Nullable<NumericVector>` mask,
      `CharacterVector free_names`, and explicit `num_threads`, `grad`,
      `gradstep`
- [ ] Up-front validation: gradstep > 0, free_names non-empty, mask
      named or NULL
- [ ] Data snapshotted with `std::vector::assign`
- [ ] Lazy canonical-vector setup on first ObjFun call
- [ ] Pure-C++ canonical-name lookup (no SEXP after construction)
- [ ] `try { ... } catch (...) { f = +Inf; g = 0 }` at the closure boundary
- [ ] `return Rcpp::XPtr<ucminf::ObjFun>(fn.release(), true)`
- [ ] Caller passes `free_names = names(param_vector)`
- [ ] Parity tests vs. R reference at `tolerance = 1e-8`
- [ ] Regression tests for unnamed mask + bad gradstep
