---
editor_options: 
  markdown: 
    wrap: 72
---

# Plan: from xsdm-devel v1 to v2

This document is the working plan for upgrading the `xsdm-devel` package
from its current version 1 to version 2. It is a planning document —
nothing in it is final until Dan signs off — and it should be read
alongside the three sources cited below. It is intentionally explicit
and a bit verbose, because the user has asked that the plan be detailed
enough that very little judgement needs to be exercised at
implementation time.

All planning artefacts produced through this process must live under
`xsdm-devel/plan/`. No file outside that folder is to be modified during
planning.

------------------------------------------------------------------------

## 1. Sources and definitions

-   **source1** — `xsdm-devel/manual/xsdmmodel/xsdmModel.pdf` (compiled
    from the `xsdmModel.Rnw` in the same folder). This document has
    already been updated by Dan to describe the v2 mathematics. It is
    the authoritative reference for the mathematics of v2, including the
    new ltsgr+iv detection link, the parameter reductions that remove
    the two structural non-identifiabilities in that case, and the
    resulting biological parameter set.
-   **source2** —
    `~/Projects/xsdmMle/ClaudeOutputs/Planning/FormalV2Plan.md`. This is
    the formal v2 plan for a *different* branch of the package
    (`xsdmMle`, which targets n\>=2). It is **not** definitive for
    `xsdm-devel` v2, but section 2 ("Major data structures") is the
    inspiration for the `model_structure` / `params_bio` / `params_math`
    arguments adopted here. Sections 3 (detection links), 4 (function
    signatures), and 5 (validation rules) are also useful reference
    material.
-   **source3** — `~/Projects/xsdmMle/ClaudeOutputs/`. Reference R/C++
    prototypes of the n\>=2 case. Most relevant here are the docstrings
    and validation blocks at the top of the main `*_nxn.R` files, plus
    `detection_link_iv.R` and `detection_link.R`. The numeric machinery
    in source3 is mostly for the n\>=2 case and is out-of-scope for v2
    of `xsdm-devel`; do not port n\>=2 code wholesale.

When this document refers to "n\>1" or "n\>=2", we mean the
matrix-population case described in source2 §2. **None of the n\>=2
numerical code paths will be implemented in v2 of `xsdm-devel`.** v2
will, however, *accept and validate* arguments that describe the n\>=2
case, and then immediately throw a clear "not implemented" error so that
v3 can later supply the missing implementation without breaking the v2
API.

------------------------------------------------------------------------

## 2. Goals and scope of v2

### 2.1 What v2 adds relative to v1

1.  **A second detection link.** v1 supports only the ltsgr-only
    detection link of source1 §3. v2 adds support for the ltsgr-and-iv
    link of source1 §4. Either link can be chosen by the user.
2.  **A new argument convention** based on `model_structure`,
    `params_bio`, and `params_math`. These objects (described in §3
    below) replace the current pattern of passing `mu`, `sigltil`,
    `sigrtil`, `o_mat`, `ctil`, `pd` directly. The new convention is
    forward-compatible with the n\>=2 work planned for v3.
3.  **Forward-compatible validators.** v2 accepts and validates
    `model_structure` / `params_bio` arguments that describe n\>=2
    models. When such arguments are valid but currently unsupported
    (n\>=2 in v2, or `ltsgriv_method = "Tuljapurkar"` for n\>=2 in v2),
    the dispatching function emits a clear "not implemented in
    xsdm-devel v2" error and stops. This is so v3 can wire in the n\>=2
    path without changing any user-facing signatures.

### 2.2 What v2 must NOT change (the **no-deltas** invariant)

For every legal v1 call, the corresponding v2 call (using the new
arguments, with link 0 selected and n=1) must return numerically
identical values. "Identical" here means equal up to a documented
tolerance of `1e-8` (absolute) for log-likelihoods and probabilities,
and exact equality for integer/logical fields. This is enforced by the
no-deltas mechanism described in §6.

This invariant excludes only:

-   the **shape** of the value returned by `math_to_bio()`, which by
    construction must change to the `params_bio` shape (list of length 2
    in the n=1 case) — but the *numerical content* of that list, when
    restricted to fields that exist in v1, must match v1.
-   the **shape** of the value returned by `bio_to_math()`'s argument;
    the return value (a math-scale named numeric vector) is unchanged in
    shape *and* content for n=1 link 0.

Anywhere a v1 → v2 reshape is unavoidable, the plan calls it out
explicitly in §5 below.

### 2.3 Out of scope for v2

-   Any actual computation under n\>=2. (Validators only.)
-   Any actual computation under `ltsgriv_method = "Tuljapurkar"`.
    (Validators only — and only for n\>=2; `"Tuljapurkar"` is forbidden
    for n=1 per source2 §2.1.)
-   Renaming or restructuring of the package (no `_1x1` / `_nxn`
    suffixes are introduced; the dispatch happens internally, driven by
    `model_structure`). This differs from source2's "rough approach" in
    §5 of that document; we agree with Dan that the source2 suggestion
    was tentative.
-   Performance work beyond what is necessary to wire in link 1.
-   Vignette rewrites beyond the minimum needed to keep examples running
    and to introduce link 1.

### 2.4 Maintaining the R / C++ duality

Every hot-path function in v1 has both a thin R wrapper around a C++
implementation (the canonical export) and a pure-R reference (named with
the `_r` suffix, kept in `R/internals.R`). v2 must preserve this duality
for every function that is rewritten:

-   The pure-R reference must be updated to take the new arguments and
    to cover both detection links.
-   The C++ kernel must be updated to take the new arguments and to
    cover both detection links.
-   The parity test `test-<func>_r_vs_cpp.R` must continue to assert
    numerical agreement at tolerance 1e-6 (existing tolerance) for both
    link 0 and link 1.

Rationale (Dan's): the R reference is the readable specification of the
math; far more users read R than C++. We do not give that up just
because v2 is an opportunity to refactor.

------------------------------------------------------------------------

## 3. New data structures

This section is a precis of the data structures used by v2 in
`xsdm-devel`. The full validation rules are catalogued in source2 §5; we
re-state only those rules that apply in our reduced (n=1 implemented,
n\>=2 validated-only) setting. **Where this plan and source2 disagree on
naming or shape, this plan wins for `xsdm-devel`.** Open naming
questions are flagged in §10.

### 3.1 `model_structure`

A list. Its length determines `n`:

-   Length `n^2 + 3 = 4` ⇒ `n = 1` (the only computationally implemented
    case in v2).
-   Length `n^2 + 4` ⇒ `n >= 2` (validated, but the dispatching function
    will then immediately error out with "not implemented in xsdm-devel
    v2").

Element layout (all positional):

| pos | name | type / value |
|----|----|----|
| 1..n\^2 | matrix-entry indices | integer vector; strictly increasing, no duplicates; entries in `1:n_env`. Empty allowed only if n\>=2. |
| n\^2 + 1 | `ltsgriv_method` | scalar character: `"Lyapunov"` or `"Tuljapurkar"`. `"Tuljapurkar"` is forbidden when n=1. |
| n\^2 + 2 | `max_lag` | non-negative integer. Validated against `dim(env_dat)[2] - 3` when env_dat is available. For n=1 link 0 it is unused (kept for forward compatibility); for n=1 link 1 it is the truncation lag in the Bartlett-weighted iv estimator (eq. (\ref{eq:iv_comp}) of source1). |
| n\^2 + 3 | detection-link indicator | integer 0 or 1. `0` ⇒ ltsgr-only link (source1 §3); `1` ⇒ ltsgr-and-iv link (source1 §4). |
| n\^2 + 4 | `anchored_entry` | length-2 integer vector. Present iff `n >= 2`. Must be absent when `n = 1`. |

Notes specific to **n=1**:

-   `model_structure[[1]]` is a strictly-increasing integer vector of
    length `p` (here `p` is the number of environmental variables that
    act on the single matrix entry — the `p` of v1). It must be
    non-empty.
-   `model_structure[[2]] = ltsgriv_method`. For n=1 it must be
    `"Lyapunov"`.
-   `model_structure[[3]] = max_lag`. Used only when the detection-link
    indicator is `1`.
-   `model_structure[[4]] = detection-link indicator (0 or 1)`.
-   No `anchored_entry`.

### 3.2 `params_bio`

For **n=1** (the only computationally implemented case in v2),
`params_bio` is a list of length 2:

-   `params_bio[[1]]`: named list describing the dependence of the
    single matrix entry on the environment. Always contains:
    -   `mu` — numeric vector of length
        `p = length(model_structure[[1]])`.
    -   `sigltil` — positive numeric vector of length `p` (Inf allowed
        in boundary-model context).
    -   `sigrtil` — positive numeric vector of length `p` (Inf allowed
        in boundary-model context). Contains additionally `o_mat` if and
        only if `p > 1`, in which case `o_mat` is a `p x p` orthogonal
        matrix.
-   `params_bio[[2]]`: named list of detection-link parameters. The
    names depend on `model_structure[[4]]`:
    -   If `0`: `pd` (in (0, 1]) and `ctil` (real). (Same names as v1.)
    -   If `1`: `pd` (in (0, 1]), `gamma` (\>=0), `betahat` (\>=0).
        These are the post-reduction parameters of source1 eq.
        (\ref{eq:probdetect_iv}); they correspond to the symbols `pd`,
        `hat_gamma`, `beta` in source1 (the `hat`s and renaming are an
        open question — see §10).

For **n \>= 2**: list of length `n^2 + 1`, structured as described in
source2 §2.2.1. `xsdm-devel` v2 only needs the validation rules; no
computation is performed.

### 3.3 `params_math`

A named numeric vector. Its length and names are uniquely determined by
`model_structure` via a new function `mathscale_names(model_structure)`
(see §5.4). For n=1 link 0 the names and order match exactly what
`make_mask_names(p)` returns today, so existing v1 math-scale vectors
continue to work without renaming.

The transformations between math and bio scales are:

| bio | math | transformation |
|----|----|----|
| `mu` | `mu1..p` | identity |
| `sigltil` | `sigltil1..p` | `log` |
| `sigrtil` | `sigrtil1..p` | `log` |
| `o_mat` (p\>=2) | `o_par1..q` | matrix log of skew-symmetric basis (existing) |
| **link 0 only** |  |  |
| `ctil` | `ctil` | identity |
| `pd` | `pd` | `logit` (math → `expit` to bio) |
| **link 1 only** |  |  |
| `gamma` | `gamma` | `log` (gamma must be \>=0) |
| `betahat` | `betahat` | `log` (betahat must be \>=0) |
| `pd` | `pd` | `logit` |

The `gamma`, `betahat` math-scale names (and their non-negativity
transform choice) are tentative; see §10.

------------------------------------------------------------------------

## 4. Inventory of affected functions

This section enumerates every exported function in v1, classifies it,
and states what (if anything) must change in v2. The classifications
are:

-   **(A)** Signature changes (now takes `model_structure` /
    `params_bio` / `params_math`); behavior must be backwards-compatible
    for n=1 link 0.
-   **(B)** New function in v2 (no v1 equivalent).
-   **(C)** No change in v2 (utility, untouched).
-   **(D)** Internal helper repurposed; not in NAMESPACE.

| function | class | notes |
|----|----|----|
| `bio_to_math` | A | takes `model_structure` + `params_bio`. |
| `build_orthogonal_matrix` | C | unchanged. |
| `convert_equivalence_class` | C | unchanged (still operates on the n=1 sub-list of `params_bio`). |
| `create_mask` | A | takes `model_structure` instead of `p`. Names are produced by `mathscale_names()`. |
| `create_param_vector_masked` | A | takes `model_structure` instead of `p`. |
| `dist_between_params` | A | takes `model_structure`; `p1`/`p2` now in `params_bio` or `params_math` form. |
| `env_data_array` | C | unchanged. |
| `expit` | C | unchanged. |
| `extract_orthogonal_matrix_parameters` | C | unchanged. |
| `habitat_suitability` | A | takes `env_list`, `model_structure`, `params_bio`. |
| `interpret_parameters` | A | takes `model_structure`, `params_bio`. |
| `like_ltsg` | C | low-level kernel; unchanged. |
| `like_neg_ltsgr` | C | unchanged (it is link-agnostic — it computes `mean(-S_t/2) + log(lambda_max)`-like quantity, which is reused by both links). Re-examined in stage L1; may turn out to need a sibling helper for the iv side. |
| `log1mexp`, `log1pexp` | C | unchanged. |
| `log_prob_detect` | A | takes `env_dat`, `model_structure`, `params_bio`. Internal dispatch on detection-link indicator. |
| `loglik_bio` | A | takes `env_dat`, `occ`, `model_structure`, `params_bio`. |
| `loglik_math` | A | takes `params_math`, `env_dat`, `occ`, `model_structure`, `mask`, `negative`, `num_threads`. |
| `make_mask_names` | A→D | superseded by `mathscale_names()`. We retain `make_mask_names()` as an internal alias for one minor cycle so old code in tests keeps working, then remove it. See §5.4 for the rationale. |
| `math_to_bio` | A | takes `model_structure` and `params_math`; returns `params_bio`. |
| `num_par` | A→D | replaced by `length(mathscale_names(model_structure))`. Kept as an internal alias during the transition. |
| `num_env_var` | A→D | becomes useless once `mathscale_names()` exists; kept internal during the transition. |
| `optimize_likelihood` | A | takes `env_dat`, `occ`, `model_structure`, `mask`, ... |
| `profile_likelihood` | A | takes `env_dat`, `occ`, `model_structure`, `mask`, `params_math_optim`, `profile_parameter`, ... |
| `start_parms` | A | takes `env_dat`, `model_structure`, `mask`, `breadth`, `num_starts`. |
| `vsp` | A | takes `env_data`, `model_structure`, `params_bio`. |
| **new in v2:** |  |  |
| `mathscale_names` | B | replaces `make_mask_names`; returns canonical names from `model_structure`. |
| `validate_model_structure` | B | exported; runs all checks of §3 / source2 §5.1 + §5.4. |
| `validate_params_bio` | B | exported; takes `model_structure`, `params_bio`, optional `env_dat`. |
| `detection_link_ltsgr_iv` (or similar) | B | the new ltsgr-and-iv detection-link worker. Internal-only initially; export decision deferred. Implements eq. (\ref{eq:probdetect_iv}) of source1. |
| `iv_estimator` | B | internal helper that estimates infinitesimal variance from a sequence of per-year `S_t` values, via the Bartlett-weighted truncated autocovariance sum (eq. (\ref{eq:iv_comp}) of source1). Internal-only. |

In all (A) cases the v1 signature is removed. We do not maintain dual
signatures: v2 is a major version, the existing branch is `version2`,
and a controlled, mass cutover is much simpler than a long deprecation
lane. Any downstream code that relies on the old signatures will need to
be updated when it is rebuilt against v2; that is acceptable.

The internal-C++ kernels (`loglik_bio_cpp`, `loglik_math_cpp`,
`log_prob_detect_cpp`, `loglik_math_xptr_cpp`,
`make_loglik_math_xptr_cpp`, `.math_to_bio_cpp`,
`.build_canonical_param_vector_cpp`) all receive analogous changes;
per-kernel detail in §5.

------------------------------------------------------------------------

## 5. Function-by-function specification

The following sub-sections list, for each affected v1 function, the v2
signature, the v2 behaviour, and the changes to its R reference and C++
kernel.

Argument-order convention (per source2 §4.1): data first, structure
arguments second, parameter arguments third, algorithmic/control
arguments last. The exception is `loglik_math`, whose first argument is
the parameter vector so it can be passed to optimizers directly.

### 5.1 `loglik_bio`

```         
loglik_bio(env_dat, occ, model_structure, params_bio,
           return_prob = FALSE, sum_log_p = TRUE,
           num_threads = RcppParallel::defaultNumThreads())
```

Behaviour:

1.  Run `validate_model_structure(model_structure, env_dat)` and
    `validate_params_bio(model_structure, params_bio, env_dat)`.
    Validators throw on any inconsistency.
2.  Determine `n` from `length(model_structure)`. If `n >= 2`: stop with
    "n\>1 case not implemented in xsdm-devel v2".
3.  Determine link from `model_structure[[n^2 + 3]]` (which is element 4
    when n=1). Dispatch:
    -   link 0: existing computation, but now sourced from
        `params_bio[[1]]` (`mu`, `sigltil`, `sigrtil`, optional `o_mat`)
        and `params_bio[[2]]` (`pd`, `ctil`).
    -   link 1: new computation per source1 §4, using `iv_estimator`.
4.  Combine per-location log detection probabilities with `occ` exactly
    as v1 does: `log_p` is the per-location log P(detect). Likelihood
    reduction is `sum(occ * log_p + (1 - occ) * log1mexp(-log_p))`.
5.  Honour `return_prob` and `sum_log_p` exactly as v1.

R reference (`loglik_bio_r` in `R/internals.R`):

-   Updated in lockstep with the exported function.
-   Continues to use `log_prob_detect_r` for the slow path and the v1 R
    arithmetic for the reduction step.

C++ kernel (`loglik_bio_cpp` and supporting code in `src/`):

-   The `BioParams` struct in `src/loglik_bio.h` gains a tagged-union or
    two parallel structs (`BioParamsLink0`, `BioParamsLink1`) —
    implementer's choice; whichever is cleanest in C++17. Follow the
    prevailing style of the existing kernel.
-   The signature of the exported `[[Rcpp::export]] loglik_bio_cpp` is
    changed to take `model_structure` and `params_bio` directly, parsed
    in C++ (model_structure is a small `Rcpp::List`, `params_bio`
    likewise), or to take a flattened set of arguments produced by
    R-side helpers (whichever is cleaner — implementer's choice;
    document either way).
-   The `xsdm::loglik_bio_tile` core in `src/loglik_bio.h` gains a
    switch on the link indicator. The link-0 branch is the existing
    code, byte-for-byte. The link-1 branch calls a new core kernel
    implementing eq. (\ref{eq:probdetect_iv}).

Hot-path consideration: the R wrapper still routes the
`sum_log_p = TRUE && !return_prob` case directly to C++. The other
branches go through `log_prob_detect()` and reduce in R, as today.

### 5.2 `log_prob_detect`

```         
log_prob_detect(env_dat, model_structure, params_bio,
                return_prob = FALSE,
                num_threads = RcppParallel::defaultNumThreads())
```

Behaviour:

1.  Validators (as in `loglik_bio`).
2.  n=1 only; link dispatch as in `loglik_bio`.
3.  Returns a numeric vector of length `n_loc`:
    -   link 0: existing `log(pd) - log1pexp(ctil + h)`, where
        `h = like_neg_ltsgr(...)`.
    -   link 1: per source1 eq. (\ref{eq:probdetect_iv}),
        `log(pd) - log1pexp(beta - 4*gammahat/W + 2*Sbar/W)` where
        `Sbar = mean(hat_S_t)` and `W = iv_estimator(hat_S_t, max_lag)`.
        The hat scaling has already been absorbed into `hat_sigL`,
        `hat_sigR`, so in code we simply read `sigltil_hat` from
        `params_bio[[1]]` and use `like_neg_ltsgr_like` on those.
4.  Honour `return_prob`.

R reference (`log_prob_detect_r`): updated; same semantics as exported.

C++ kernel: new branch in `log_prob_detect_cpp` and the underlying
`src/log_prob_detect.h` core. The link-1 branch needs a per-location iv
computation; this is a new C++ helper `xsdm::iv_estimator_tile`, with a
companion R reference `iv_estimator_r` for parity tests.

### 5.3 `loglik_math`

```         
loglik_math(params_math, env_dat, occ, model_structure,
            mask = NULL, negative = TRUE,
            num_threads = RcppParallel::defaultNumThreads())
```

Behaviour:

1.  Validate `model_structure`. Validate `params_math` and `mask`
    together: union of names must equal
    `mathscale_names(model_structure)` exactly, `mask` may take `Inf`
    for `sigltil*`, `sigrtil*`, `pd` and (link 1 only) for `gamma`,
    `betahat`.
2.  `params_math` may arrive unnamed from a C++ optimizer callback. Same
    recovery logic as v1: assign canonical free-parameter names if
    length matches.
3.  Combine `params_math` with `mask` into a full math-scale vector.
4.  Convert to `params_bio` via
    `math_to_bio(model_structure, full_vec)`.
5.  Call `loglik_bio` with the same `model_structure`. Apply `negative`.

R reference: updated; routes through `math_to_bio_r` and `loglik_bio_r`.

C++ kernel: `loglik_math_cpp` now reads `model_structure` and dispatches
internally. The XPtr factory `make_loglik_math_xptr_cpp` is updated
analogously: the closure captures `model_structure` along with the data,
mask, and threading config.

### 5.4 `mathscale_names`, and the retirement of `make_mask_names` / `num_par` / `num_env_var`

```         
mathscale_names(model_structure)
```

Returns a character vector of canonical math-scale names for the model
described by `model_structure`. For n=1 link 0:
`c("mu1", ..., "mup", "sigltil1", ..., "sigltilp", "sigrtil1", ..., "sigrtilp", "ctil", "pd", "o_par1", ..., "o_par_q")`
— **same names and order as `make_mask_names(p)` returns today**, so
existing math-scale vectors interoperate unchanged. For n=1 link 1: as
above except `ctil` is replaced by `gamma`, `betahat`. For n\>=2: per
source2; not implemented in v2 except for validation.

`make_mask_names(p)`, `num_par(p)`, `num_env_var(n)` remain as
*non-exported* aliases during the transition (one stage), to keep the
test suite compiling. They are removed before v2 release.

### 5.5 `math_to_bio` / `bio_to_math`

```         
math_to_bio(model_structure, params_math)        # → params_bio
bio_to_math(model_structure, params_bio)         # → named numeric (params_math)
```

`math_to_bio` is the function whose return shape changes most visibly:
v1 returned a flat list with names `mu`, `sigltil`, `sigrtil`, `ctil`,
`pd`, `o_mat`. v2 returns `params_bio` (a list of length 2). The link-0
n=1 case is exactly the same parameters in a different shape:
`list(list(mu, sigltil, sigrtil, [o_mat]), list(pd, ctil))`. The
`bio_to_math` direction is unchanged in shape (still returns a named
numeric); the only difference is that the input is now `params_bio`.

R references: `math_to_bio_r` and the (non-existent v1 equivalent of)
`bio_to_math_r` are updated. Add a `bio_to_math_r` for parity testing
since today there is no R reference for `bio_to_math`.

C++ kernel: `.math_to_bio_cpp` becomes
`.math_to_bio_cpp(model_structure, params_math)` and returns the new
shape. `.build_canonical_param_vector_cpp` takes `model_structure`
instead of `p`.

### 5.6 `start_parms`

```         
start_parms(env_dat, model_structure, mask = NULL, breadth = 1,
            num_starts = 100L)
```

Behaviour:

1.  Validate `model_structure`.
2.  n=1 only; dispatch on detection-link indicator.
3.  Build a `range_df` whose row-names are exactly
    `mathscale_names(model_structure)` (minus any `mask` names). For
    link 0 the existing data-driven heuristics apply unchanged. For link
    1 we need new heuristics for `gamma` and `betahat` ranges and a
    different centering rule for `pd` / no `ctil` row. Heuristic details
    are listed in §10 (open question).
4.  Sobol' design via `pomp::sobol_design`, identical to v1.

R reference: not currently maintained; `start_parms` is a thin wrapper
around `get_range_df_` and `get_start_parms_`. We do not add a C++
implementation either. The function stays pure R.

### 5.7 `optimize_likelihood`

```         
optimize_likelihood(env_dat, occ, model_structure, mask = NULL,
                    num_starts = 100L, breadth = 1,
                    parallel = FALSE,
                    num_threads = RcppParallel::defaultNumThreads(),
                    control = list(), verbose = FALSE)
```

Behaviour: drop-in replacement for v1, plus the dispatch through
`model_structure`. `loglik_xptr` is now built from
`model_structure`-aware factory; otherwise the multi-start scaffolding
is unchanged.

### 5.8 `profile_likelihood`

```         
profile_likelihood(env_dat, occ, model_structure, mask = NULL,
                   params_math_optim, profile_parameter,
                   increment_left = 0.1, increment_right = increment_left,
                   num_steps_left = 20L, num_steps_right = num_steps_left,
                   alpha = 0.95,
                   num_threads = RcppParallel::defaultNumThreads(),
                   control = list(), verbose = FALSE)
```

Note: `optim_param_vector` is renamed `params_math_optim` to align with
source2 §4.9. The behaviour is identical for link 0 n=1; link 1 simply
uses a different set of math-scale free parameter names (no `ctil`, but
`gamma` and `betahat`).

### 5.9 `dist_between_params`

```         
dist_between_params(model_structure, p1, p2,
                    mask = NULL, give_closest_rep = FALSE)
```

Behaviour: as v1, with the dispatcher reading the link indicator out of
`model_structure`. The Hungarian distance treats `mu`, the widths, and
`o_mat` identically across links. The link-specific scalar parameters
(`ctil` for link 0; `gamma`, `betahat` for link 1) contribute additively
to the squared distance in the same simple way `ctil` and `pd` did in
v1.

Rationale: the equivalence-class symmetries in v1 are sign flips and
permutations of `o_mat` columns (with corresponding column-permutations
of `sigltil` / `sigrtil`). These symmetries do not interact with the
detection-link parameters, so the link 1 case factors cleanly.

Open question: do we need a different distance metric on the inverted
widths in link 1, given the additional `alpha` reduction? My current
read of source1 says *no* — the widths in link 1 are `hat_sig`
(post-reduction) and they enter the Hungarian cost identically to link
0.

### 5.10 `habitat_suitability` / `vsp`

```         
habitat_suitability(model_structure, params_bio, env_list,
                    output = "", overwrite = FALSE,
                    return_prob = TRUE, threads = 0L, wopt = list())

vsp(env_data, model_structure, params_bio,
    return_raster = FALSE)
```

Behaviour: routes to the same C++ kernel as `log_prob_detect` (under the
hood `habitat_suitability` already calls `log_prob_detect_cpp`
directly). The hot-tile loop in `habitat_suitability` builds a per-tile
`env_dat` chunk and forwards `model_structure` and `params_bio` to the
kernel.

`vsp` is a thin wrapper around `log_prob_detect` and stays that way.

### 5.11 `interpret_parameters`

```         
interpret_parameters(model_structure, params_bio, plot_indices,
                     plot_lims = NULL, env_dat = NULL, occ = NULL,
                     breadth = 1, ...)
```

Behaviour: the diagnostic plots interpret the *growth-environment*
function. That function does not depend on the choice of detection link,
so the inferred contour shape is computed identically in both link 0 and
link 1. The plotting code reads `params_bio[[1]]` for `mu`, `sigltil` /
`hat_sigltil`, `sigrtil` / `hat_sigrtil`, `o_mat`. The link-specific
parameters are not used for the plot.

There is one subtlety: source1 §5 (link 0) and §6 (link 1) explain that
the inferred function differs by an unknown affine transformation across
links, but the *shape* (contours) is what is biologically meaningful.
Since this function only plots the shape, the implementation does not
need to know the link.

### 5.12 `validate_model_structure` and `validate_params_bio`

These are new exported helpers. Their checks are exactly source2 §5.1,
§5.2, §5.3, §5.4, restricted to the union of n=1 and n\>=2 cases. The
n=1 sections override source2 only insofar as we explicitly accept the
link 1 parameter names `gamma` / `betahat` (or whatever names §10
settles).

Both validators take optional `env_dat`. Cross-consistency checks that
require `env_dat` (like `max_lag <= n_time - 3` and "indices in
`model_structure[[k]]` lie in `1:n_env`") run only when env_dat is
supplied; otherwise the function emits a `message()` saying that those
checks were skipped, mirroring the "DAN: ..." comment in source2 §5.4.

Validator failures use `stop()` with a clear message that names the
offending element. No partial validation: validators should run all
applicable checks and only abort on the first failure (per source2's
"raise a clear error at the first failure").

------------------------------------------------------------------------

## 6. The "no-deltas" verification component

The user's invariant is "for every legal v1 call, the equivalent v2 call
returns the same value, up to documented tolerance, in the n=1 link 0
case". The user proposed a per-function pair of scripts
(`<func>_nodeltas_pre.R` and `<func>_nodeltas_post.R`) producing three
RData files. That works, but I think a slightly tightened version is
cleaner:

### 6.1 Recommended structure

A single directory `xsdm-devel/plan/nodeltas/` (under plan/, since this
infrastructure lives only during the migration).

For each function `<func>` in class (A) of §4, there is a single
sub-folder:

```         
xsdm-devel/plan/nodeltas/<func>/
  ├── pre.R          # run on v1 (or with v1 functions in scope)
  ├── post.R         # run on v2
  ├── snapshot.rds   # produced by pre.R, consumed by post.R
  └── notes.md       # human notes about case coverage
```

-   `pre.R`:
    -   Sets a fixed seed (`set.seed(20260504)` is fine — today's date
        as YYYYMMDD).

    -   Calls a small generator that produces a list of test cases, each
        being a list of v1-style arguments. Cases must cover:

        -   the small-n regime (e.g. 50 locations, 20 timepoints, p=1) —
            used so failures are easy to debug;
        -   the realistic regime (size of `example_1$env_array`, p=2);
        -   boundary models: at least one case per `Inf` mask slot
            (sigl=Inf, sigr=Inf, pd=Inf) when the function under test
            allows masks;
        -   the non-default `mask` case (some parameters fixed) where
            applicable;
        -   the `return_prob = TRUE` and `sum_log_p = FALSE` flag
            combinations where applicable.

    -   For each case, calls v1 `<func>` and stores the result.

    -   Calls a converter `to_v2_args(<v1-args>)` (see §6.2) to produce
        the equivalent v2 arguments.

    -   Writes one rds bundle:

        ``` r
        saveRDS(list(cases_v1 = ..., cases_v2 = ..., results_v1 = ...,
                     metadata = list(R_version = ..., xsdm_version = "1.0.0",
                                     generated_at = Sys.time())),
                file = "snapshot.rds")
        ```
-   `post.R`:
    -   Loads `snapshot.rds`.
    -   For each case, calls v2 `<func>` with `cases_v2[[i]]`.
    -   Compares to `results_v1[[i]]` with
        `all.equal(..., tolerance = 1e-8)`, short-circuiting and
        printing a clear diff on failure.
    -   Reports a one-line PASS/FAIL summary at the end.

This is one snapshot file per function (vs. three in the user's
proposal), which is cleaner without losing any of the testing value.

### 6.2 The v1→v2 converter

This is a small, single file `plan/nodeltas/v1_to_v2.R` with one
converter per function:

``` r
to_v2_args_loglik_bio <- function(v1_args) {
  p <- length(v1_args$mu)
  list(
    env_dat = v1_args$env_dat,
    occ = v1_args$occ,
    model_structure = list(
      seq_len(p),         # element 1: matrix-entry indices
      "Lyapunov",         # ltsgriv_method
      0L,                 # max_lag (unused for link 0 n=1)
      0L                  # detection-link indicator: 0 = ltsgr-only
    ),
    params_bio = list(
      list(mu = v1_args$mu, sigltil = v1_args$sigltil,
           sigrtil = v1_args$sigrtil,
           o_mat = if (p > 1) v1_args$o_mat else NULL),  # NULL slot omitted
      list(pd = v1_args$pd, ctil = v1_args$ctil)
    ),
    return_prob = v1_args$return_prob,
    sum_log_p = v1_args$sum_log_p,
    num_threads = v1_args$num_threads
  )
}
```

(omitting the `o_mat` slot when p==1 per source2 §2.2.2; the validator
should reject `o_mat` when p==1.)

The converters are pure functions of the v1 arguments. They are tested
by `pre.R` itself (calling the v2 function once on `cases_v2[[1]]` is
*not* done in `pre.R` — `pre.R` runs against v1, where v2 doesn't
exist).

### 6.3 When the no-deltas check runs

For each step in §7:

1.  **Before** writing any code for the step, the engineer (or Claude)
    creates `pre.R`, `post.R`, and the converter for every function
    touched by the step. `pre.R` is run *now*, against v1, on a
    freshly-checked-out copy of the v1 code (e.g. on a worktree at
    `main`). The resulting `snapshot.rds` is committed to
    `plan/nodeltas/<func>/`.
2.  The engineer writes `post.R` and the v2 implementation in the same
    commit (or chain of commits).
3.  After implementation, the engineer runs `post.R` against the v2
    code. It must PASS.
4.  If `post.R` fails, the engineer fixes the v2 code (or, very rarely,
    the converter — never the snapshot).

Why I prefer this over `testthat`-style tests: the snapshots are large
and not particularly stable to refactors of *test* infrastructure (they
are RData, not source). Putting them in `plan/nodeltas/` rather than
`tests/testthat/` makes it clear they're a one-time migration scaffold,
to be deleted (or compressed into a smaller permanent regression test)
once v2 ships.

Alternatively, a single permanent regression test
`tests/testthat/test-v1_v2_parity.R` could be retained that loads a
single, smaller snapshot covering the realistic-regime case. This is
discussed in §10 as an open question.

### 6.4 What about the `_r_vs_cpp` parity tests?

Those tests live in `tests/testthat/` and are unrelated to no-deltas:
they assert that the R reference and the C++ kernel produce the same
answer at tolerance 1e-6. Because both are being updated together, the
parity tests must be updated alongside the implementation in each step,
and they will continue to run in `devtools::test()` long after v2 ships.

The two mechanisms are complementary: `_r_vs_cpp` ensures internal
consistency within v2; `nodeltas` ensures cross-version consistency
between v1 and v2.

------------------------------------------------------------------------

## 7. Staged plan

Each stage is a self-contained chunk of work that can be implemented,
tested, reviewed, and merged independently. Stages are mostly ordered by
dependency; see the dependency map at the end of this section.

### Stage 0: planning artefacts and scaffolding (this branch only)

-   **What**: write `plan1.md` (this file), commit. Add an empty
    `plan/nodeltas/` directory (with a `.gitkeep`) and an empty
    `plan/nodeltas/v1_to_v2.R` skeleton.
-   **Tests / checkpoints**: none beyond "the document compiles in a
    reader's head". Dan reviews and signs off, possibly producing
    `plan2.md` to capture revisions.
-   **Branch / commit / merge**: branch is `version2` (already exists).
    Single commit "plan: add v2 plan1.md" on `version2`. No merge to
    `main`.

### Stage 1: validators

-   **What**: add
    `validate_model_structure(model_structure, env_dat = NULL)` and
    `validate_params_bio(model_structure, params_bio, env_dat = NULL)`
    in `R/validate_model_structure.R` and `R/validate_params_bio.R`. Add
    internal helpers as needed. Wire `roxygen2` so they show up in
    `NAMESPACE`.
-   **Tests / checkpoints**:
    -   `tests/testthat/test-validate_model_structure.R` covers all
        checks of §3.1 / source2 §5.1 (every error path on its own, plus
        a minimum of happy-path n=1 link 0 / n=1 link 1 / n\>=2 cases).
    -   `tests/testthat/test-validate_params_bio.R` likewise.
-   **No-deltas**: not applicable (these are new functions; no v1
    equivalent).
-   **Branch / commit / merge**: sub-branch `version2-stage1-validators`
    off `version2`; merge back to `version2` when green.

### Stage 2: `mathscale_names` and the canonical names plumbing

-   **What**: add `mathscale_names(model_structure)` in
    `R/mathscale_names.R`. Implement the n=1 link 0 / link 1 cases and a
    best-effort n\>=2 case (the n\>=2 case is needed by
    `validate_params_bio` to count expected entries; doesn't need to be
    performance-critical). Re-implement `make_mask_names`, `num_par`,
    `num_env_var` as internal thin wrappers around `mathscale_names` so
    the rest of the code keeps compiling.
-   **Tests / checkpoints**:
    -   `tests/testthat/test-mathscale_names.R` covers all branches.
    -   Existing tests for `make_mask_names`, `num_par`, `num_env_var`
        continue to pass.
-   **No-deltas**: trivial; covered by the existing tests for
    `make_mask_names` etc., which now exercise the new implementation
    through the alias.
-   **Branch / commit / merge**: sub-branch
    `version2-stage2-mathscale-names`; merge back to `version2`.

### Stage 3: detection-link 1 worker functions (link-1 math infrastructure)

-   **What**: add `iv_estimator_r` (R reference) and
    `xsdm::iv_estimator_*` (C++) implementing the Bartlett-weighted
    truncated autocovariance sum of source1 eq. (\ref{eq:iv_comp}). Add
    `detection_link_iv_r` (R reference, possibly identical to
    `detection_link_iv` from source3) and the C++ analogue. These are
    not yet wired into `log_prob_detect` — Stage 5 does that. They live
    as internal helpers.
-   **Tests / checkpoints**:
    -   `tests/testthat/test-iv_estimator.R`: closed-form checks against
        short, hand-computed series. R-vs-C++ parity at tolerance 1e-6.
    -   `tests/testthat/test-detection_link_iv.R`: matches the source3
        `detection_link_iv` function's behaviour. R-vs-C++ parity.
-   **No-deltas**: not applicable (new helpers).
-   **Branch / commit / merge**: sub-branch
    `version2-stage3-link1-helpers`; merge back to `version2`.

### Stage 4: `math_to_bio` / `bio_to_math` reshape

-   **What**: change the signature of `math_to_bio` to take
    `model_structure, params_math` and return a `params_bio` list.
    `bio_to_math` likewise (input is now a `params_bio`, output
    unchanged in shape). Update `.math_to_bio_cpp` and
    `.build_canonical_param_vector_cpp`. Update R references. Add
    `bio_to_math_r` to the references list.
-   **No-deltas**: per §6, run `pre.R` for `math_to_bio` and
    `bio_to_math` *before* this stage starts. The shape difference is
    handled in the converter and assertion: the snapshot stores the
    v1-shape output, and `post.R` calls
    `flatten_params_bio_to_v1_list(v2_out)` (a small helper in
    `plan/nodeltas/v1_to_v2.R`) before comparing.
-   **Tests / checkpoints**:
    -   `test-math_to_bio.R`, `test-math_to_bio_cpp.R`,
        `test-math_to_bio_r_vs_cpp.R` updated to the new signature.
    -   `test-bio_to_math.R` updated.
    -   `plan/nodeltas/math_to_bio/post.R` PASS.
    -   `plan/nodeltas/bio_to_math/post.R` PASS.
-   **Branch / commit / merge**: sub-branch `version2-stage4-math-bio`;
    merge back to `version2`.

### Stage 5: `log_prob_detect`

-   **What**: change the signature; add the link-1 branch (using the
    Stage 3 helpers); update R reference; update C++ kernel.
-   **No-deltas**: pre-run before this stage. Snapshot covers the link-0
    n=1 case at p=1 and p=2, each with several seeds.
-   **Tests / checkpoints**:
    -   `test-log_prob_detect.R`, `test-log_prob_detect_cpp.R`,
        `test-log_prob_detect_r_vs_cpp.R` updated.
    -   Add `test-log_prob_detect_link1.R` covering link-1 specific
        cases (parity with the standalone `detection_link_iv` helper
        from Stage 3).
    -   `plan/nodeltas/log_prob_detect/post.R` PASS.
-   **Branch / commit / merge**: sub-branch
    `version2-stage5-log-prob-detect`; merge back to `version2`.

### Stage 6: `loglik_bio`

-   **What**: change the signature; thread the link dispatch through;
    update R reference; update C++ kernel.
-   **No-deltas**: pre-run before this stage. Snapshot covers all four
    flag combinations of `(return_prob, sum_log_p)` and both p=1 and
    p=2.
-   **Tests / checkpoints**:
    -   `test-loglik_bio.R`, `test-loglik_bio_cpp.R`,
        `test-loglik_bio_r_vs_cpp.R` updated.
    -   Add `test-loglik_bio_link1.R`.
    -   `plan/nodeltas/loglik_bio/post.R` PASS.
-   **Branch / commit / merge**: sub-branch
    `version2-stage6-loglik-bio`; merge back to `version2`.

### Stage 7: `loglik_math`, `create_param_vector_masked`, `create_mask`, XPtr factory

-   **What**: change all four signatures (the three R helpers plus the
    C++ XPtr factory) to take `model_structure`. Update R references and
    C++ kernels. The XPtr factory captures `model_structure` in the
    closure.
-   **No-deltas**: pre-run for `loglik_math`,
    `create_param_vector_masked`, `create_mask`. Snapshot for
    `loglik_math` covers `mask = NULL`, mask with one fixed mu, mask
    with Inf for sigl/sigr/pd.
-   **Tests / checkpoints**:
    -   The corresponding `test-*` files updated.
    -   `plan/nodeltas/loglik_math/post.R` PASS.
    -   `plan/nodeltas/create_*/post.R` PASS.
-   **Branch / commit / merge**: sub-branch
    `version2-stage7-loglik-math`; merge back to `version2`.

### Stage 8: `start_parms`, `optimize_likelihood`

-   **What**: update both for the new arguments. The `start_parms`
    heuristic for link 1 is new — see §10 for the open question on what
    ranges to use for `gamma` and `betahat`.
-   **No-deltas**: pre-run. `start_parms` snapshots are deterministic
    given the seed used inside `pomp::sobol_design`; the snapshot stores
    the entire returned tibble. `optimize_likelihood` is
    *non-deterministic* in the sense that floating-point noise across
    optimizer iterations can change the last digit of the optimum. The
    no-deltas check therefore uses tolerance `1e-4` for `loglik` and
    `1e-2` for the parameter-vector components, with a documented caveat
    in `plan/nodeltas/<func>/notes.md`.
-   **Tests / checkpoints**:
    -   `test-start_parms.R`, `test-optimize_likelihood.R`,
        `test-get_start_parms.R`, `test-get_range_df.R`,
        `test-num_starts_validation.R`, `test-optimize_helpers.R`
        updated.
    -   `plan/nodeltas/start_parms/post.R` PASS.
    -   `plan/nodeltas/optimize_likelihood/post.R` PASS (with the
        relaxed tolerances above).
-   **Branch / commit / merge**: sub-branch
    `version2-stage8-starts-opt`; merge back to `version2`.

### Stage 9: `profile_likelihood`

-   **What**: update for the new arguments. Internal logic is unchanged.
-   **No-deltas**: pre-run with one small case (1 step left, 1 step
    right) and one realistic-size case (the example from
    `?profile_likelihood`). Same loose tolerances as Stage 8.
-   **Tests / checkpoints**:
    -   `test-profile_likelihood.R` updated.
    -   `plan/nodeltas/profile_likelihood/post.R` PASS.
-   **Branch / commit / merge**: sub-branch `version2-stage9-profile`;
    merge back to `version2`.

### Stage 10: `dist_between_params`

-   **What**: update for the new arguments. The link-1 branch is new but
    small (different scalar parameters).
-   **No-deltas**: pre-run. Tolerance 1e-8 (deterministic).
-   **Tests / checkpoints**:
    -   `test-dist_between_params.R`, `test-solve_lsap_cpp.R` updated.
    -   `plan/nodeltas/dist_between_params/post.R` PASS.
-   **Branch / commit / merge**: sub-branch `version2-stage10-dist`;
    merge back to `version2`.

### Stage 11: `habitat_suitability` and `vsp`

-   **What**: update both for the new arguments. The hot-tile loop in
    `habitat_suitability` forwards `model_structure` and `params_bio` to
    the C++ kernel updated in Stage 5.
-   **No-deltas**: pre-run with one small
    (`example_1$bio01[1:10, 1:10]`) raster case to keep the snapshot
    small.
-   **Tests / checkpoints**:
    -   `test-habitat_suitability.R`, `test-vsp.R` updated.
    -   `plan/nodeltas/habitat_suitability/post.R` PASS.
    -   `plan/nodeltas/vsp/post.R` PASS.
-   **Branch / commit / merge**: sub-branch `version2-stage11-rasters`;
    merge back to `version2`.

### Stage 12: `interpret_parameters`

-   **What**: update signature. Internal logic unchanged in shape; reads
    fields out of `params_bio[[1]]` instead of a flat list.
-   **No-deltas**: not applicable in the strictest sense (the function
    is almost entirely side effects via plotting). Instead we add a
    `tests/manual/interpret_parameters_v1_v2_visual.R` script that draws
    the v1 plot and the v2 plot side by side; Dan eyeballs them.
    Internally, `interpret_parameters` does compute the contour `z`
    matrix; we expose a small internal `compute_contour_grid_` and add a
    `plan/nodeltas/interpret_parameters/post.R` that asserts numerical
    equality of `z` between v1 and v2 calls.
-   **Tests / checkpoints**:
    -   `test-interpret_parameters.R` updated.
    -   `plan/nodeltas/interpret_parameters/post.R` PASS.
-   **Branch / commit / merge**: sub-branch
    `version2-stage12-interpret`; merge back to `version2`.

### Stage 13: clean up `make_mask_names` / `num_par` / `num_env_var` shims

-   **What**: remove the internal aliases. Confirm no remaining callers.
    Update test files that still use them. Remove from NAMESPACE.
-   **Tests / checkpoints**: full `devtools::check()` and
    `devtools::test()` pass.
-   **Branch / commit / merge**: sub-branch `version2-stage13-aliases`;
    merge back to `version2`.

### Stage 14: examples, vignettes, docs, NEWS

-   **What**: update every example block in `R/*.R` (most use
    `example_1$par_vec` or analogous and need to be reshaped). Update
    `vignettes/` (we did not survey vignette contents; assume one or two
    vignettes need adjustment). Bump `DESCRIPTION` to `Version: 2.0.0`.
    Add a NEWS.md entry summarising the breaking-API change and the
    ltsgr+iv detection link.
-   **Tests / checkpoints**:
    -   `devtools::check(vignettes = TRUE)` passes.
    -   `devtools::test()` passes.
    -   `plan/nodeltas/*/post.R` all pass (final aggregate run).
-   **Branch / commit / merge**: sub-branch `version2-stage14-docs`;
    merge back to `version2`.

### Stage 15: integration regression and final review

-   **What**: from the now-fully-updated `version2`, run a single final
    regression that calls every public function in v2 with the v1 test
    data and a reasonable v2 wrapping. This is essentially the union of
    all `plan/nodeltas/*/post.R` runs, executed end-to-end via a single
    `plan/nodeltas/run_all.R` driver. Document any deltas that exceed
    documented tolerances.
-   **Tests / checkpoints**: as above; produce
    `plan/nodeltas/final_report.md`.
-   **Branch / commit / merge**: this stage produces only documentation
    artefacts; commit on `version2`.

### Dependency map (informal)

```         
Stage 0 ─→ Stage 1 ─→ Stage 2 ─→ Stage 3 ─→ Stage 4 ─→ Stage 5 ─→ Stage 6 ─→ Stage 7
                                                                              ↓
              ┌───────────────────────────────────────────────────────────────┤
              ↓                                                                ↓
           Stage 8                                                          Stage 9
              ↓                                                                ↓
              └─→ Stage 10 ─→ Stage 11 ─→ Stage 12 ─→ Stage 13 ─→ Stage 14 ─→ Stage 15
```

Stage 1 (validators) and Stage 3 (link-1 helpers) could in principle be
done in parallel since neither depends on the other. We do them
sequentially in the stage list to keep the linear narrative simple, but
nothing in Stage 3 requires Stage 1, so a fast-track parallel approach
is acceptable.

------------------------------------------------------------------------

## 8. Branch / commit / merge conventions

-   **Branch:** `version2` is the trunk for all v2 work (already
    exists). Each stage gets a sub-branch named
    `version2-stage<N>-<short-tag>`, branched off `version2`.
    Sub-branches are short-lived: typically one sitting of work, one
    merge back, then deleted locally. Force-pushes are allowed on
    sub-branches up until the merge; never on `version2` and never on
    `main`.
-   **Commits:** small and atomic, with informative subject lines
    starting with the stage tag (e.g.
    `stage5: log_prob_detect signature & R wrapper`,
    `stage5: log_prob_detect link-1 C++ branch`,
    `stage5: log_prob_detect parity tests for link-1`). Avoid the "fixed
    everything" mega-commit pattern.
-   **Merging back to `version2`:** standard `git merge --no-ff` so each
    stage shows up as a merge commit. This is helpful for `git log`
    after the fact.
-   **`main`:** untouched until v2 is ready to release. Final merge to
    `main` is a single PR that goes through `devtools::check()`,
    `plan/nodeltas/run_all.R`, and Dan's review. After merge, `version2`
    branch can be deleted.
-   **Tags:** tag `v1.0.0` on `main` *before* any v2 work merges back,
    if not already tagged. Tag `v2.0.0` on the v2-merge commit.
-   **No-deltas snapshots are committed** to `version2` as part of each
    stage. They are versioned data; do not regenerate them after the
    initial commit unless the v1 reference itself changes.

------------------------------------------------------------------------

## 9. Testing summary

-   **Unit tests** (`tests/testthat/`): every `(A)`-class function has
    its v2 unit tests; every `(B)`-class function has new unit tests;
    every `*_r_vs_cpp` parity test is updated.
-   **No-deltas tests** (`plan/nodeltas/*/post.R`): one per `(A)`-class
    function; they pass at the documented tolerances.
-   **Vignette compilation** (`devtools::check(vignettes = TRUE)`):
    green in stage 14 onward.
-   **R CMD check** (`devtools::check()`): green from stage 13 onward.
    Until then, expected warnings/errors arise from the alias shims and
    half-converted vignettes; the engineer documents these in commit
    messages.
-   **Manual visual inspection**: only `interpret_parameters` requires
    this. A short visual-comparison script lives in
    `tests/manual/interpret_parameters_v1_v2_visual.R`.

------------------------------------------------------------------------

## 10. Open questions and things to keep developing

These are the items that still need a decision from Dan before
implementation begins (or, in some cases, before the corresponding stage
begins).

### 10.1 Link-1 parameter naming

`xsdm-devel` v2 needs concrete names for the link-1 detection-link
parameters in both the bio-scale `params_bio[[2]]` list and the
math-scale `params_math` vector. Source1 §4 calls the post-reduction
parameters `pd`, `beta`, and `hat_gamma`. Source2 §2.2.1 (n\>=2 case)
calls them `pd`, `alpha`, `beta`, `gamma` (pre-reduction). The user's
prose in `cp1.txt` calls them `pd`, `gamma`, and `betahat`.

My recommendation:

-   Bio-scale (post-reduction, used in `params_bio[[2]]` for n=1): `pd`,
    `gammahat`, `betahat`. The `hat`s match source1's notation for
    "after the reduction was applied"; they distinguish these from the
    pre-reduction `gamma` and (more importantly) the pre-reduction
    `beta` which is used for the n\>=2 input convention in source2 (and
    in v3).
-   Math-scale: `pd`, `gammahat`, `betahat` (i.e. `gammahat` is on the
    log scale, `betahat` is on the log scale, `pd` is on the logit
    scale).

But — alternative, and arguably cleaner because it does not put `hat` in
identifier names — drop the `hat` and just use `gamma`, `beta`, since
for n=1 the only entries one ever sees are post-reduction. Dan to
decide.

### 10.2 Should `start_parms` produce link-1 ranges?

We need data-driven heuristics for `gamma` and `betahat` ranges. A
defensible default:

-   `betahat` math-scale range `log(c(1e-2, 1, 1e2)) = c(-4.6, 0, 4.6)`
    (i.e. centred at `betahat = 1`, breadth ±2 orders of magnitude).
-   `gammahat` math-scale range derived from the empirical scale of
    `Sbar` and `W` at the centred starting point: pick `gammahat` such
    that `4 * gammahat / W ≈ Sbar / W`, then take ± a factor.

These are guesses; before Stage 8 we should run a quick sensitivity
analysis on a virtual species to confirm reasonable optimization
behaviour from these ranges.

### 10.3 Do we keep a permanent v1↔v2 regression test?

§6.3 suggests a single small `tests/testthat/test-v1_v2_parity.R` that
loads a compact "permanent" snapshot. Pros: catches accidental v2
behavioural drift in CI long after v2 ships. Cons: maintenance burden,
and once v3 lands we'll be carrying around a v1 snapshot forever.

My weak preference: keep one tiny permanent snapshot
(`tests/testthat/_snaps/v1_v2_parity.rds`, ≤ 50KB) for the example_1
realistic case only. Delete all other snapshots after the v2 release.

### 10.4 Validator behaviour when `env_dat` is missing

Source2 §5.4 has a "DAN: ..." comment about what to do when `env_dat`
isn't supplied at validation time (some checks become impossible). My
current spec says: emit a `message()` and skip. Alternatives: (a)
require `env_dat`; (b) warn instead of message; (c) make `env_dat`
optional but loud-fail when validation that *would have caught* a
problem turns out to have been needed downstream. Recommend (a) for
almost all paths — the validators are usually called from inside a
function that already has `env_dat` in hand — and (b) for the few
`bio_to_math` / `math_to_bio` calls that genuinely don't have `env_dat`.

### 10.5 Detection-link indicator: integer or string?

Source2 §5.1.6 has another "DAN: ..." comment: 0/1 vs `"ltsgr"` /
`"ltsgr_iv"`. Recommend descriptive strings: `"ltsgr"` and `"ltsgr_iv"`.
Strings make `model_structure` self-documenting and remove the
implicit-coercion hazard of integer 0/1 vs logical FALSE/TRUE. Cost:
every dispatch site does a string compare instead of an integer compare;
trivial. Spec note: replace "integer 0/1" with "character `\"ltsgr\"` or
`\"ltsgr_iv\"`" throughout this plan if Dan agrees.

### 10.6 Should `loglik_bio` and friends accept a "pre-validated" flag?

Validators run on every call. For the multi-start optimizer (\~100
starts × thousands of iterations each), validating `model_structure` and
`params_bio` at every loglik call is wasted work. Two options:

-   Accept it: profile shows the validation cost is \< 1% of total.
-   Add a `validate = TRUE` flag whose default is `TRUE` but which the
    optimizer scaffolding sets to `FALSE` after a one-time validation at
    the top of `optimize_likelihood`.

I lean towards "add the flag" because it's cheap to add and lets us keep
the validators thorough without worrying about their cost.

### 10.7 What about `like_neg_ltsgr`?

This subsection corrects an imprecision in earlier drafts of the plan
and flags an implementation choice that Stage 5 will have to make.

**What `like_neg_ltsgr` actually does today.** Despite the name, it does
*not* return `-ltsgr + log(lambda_max)`. Per location it returns

$$h \;=\; \tfrac{1}{2}\,\mean_{t}\!\left[\left(\frac{u_{t}}{\tilde\sigma(u_{t})}\right)^{2}\right]
       \;=\; \tfrac{1}{2}\,\mean{\tilde S_{t}},$$

i.e. half the per-location time-mean of `S_t` from source1 eq.
(\ref{eq:St}), evaluated with whatever widths the caller supplies. In v1
those widths are the link-0 post-reduction `\tilde\sigma_{L/R}`, and the
returned `h` is precisely the quantity fed into the link-0 detection
formula as `log(pd) - log1pexp(ctil + h)`. So `h` is "the part of
`-ltsgr` that survives the link-0 reduction", not `-ltsgr` itself.

**What link 1 needs from this kernel area.** Source1 eq.
(\ref{eq:probdetect_iv}) (in the new tilde-only notation, post Stage 0
notation cleanup) requires two per-location quantities:

-   **(a)** `\mean{\tilde S_t}` — the per-location time mean, same shape
    as what `like_neg_ltsgr` already produces, just evaluated on the
    link-1 post-reduction widths `\tilde\sigma = \sqrt{\alpha}\,\sigma`.
-   **(b)** `\tilde W` — the long-run variance of `\tilde S_t`,
    estimated via the Bartlett-weighted truncated autocovariance sum
    (eq. (\ref{eq:iv_comp})). This is **new**: `like_neg_ltsgr` only
    returns the mean, not the variance, and computing `\tilde W`
    requires access to the full per-year time series of `\tilde S_t`.

So `like_neg_ltsgr`'s **R-level signature and return value do not need
to change** to support link 1 — its existing output covers (a) directly.
What's new is (b): we need a sibling helper that produces the
per-location long-run variance of `\tilde S_t`. (My earlier "leave it
alone" recommendation stands at the user-facing level; class C in §4 is
unchanged.)

**The implementation choice (Stage 5).** The sibling helper can be
implemented in two ways, and Stage 5 must pick one:

-   **Path A — independent kernel.** Add a new C++ helper
    `iv_estimator_tile` that walks `env_dat` from scratch and produces
    `\tilde W` per location. `like_neg_ltsgr` and the underlying
    `like_ltsg` C++ kernel are byte-for-byte unchanged. *Cost:* link 1
    loglik evaluations sweep `env_dat` twice per call (once for the mean
    via `like_neg_ltsgr`, once for the variance via
    `iv_estimator_tile`).
-   **Path B — single-pass refactor.** Either (i) introduce a
    lower-level kernel that exposes the per-year `\tilde S_t` grid, on
    top of which both the existing `like_neg_ltsgr` and the new iv
    estimator are thin reductions; or (ii) extend the C++ core of
    `like_neg_ltsgr` so that, when given a "with-iv" flag, it returns
    both `mean` and `W` in a single pass. *Cost:* more refactor work;
    requires care so the link-0 hot path is not slowed by the link-1
    capability.

Path A is simpler and is the safe default. Path B is worth the effort
only if profiling in Stage 5 shows that the second sweep is a material
fraction of optimizer time. Recommendation: start with Path A, benchmark
a realistic optimization, and switch to Path B only if the double-sweep
cost is meaningful.

Whichever path is chosen, `like_neg_ltsgr` itself stays class C and
needs no no-deltas snapshot of its own. The new sibling kernel(s) get
their own R-vs-C++ parity test (per §2.4) and a closed-form unit test on
a hand-computed short series (per Stage 3 of §7).

### 10.8 Is there a use-case for fitting *both* links and selecting the

best by AIC?

Not in v2, but the API is shaped so this is straightforward to add
later: just call `optimize_likelihood` with two different
`model_structure` objects and compare. No plan changes needed.

### 10.9 `o_mat` slot in n=1 `params_bio[[1]]` when p=1

Source2 §2.2.2 says `o_mat` must be omitted when `p == 1`. Our
validators must enforce this; the existing `bio_to_math` / `math_to_bio`
code-paths already handle "no `o_par` entries when p=1", but the
validator must reject a list that *contains* an `o_mat` slot when p=1.
Implementation note for stage 1.

### 10.10 Should `model_structure` be an S3 class?

Right now we treat `model_structure` as a plain list. We could promote
it to an S3 class `xsdm_model_structure` with a constructor and a
`print.xsdm_model_structure` method, which would (a) make printing nicer
and (b) make it harder for users to accidentally pass the wrong list
shape. Cons: a small style decision that ripples through every function.
Recommend deferring to v2.x or v3; not blocking for v2.0.

### 10.11 Forward-compatibility edge cases for the `n>=2` validator path

Source2 §5.2 specifies a richer `params_bio` shape for `n>=2`. Our
validator must accept that shape and then (in the dispatch step of each
function) immediately error out with "n\>1 case not implemented".
Question: do we want the *validator* itself to error on n\>=2, so users
discover the limitation early? Recommend **no** — let
`validate_params_bio` accept n\>=2 as structurally valid, and let the
dispatch step of each function be the place where "not implemented" is
signalled. This keeps the validators reusable in v3 without
modification.

### 10.12 What does "boundary model" look like in link 1?

In v1, boundary models are encoded by setting `mask` entries for
`sigltil*` / `sigrtil*` / `pd` to `Inf`. In link 1 we have `gammahat` /
`betahat` / `pd` as the new detection-link parameters. For `pd` the
existing `Inf`-on-math-scale (= 1 on bio-scale) convention works
unchanged. For `gammahat` and `betahat`: do we *want* boundary models
that send these to `+Inf`? Probably not — their bio-scale interpretation
as multiplicative scale parameters makes infinite-bound boundaries
pathological. Recommend rejecting `Inf` for these in `mask` validation.
Confirm during Stage 7 implementation.

### 10.13 Performance regression watchlist

Items to keep an eye on, mostly logged here so they don't get lost:

-   Validators called inside the optimizer hot path (see §10.6).
-   The link-1 iv estimator: for a 39-year time series at 2000
    locations, we're doing one autocovariance sum per location per
    loglik evaluation. Order is `n_loc * max_lag * ts_length`, which for
    `max_lag = 4` and `ts_length = 39` is \~300k ops per loglik —
    negligible compared to the n_loc \* ts_length matrix exponentiation
    already done in the link-0 path. Should be fine; sanity-check in
    Stage 5.
-   The XPtr factory in `loglik_math_xptr_cpp` captures more state (the
    parsed `model_structure`) — needs care to ensure no per-call R
    round-trip.

------------------------------------------------------------------------

## 11. Summary checklist

A condensed tick-list of what "v2 is done" looks like:

-   [ ] All exported v1 functions have v2 signatures per §5.
-   [ ] No exported v1 signature survives.
-   [ ] `model_structure` and `params_bio` validators are exported and
    used by every dispatching function.
-   [ ] R reference + C++ kernel duality preserved for every hot-path
    function.
-   [ ] All `*_r_vs_cpp` parity tests pass at tolerance 1e-6.
-   [ ] All `plan/nodeltas/*/post.R` scripts pass at their documented
    tolerances.
-   [ ] `devtools::check()` passes.
-   [ ] `devtools::check(vignettes = TRUE)` passes.
-   [ ] DESCRIPTION bumped to `Version: 2.0.0`.
-   [ ] NEWS.md describes the breaking change and the new detection
    link.
-   [ ] Vignettes and examples are reshaped to use `model_structure` /
    `params_bio`.
-   [ ] All Stage-1-through-Stage-15 sub-branches are merged back to
    `version2`.
-   [ ] `version2` is merged to `main` in a single reviewed PR.
-   [ ] Tag `v2.0.0` lands on `main`.

------------------------------------------------------------------------

*End of plan1.md.*
