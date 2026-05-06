---
editor_options: 
  markdown: 
    wrap: 72
---

# Plan: from xsdm-devel v1 to v2 (revision 2)

This is the second revision of the working plan for upgrading the
`xsdm-devel` package from its current version 1 to version 2. It
supersedes plan1.md and is intended to stand on its own; it does not
refer back to plan1.md, but it does still refer to the three sources
listed in §1.

The plan is intentionally explicit and a bit verbose, because it has
been requested that the plan be detailed enough that very little
judgement needs to be exercised at implementation time.

All planning artefacts must live under `xsdm-devel/plan/`. No file
outside that folder is to be modified during planning, except where this
plan explicitly authorises a notational update of source1 (which has
already been carried out for sections 4 and 6 of source1 prior to this
revision).

------------------------------------------------------------------------

## 1. Sources and definitions

-   **source1** — `xsdm-devel/manual/xsdmmodel/xsdmModel.pdf` (compiled
    from `xsdmModel.Rnw` in the same folder). Updated by Dan to describe
    the v2 mathematics, including (since the latest notational pass) a
    tilde-only convention for the post-reduction parameters and
    time-series of the link that takes both ltsgr and iv:
    `\tilde{\gamma}`, `\tilde{\sigma}_{L/R}`, `\tilde{S}_t`,
    `\tilde{W}`. Source1 is the authoritative reference for the
    mathematics of v2.
-   **source2** —
    `~/Projects/xsdmMle/ClaudeOutputs/Planning/FormalV2Plan.md`. This is
    the formal v2 plan for a *different* branch of the package
    (`xsdmMle`, which targets n\>=2). It is **not** definitive for
    `xsdm-devel` v2, but section 2 ("Major data structures") is the
    inspiration for the `model_structure` / `params_bio` / `params_math`
    arguments adopted here. Sections 3 (detection links), 4 (function
    signatures), and 5 (validation rules) are useful reference material;
    this document overrides source2 wherever they disagree.
-   **source3** — `~/Projects/xsdmMle/ClaudeOutputs/`. Reference R/C++
    prototypes of the n\>=2 case. Most relevant here are the docstrings
    and validation blocks at the top of the main `*_nxn.R` files, plus
    `detection_link.R`. The numeric machinery in source3 is mostly for
    n\>=2 and is out-of-scope for v2 of `xsdm-devel`; do not port n\>=2
    numerics wholesale. The validation code in source3 is, however, a
    useful starting point for the validators in v2 of `xsdm-devel`,
    since v2 validators must accept and check arbitrary
    n.  

When this document refers to "n\>1" or "n\>=2" we mean the
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
    `params_bio`, and `params_math`. These objects (described in §4
    below) replace the current pattern of passing `mu`, `sigltil`,
    `sigrtil`, `o_mat`, `ctil`, `pd` directly. The new convention is
    forward-compatible with the n\>=2 work planned for v3.
3.  **Forward-compatible validators and `mathscale_names`.** The
    validators (`validate_model_structure`, `validate_params_bio`) and
    `mathscale_names` are required to handle arbitrary n in v2, not just
    n=1. The intent is that these three functions become in v2 exactly
    what they will need to be in v3 when the n\>=2 case is implemented
    numerically — so they should not need to change at the v2-to-v3
    transition. This makes the v2 data-structure contract a stable
    public surface that v3 can build on.
4.  **A new function `Sttil_mean_lrv`** that computes per-location time
    means of `\tilde{S}_t` (and, optionally, the long-run variance of
    `\tilde{S}_t`). This function replaces and extends the existing
    `like_neg_ltsgr` (see §5 and §7).

### 2.2 What v2 must NOT change (the **no-deltas** invariant)

For every legal v1 call, the corresponding v2 call (using the new
arguments, with link `"ltsgr"` selected and n=1) must return numerically
identical values. "Identical" here means equal up to a documented
tolerance of `1e-8` (absolute) for log-likelihoods and probabilities,
and exact equality for integer/logical fields. This is enforced by the
no-deltas mechanism described in §6.

This invariant excludes only:

-   the **shape** of the value returned by `math_to_bio()`, which by
    construction must change to the `params_bio` shape (list of length 2
    in the n=1 case) — but the *numerical content* of that list, when
    restricted to fields that exist in v1, must match v1.
-   the **shape** of `bio_to_math()`'s argument; the return value (a
    math-scale named numeric vector) is unchanged in shape *and* content
    for n=1 link `"ltsgr"`.

Anywhere a v1 → v2 reshape is unavoidable, the plan calls it out
explicitly in §5.

### 2.3 Out of scope for v2

-   Any actual computation under n\>=2. (Validators and
    `mathscale_names` only.)
-   Any actual computation under `ltsgr_iv_method = "Tuljapurkar"`.
    (Validators only — and only for n\>=2; `"Tuljapurkar"` is forbidden
    for n=1 per source2 §2.1.)
-   Renaming or restructuring of the package along `_1x1`/`_nxn` lines.
    The dispatch happens internally, driven by `model_structure`. This
    differs from the "rough approach" tentatively suggested in source2
    §5.
-   Vignette rewrites beyond the minimum needed to keep examples running
    and to introduce the new detection link.

### 2.4 Maintaining the R / C++ duality

Every hot-path function in v1 has both a thin R wrapper around a C++
implementation (the canonical export) and a pure-R reference (with the
`_r` suffix, kept in `R/internals.R`). v2 must preserve this duality for
every function that is rewritten:

-   The pure-R reference must be updated to take the new arguments and
    to cover both detection links.
-   The C++ kernel must be updated to take the new arguments and to
    cover both detection links.
-   The parity test for each duality pair must continue to assert
    numerical agreement at tolerance 1e-6 (existing tolerance) for both
    the `"ltsgr"` and `"ltsgr_iv"` link types.

Rationale: the R reference is the readable specification of the math;
far more users read R than C++. We do not give that up just because v2
is an opportunity to refactor.

------------------------------------------------------------------------

## 3. Notation conventions used in this plan

To reduce ambiguity:

-   We refer to the two detection links by the canonical strings
    `"ltsgr"` (the v1 ltsgr-only link of source1 §3) and `"ltsgr_iv"`
    (the new ltsgr-and-iv link of source1 §4).
-   Math symbols follow source1 in its current form. In particular the
    *final* post-reduction link-`"ltsgr_iv"` parameters are
    `\tilde{\gamma}`, `\tilde{\sigma}_{L/R}`, `\tilde{S}_t`,
    `\tilde{W}`, with `\hat{\gamma}` being the *intermediate* reduction
    (used inside derivations only, not exposed at the user level).
-   In R code, `\tilde{\gamma}` is named `gammatil`. This and the other
    bio-scale parameter names (and their math-scale conversions) are
    fixed in §4.

------------------------------------------------------------------------

## 4. New data structures

### 4.1 `model_structure`

A list whose elements have positional meaning. Its length determines
`n`:

-   Length `n^2 + 3 = 4` ⇒ `n = 1` (the only computationally implemented
    case in v2).
-   Length `n^2 + 4` ⇒ `n >= 2` (validated, but the dispatching function
    will then immediately error out with "not implemented in xsdm-devel
    v2").

The positional layout is:

| pos | informal name | type / value |
|-------------|-------------|---------------------------------------------|
| 1..n\^2 | matrix-entry indices | integer vector; strictly increasing, no duplicates; entries in `1:n_env`. Empty integer vector allowed only if n\>=2. These cannot *all* be empty. |
| n\^2 + 1 | `ltsgr_iv_method` | scalar character: `"Lyapunov"` or `"Tuljapurkar"`. `"Tuljapurkar"` is forbidden when n=1. |
| n\^2 + 2 | `max_lag` | non-negative integer. Validated against `dim(env_dat)[2] - 3` when env_dat is supplied. For n=1 link `"ltsgr"` it is unused (kept for forward compatibility); for n=1 link `"ltsgr_iv"` it is the truncation lag of the Bartlett-weighted long-run-variance estimator (eq. (\ref{eq:iv_comp}) of source1). |
| n\^2 + 3 | `detection_link_type` | scalar character: `"ltsgr"` (the link of source1 §3) or `"ltsgr_iv"` (the link of source1 §4). |
| n\^2 + 4 | `anchored_entry` | length-2 integer vector. Present **iff** n\>=2. **Must be absent or present but equal to NA when n=1; presence and non-NA at n=1 is a validation error.** |

Notes specific to **n=1**:

-   `model_structure[[1]]` is a strictly-increasing integer vector of
    length `p` (the number of environmental variables that act on the
    single matrix entry — the `p` of v1). It must be non-empty.
-   `model_structure[[2]] = ltsgr_iv_method`. For n=1 it must be
    `"Lyapunov"`.
-   `model_structure[[3]] = max_lag`. Used only when
    `detection_link_type == "ltsgr_iv"`.
-   `model_structure[[4]] = detection_link_type`.
-   No `anchored_entry` element, or if present it must be NA. If element
    5 is present and non-NA at n=1, validation fails with a clear
    message.

The list elements may be named (`"ltsgr_iv_method"`, `"max_lag"`,
`"detection_link_type"`, `"anchored_entry"`) or unnamed; both forms are
accepted. The validator does not require names. Naming is a useful
convention for readability and is recommended in user-facing
constructors and examples.

The validator must check the names if they are present. If one or more
of them is present but not equal to the positionally expected name,
throw an error. So the positions drive what kind of information is
expected in each spot, and any given name must follow that if it is
present.

### 4.2 `params_bio`

For **n=1** (the only computationally implemented case in v2),
`params_bio` is a list of length 2:

-   `params_bio[[1]]`: named list describing the dependence of the
    single matrix entry on the environment. Always contains:
    -   `mu` — numeric vector of length
        `p = length(model_structure[[1]])`.
    -   `sigltil` — positive numeric vector of length `p` (Inf allowed
        in boundary-model context; see §10).
    -   `sigrtil` — positive numeric vector of length `p` (Inf allowed
        in boundary-model context).
    -   Contains additionally `o_mat` if and only if `p > 1`, in which
        case `o_mat` is a `p × p` orthogonal matrix. **When p=1, `o_mat`
        must be absent**; presence is a validation error.
-   `params_bio[[2]]`: named list of detection-link parameters. Names
    depend on `detection_link_type`:
    -   If `"ltsgr"`: `pd` (in (0, 1]) and `ctil` (real). Same names as
        v1.
    -   If `"ltsgr_iv"`: `pd` (in (0, 1]), `gammatil` (≥0), `beta` (≥0).
        These are the post-reduction parameters of source1 eq.
        (\ref{eq:probdetect_iv}) in the current tilde-only notation.

For **n \>= 2**: list of length `n^2 + 1`, structured as described in
source2 §2.2.1. `xsdm-devel` v2 only needs the validation rules; no
computation is performed.

### 4.3 `params_math`

A named numeric vector. Its length and names are uniquely determined by
`model_structure` via `mathscale_names(model_structure)` (see §5.4). For
n=1 link `"ltsgr"` the names and order match exactly what the v1
function `make_mask_names(p)` returns today, so existing v1 math-scale
vectors continue to interoperate without renaming.

The transformations between math and bio scales are:

| bio | math | transformation from math to bio scales |
|-------------------------|------------------------|-----------------------|
| `mu` | `mu1..p` | identity |
| `sigltil` | `sigltil1..p` | `exp` |
| `sigrtil` | `sigrtil1..p` | `exp` |
| `o_mat` (p\>=2) | `o_par1..q` | matrix exponential of skew-symmetric matrix the lower-triangle of which is the o_par1..q values |
| **link `"ltsgr"` only** |  |  |
| `ctil` | `ctil` | identity |
| `pd` | `pd` | `expit` |
| **link `"ltsgr_iv"` only** |  |  |
| `gammatil` | `gammatil` | `exp` |
| `beta` | `beta` | `exp` |
| `pd` | `pd` | `expit` |

The names are deliberately the same on both scales — `gammatil` is
called `gammatil` everywhere, and `beta` is called `beta` everywhere —
with the math-scale value being the log of the bio-scale value (for
`gammatil`, `beta`, the `sig` widths) or the logit (for `pd`). This
means no `_log` / `_logit` suffixes appear in identifier names.

### 4.4 Validity of `model_structure` and `params_bio`

Both objects are subject to a number of structural constraints. The full
check list is essentially that of source2 §5 (with the following
overrides specified by this document):

-   `detection_link_type` is a string `"ltsgr"` or `"ltsgr_iv"` (not a
    0/1 integer).
-   For n=1 the n\^2+4 element (`anchored_entry`) **must be absent**;
    its presence is a validation error.
-   For n=1 with p=1, `params_bio[[1]]$o_mat` **must be absent**; its
    presence is a validation error.
-   For n=1 link `"ltsgr_iv"`, `params_bio[[2]]` has names `pd`,
    `gammatil`, `beta` (not `alpha`/`beta`/`gamma` as in source2's n\>=2
    spec, because in n=1 the parameter reductions remove `alpha` and
    `\hat{\gamma}` is folded into `\tilde{\gamma}`).
-   The validators are **required to handle arbitrary n** so that they
    can be reused without modification in v3. See §5.12 for details.

------------------------------------------------------------------------

## 5. Inventory and function-by-function spec

This section enumerates **every function in the package** — exported R
user-facing functions, internal R helpers and references, Rcpp-exported
C++ entry points, and C++ namespace-scope inline helpers — and
classifies each by what happens to it in v2.

### 5.1 Classification

-   **(A)** Signature change. The function continues to exist in v2 but
    takes different arguments — typically `model_structure` /
    `params_bio` / `params_math` instead of the v1 individual parameter
    arguments. Behavior must be numerically backwards- compatible on the
    n=1 link `"ltsgr"` overlap (the no-deltas invariant of §2.2).
-   **(B)** New function in v2 (no v1 equivalent).
-   **(C)** No control-flow / behavioural change. The function may
    receive cosmetic-only adjustments during Stage v1px (renames for
    identifier consistency, doc fixes, dead-code removal), but the
    signature, semantics, and computed result are unaffected by v2. For
    exported functions this means the v1 user-facing API of the function
    survives v2 unchanged.
-   **(D)** Deleted. Removed from the package source as part of v2. A→D
    means the function exists in v1, is touched (or shadowed) during the
    v2 staging, and is removed before v2 ships.

The four tables below cover, in turn, exported R user-facing functions
(Tier 1), internal R helpers and pure-R references (Tier 2),
Rcpp-exported C++ entry points called by R wrappers (Tier 3), and C++
namespace-scope inline helpers (Tier 4). For each row, "stage" is the v2
stage in which the change happens.

#### Tier 1: exported R user-facing functions (current NAMESPACE)

| function | class | stage | notes |
|---------------------|---------|-------|--------------------------------|
| `bio_to_math` | A | v1px + 4 | Doc-and-error-handling cleanup in Stage v1px (§7.1.b); signature changes to `(model_structure, params_bio)` in Stage 4. |
| `build_orthogonal_matrix` | C | — | Not changed. (Possible Stage v1px cosmetic tweak only.) |
| `convert_equivalence_class` | C | — | Not changed. Continues to operate on the n=1 sub-list of `params_bio` (same shape as the v1 list). |
| `create_mask` | A→D | 7 | Deleted in Stage 7. Job folded into the new internal helper `assemble_params_math_()`. |
| `create_param_vector_masked` | A→D | 7 | Deleted in Stage 7. Replaced by `assemble_params_math_(model_structure, params_math, mask)`, which wraps the C++ kernel `.build_canonical_param_vector_cpp` (whose signature changes in Stage 7 to take `model_structure`). |
| `dist_between_params` | A | 10 | Takes `model_structure`; `p1`/`p2` in `params_bio` or `params_math` form. |
| `env_data_array` | C | — | Not changed. |
| `expit` | C | — | Not changed. |
| `extract_orthogonal_matrix_parameters` | C | — | Not changed. |
| `habitat_suitability` | A | 11 | Takes `env_list`, `model_structure`, `params_bio`. |
| `interpret_parameters` | A | 12 | Takes `model_structure`, `params_bio`, plus an optional internal `compute_contour_grid_()` extracted from `compute_y` for no-deltas testing. |
| `like_ltsg` | C | — | The low-level Rcpp-exported kernel. Not changed in v2 (signature unchanged); see Tier 3 row for the conditional-on-Stage-3a-Path-B case. |
| `like_neg_ltsgr` | A→D | 3a + 3b | Stage 3a leaves it alone (new `Sttil_mean_lrv` lands alongside). Stage 3b deletes it and rewrites the few internal callers. |
| `log1mexp` | C | — | Not changed. |
| `log1pexp` | C | — | Not changed. |
| `log_prob_detect` | A | 5 | Takes `env_dat`, `model_structure`, `params_bio`. Internal dispatch on `detection_link_type`. |
| `loglik_bio` | A | 6 | Takes `env_dat`, `occ`, `model_structure`, `params_bio`. |
| `loglik_math` | A | 7 | Takes `params_math`, `env_dat`, `occ`, `model_structure`, `mask`, ... |
| `make_mask_names` | A→D | 7 | Stays exported and unchanged through Stage 6, becomes a thin wrapper around `mathscale_names` from Stage 2 onward, deleted in Stage 7. |
| `math_to_bio` | A | 4 | Takes `model_structure` and `params_math`; returns `params_bio` (shape change). |
| `num_par` | A→D | 7 | Same lifecycle as `make_mask_names`. |
| `num_env_var` | A→D | 7 | Same lifecycle as `make_mask_names`. |
| `optimize_likelihood` | A | 8 | Takes `env_dat`, `occ`, `model_structure`, `mask`, ... |
| `profile_likelihood` | A | 9 | Takes `env_dat`, `occ`, `model_structure`, `mask`, `params_math_optim`, `profile_parameter`, ... |
| `start_parms` | A | 8 | Takes `env_dat`, `model_structure`, `mask`, `breadth`, `num_starts`. |
| `vsp` | A | 11 | Takes `env_data`, `model_structure`, `params_bio`. |
| **new in v2:** |  |  |  |
| `mathscale_names` | B | 2 | Returns canonical math-scale names from `model_structure`. **Required to work for arbitrary n.** |
| `validate_model_structure` | B | 1 | Exported. **Required to work for arbitrary n.** |
| `validate_params_bio` | B | 1 | Exported. Takes `model_structure`, `params_bio`, optional `env_dat`. **Required to work for arbitrary n.** |
| `Sttil_mean_lrv` | B | 3a | New replacement for `like_neg_ltsgr`. R wrapper + C++ kernel + R reference (see Tiers 2 and 3). See §5.13. |

#### Tier 2: internal R helpers and pure-R references

These are the non-exported R functions inside the package. They split
into "pure-R references" (named `<func>_r`, used only by `_r_vs_cpp`
parity tests) and "private helpers" (used internally by other R
functions).

| function | location | class | stage | notes |
|-----------------|--------------|--------|-------|--------------------------|
| `log_prob_detect_r` | `R/internals.R` | A | 5 | Pure-R reference for `log_prob_detect`. Tracks the exported function's signature and link dispatch. |
| `loglik_bio_r` | `R/internals.R` | A | 6 | Pure-R reference for `loglik_bio`. |
| `loglik_math_r` | `R/internals.R` | A | 7 | Pure-R reference for `loglik_math`. |
| `math_to_bio_r` | `R/internals.R` | A | 4 | Pure-R reference for `math_to_bio`. |
| `create_param_vector_masked_r` | `R/internals.R` | A→D | 7 | Deleted alongside the public counterpart. Its validation logic moves into `assemble_params_math_()`. |
| `like_neg_ltsgr_r` | `R/like_neg_ltsgr_r.R` | A→D | 3b | Deleted alongside `like_neg_ltsgr`. |
| `like_neg_ltsgr_cpp` | `R/like_neg_ltsgr_cpp.R` | A→D | 3b | Internal alias to `like_neg_ltsgr` (kept only for back-compatibility with internal callers and tests). Deleted in Stage 3b. |
| `dist_between_params_r` | `R/internals.R` | A | 10 | Hungarian-algorithm pure-R reference for `dist_between_params`. |
| `distance_between_params_r` | `R/internals.R` | C | — | Brute-force (enumerate all flips/perms) reference, used to cross-check `dist_between_params_r` on small examples. Operates on biological-scale lists; not affected by the math-scale signature changes. |
| `check_env_array` | `R/internals.R` | C | — | Internal env_dat shape validator. Not changed. |
| `logit` | `R/internals.R` | C | — | Internal utility. Not changed. |
| `permutations` | `R/internals.R` | C | — | Internal utility (used by `distance_between_params_r`). Not changed. |
| `auto_plot_lims_` | `R/interpret_parameters.R` | A | 12 | Internal helper for `interpret_parameters`. Its `param_list` argument becomes `params_bio`. |
| `compute_contour_grid_` | `R/interpret_parameters.R` | B | 12 | New internal helper extracted from the existing `compute_y` closure inside `interpret_parameters`, so the contour grid `z` is computable without drawing (needed for the no-deltas snapshot). |
| `get_range_df_` | `R/start_parms.R` | A | 8 | Range builder. Takes `model_structure`; gains link-`"ltsgr_iv"` heuristics for `gammatil` / `beta`. |
| `get_start_parms_` | `R/start_parms.R` | C | — | Sobol-design wrapper. Takes a ranges data-frame whose row-names are already the right canonical names; logic unchanged. (Possible v1px rename `numstarts` → `num_starts`.) |
| `profile_one_side_` | `R/profile_likelihood.R` | A | 9 | Internal helper. Threads `model_structure` through. |
| `optimize_loglik_math_` | `R/optimize_likelihood.R` | A | 8 | Internal single-start runner. Threads `model_structure` through; constructs the XPtr factory with `model_structure` captured. |
| `resolve_xptr_grad_control_` | `R/optimize_likelihood.R` | C | — | Small ctrl-merging helper; not affected. |
| `xsdmStartupMessage` | `R/zzz.R` | C | — | Package-startup banner; not affected. |
| **new in v2:** |  |  |  |  |
| `assemble_params_math_` | `R/internals.R` | B | 7 | Replaces `create_param_vector_masked` (and absorbs `create_mask`'s job). Validates that names of `params_math` and `mask` together equal `mathscale_names(model_structure)` with no overlap and no missing slots; delegates the merge to the C++ kernel. |
| `Sttil_mean_lrv_r` | `R/internals.R` | B | 3a | Pure-R reference for `Sttil_mean_lrv`. Used in the parity test. |

#### Tier 3: Rcpp-exported C++ entry points (callable from R)

These are the C++ functions surfaced to R via `[[Rcpp::export]]`. The R
wrappers in Tier 1/Tier 2 ultimately call into one of these.

| C++ function | location | class | stage | notes |
|-----------------|---------------|--------|--------|---------------------|
| `.build_orthogonal_matrix_cpp` | `src/expm_skew.cpp` | C | — | Not changed. |
| `like_ltsg` | `src/like_ltsg.cpp` | C\* | 3a? | Not changed by default. **Conditional**: if Stage 3a chooses Path B (single-pass refactor) for `Sttil_mean_lrv`, this kernel may be refactored to expose the per-year `S_t` grid as a first-class output. Default is Path A, in which case this is class C. |
| `log_prob_detect_cpp` | `src/log_prob_detect.cpp` | A | 5 | Takes `model_structure` + `params_bio`; link dispatch. |
| `loglik_bio_cpp` | `src/loglik_bio.cpp` | A | 6 | Takes `model_structure` + `params_bio`. |
| `loglik_math_cpp` | `src/loglik_math.cpp` | A | 7 | Takes `model_structure`. |
| `make_loglik_math_xptr` | `src/loglik_math_xptr.cpp` | D | v1px or 7 | Older R-callback variant of the XPtr factory. Used only by its own parity test (`tests/testthat/test-loglik_math_xptr_parity.R`); the comment in `R/profile_likelihood.R:69` notes that it has a bug that prevents it from working in `profile_likelihood`. Recommended for deletion in Stage v1px (or, at the latest, Stage 7) along with its parity test. |
| `make_loglik_math_xptr_cpp` | `src/loglik_math_xptr_cpp.cpp` | A | 7 | The pure-C++ XPtr factory used by `optimize_likelihood`. Captures `model_structure` in the closure; closure body uses `validate = FALSE` semantics. |
| `.math_to_bio_cpp` | `src/math_to_bio.cpp` | A | 4 | Takes `model_structure`. Return shape changes to match `params_bio`. |
| `.build_canonical_param_vector_cpp` | `src/math_to_bio.cpp` | A | 7 | Takes `model_structure` instead of `p`. The R-side helper `assemble_params_math_()` is what calls this. |
| `.solve_lsap_cpp` | `src/solve_lsap.cpp` | C | — | Hungarian-algorithm linear-sum-assignment kernel. Not changed. |
| **new in v2:** |  |  |  |  |
| `Sttil_mean_lrv_cpp` | `src/Sttil_mean_lrv.cpp` (new) | B | 3a | C++ kernel for the new replacement of `like_neg_ltsgr`. Returns time-mean-of-`\tilde{S}_t` per location; optionally also long-run variance via Bartlett-weighted truncated autocovariance sum (single sweep when feasible). |

#### Tier 4: C++ namespace-scope inline helpers (header-only)

These are the inline `xsdm::*` functions and structs in `src/*.h`,
called from the Tier 3 entry points and from each other. They are not
directly callable from R.

| C++ entity | location | class | stage | notes |
|---------------|-----------|--------|--------|----------------------------|
| `xsdm::BioParams` (struct) | `src/math_to_bio.h` | A | 6 | Extended to hold the link-`"ltsgr_iv"` parameters as well; either tagged-union or two-struct approach (implementer's choice during Stage 6). |
| `xsdm::loglik_bio_tile` | `src/loglik_bio.h` | A | 6 | Switches on link type; `"ltsgr"` branch is the existing code, `"ltsgr_iv"` branch is new. |
| `xsdm::loglik_math_eval` | `src/loglik_math.h` | A | 7 | Takes `model_structure`-derived data. |
| `xsdm::math_to_bio_apply` | `src/math_to_bio.h` | A | 4 | Takes `model_structure`-derived dispatch info. |
| `xsdm::expit_cpp` | `src/math_to_bio.h` | C | — | Not changed. |
| `xsdm::log1mexp_cpp` | `src/loglik_bio.h` | C | — | Not changed. |
| `xsdm::log1pexp` (inline in cpp) | `src/log_prob_detect.h` | C | — | Not changed. |
| `xsdm::num_par_cpp` | `src/math_to_bio.h` | A→D | 7 | Replaced by a `mathscale_names`-derived count, computed R-side once per call; deleted alongside its R counterpart. |
| `xsdm::build_orthogonal_matrix_cpp` | `src/expm_skew.h` | C | — | Not changed. |
| `xsdm::expm_pade13`, `mm_cm`, `axpy_mat`, `set_scaled_identity`, `inf_norm`, `gauss_solve` | `src/expm_skew.h` | C | — | Low-level matrix-exponential utilities. Not changed. |
| **new in v2:** |  |  |  |  |
| `xsdm::Sttil_mean_lrv_tile` | `src/Sttil_mean_lrv.h` (new) | B | 3a | C++ kernel header for `Sttil_mean_lrv_cpp`. Computes per-location mean and (optionally) long-run variance of `\tilde{S}_t`. |

#### Files (data, generated, build)

For completeness, the package also contains data files, auto- generated
code, and build artefacts. These are not "functions" but are listed here
so the inventory is genuinely exhaustive:

-   `R/example_1.R`, `R/example_2.R`, `R/example_3.R` — data-doc files.
    Class **C** (the data is unchanged); the example objects they
    document have their parameter-vector / parameter-list fixtures
    reshaped to match v2 in Stage 13.
-   `R/RcppExports.R`, `src/RcppExports.cpp` — auto-generated by
    `Rcpp::compileAttributes()`. Regenerated automatically as Tier 3
    signatures change; not manually edited. Class **C** in the sense
    that no human edits are required.
-   `R/xsdm-package.R` — package-level docs; class **C** (cosmetic
    updates in Stage 13 if needed).
-   `R/zzz.R` — package-load hooks; class **C**.
-   `src/Makevars`, `src/Makevars.win`, `src/vendor/` — build
    configuration and vendored xtensor/xsimd/xtl headers. Class **C**.

#### Summary counts

After the v2 staging is complete:

-   **Tier 1 (exported R)**: 26 v1 entries → 24 v2 entries. Six
    deletions (`make_mask_names`, `num_par`, `num_env_var`,
    `create_mask`, `create_param_vector_masked`, `like_neg_ltsgr`) and
    four additions (`mathscale_names`, `validate_model_structure`,
    `validate_params_bio`, `Sttil_mean_lrv`). The remaining 20 are kept,
    of which most receive a v2 signature change.
-   **Tier 2 (internal R)**: roughly 19 v1 entries → 18 v2 entries.
    Three deletions (`create_param_vector_masked_r`, `like_neg_ltsgr_r`,
    `like_neg_ltsgr_cpp`) and three additions (`assemble_params_math_`,
    `compute_contour_grid_`, `Sttil_mean_lrv_r`).
-   **Tier 3 (Rcpp-exported)**: 10 v1 entries → 10 or 11 v2 entries
    depending on the Stage 7 decision about `make_loglik_math_xptr` (the
    older R-callback variant), and Stage 3a adds one
    (`Sttil_mean_lrv_cpp`).
-   **Tier 4 (C++ namespace-scope)**: roughly 12 v1 entries → 12 v2
    entries, with one deletion (`xsdm::num_par_cpp`) and one addition
    (`xsdm::Sttil_mean_lrv_tile`).

In all (A) cases the v1 signature is removed. We do not maintain dual
signatures: v2 is a major version, the existing branch is `version2`,
and a controlled mass cutover is much simpler than a long deprecation
lane.

### 5.2 The `validate = TRUE` argument convention

Many functions in v2 perform validation of `model_structure` and/or
`params_bio` as their first action. Validation is cheap relative to a
single hot-path log-likelihood evaluation, but it is not free, and
inside the multi-start optimizer (≈100 starts × thousands of iterations
each) repeated re-validation is wasteful. We therefore adopt the
following convention:

-   Every (A)-class function listed above accepts an additional argument
    `validate = TRUE`. (Where the function already has many arguments,
    this one is added at the end of the algorithmic- arguments group,
    before `num_threads` if present.)
-   When `validate = TRUE`, the function runs all its applicable
    structural validation checks at entry.
-   When `validate = FALSE`, the function trusts its inputs and skips
    the structural checks. It still performs cheap shape / type coercion
    that is needed for correct execution (e.g. `as.integer` on `occ`),
    but does not run `validate_model_structure`, `validate_params_bio`,
    `mathscale_names`-name-set checks, or similar.
-   When function `f` calls another (A)-class function `g` internally,
    `f` always calls `g` with `validate = FALSE`, having already
    validated its own inputs. This applies recursively.

Concretely, that means:

-   `optimize_likelihood` validates once at the top, then passes
    `validate = FALSE` to every internal call to `loglik_math`.
-   `loglik_math` (when `validate = FALSE`) trusts its caller's
    validation, then calls `math_to_bio` and `loglik_bio` with
    `validate = FALSE`.
-   `loglik_bio` (when `validate = FALSE`) calls `log_prob_detect` with
    `validate = FALSE` for the slow path.
-   Likewise `profile_likelihood` validates once and calls `loglik_math`
    with `validate = FALSE` thereafter.
-   `habitat_suitability` validates once and forwards `validate = FALSE`
    to the per-tile kernel.
-   And so on for every other dispatcher → kernel chain.

The XPtr factory `make_loglik_math_xptr_cpp` captures `validate = FALSE`
in the closure produced for the optimizer hot path; the factory itself
validates at construction time.

This convention should be present in the documentation of every
(A)-class function. The default `TRUE` keeps the user-facing API safe;
the `FALSE` value is an internal optimization that informed callers can
request.

### 5.3 `loglik_bio`

```         
loglik_bio(env_dat, occ, model_structure, params_bio,
           return_prob = FALSE, sum_log_p = TRUE,
           validate = TRUE,
           num_threads = RcppParallel::defaultNumThreads())
```

Behaviour:

1.  If `validate`, run
    `validate_model_structure(model_structure,     env_dat)` and
    `validate_params_bio(model_structure, params_bio,     env_dat)`.
    Validators throw on any inconsistency.
2.  Determine `n` from `length(model_structure)`. If `n >= 2`: stop with
    "n\>1 case not implemented in xsdm-devel v2".
3.  Read `detection_link_type` from `model_structure`. Dispatch:
    -   `"ltsgr"`: existing v1 computation, sourced from
        `params_bio[[1]]` (`mu`, `sigltil`, `sigrtil`, optional `o_mat`)
        and `params_bio[[2]]` (`pd`, `ctil`).
    -   `"ltsgr_iv"`: new computation per source1 §4 (using
        `Sttil_mean_lrv` with `compute_lrv = TRUE`).
4.  Reduce per-location log detection probabilities with `occ` exactly
    as v1 does (`sum(occ * log_p + (1 - occ) * log1mexp(-log_p))`),
    honouring `return_prob` and `sum_log_p` exactly as v1.

R reference (`loglik_bio_r` in `R/internals.R`): updated in lockstep,
covers both link types, and routes through `log_prob_detect_r` for the
slow path.

C++ kernel: `loglik_bio_cpp` and `xsdm::loglik_bio_tile` are extended to
switch on the link type. The `"ltsgr"` branch is the existing code,
byte-for-byte.

### 5.4 `log_prob_detect`

```         
log_prob_detect(env_dat, model_structure, params_bio,
                return_prob = FALSE,
                validate = TRUE,
                num_threads = RcppParallel::defaultNumThreads())
```

Behaviour:

1.  If `validate`, run validators.
2.  n=1 only; dispatch on `detection_link_type`.
3.  Returns a numeric vector of length `n_loc`:
    -   `"ltsgr"`: call `Sttil_mean_lrv(..., compute_lrv = FALSE)` to
        get per-location `mean(\tilde{S}_t)`; then compute
        `h = 0.5 * mean_Sttil`; then return `log(pd) - log1pexp(ctil + h)`.
    -   `"ltsgr_iv"`: call
        `Sttil_mean_lrv(..., compute_lrv = TRUE,     max_lag = ...)` to
        get both `mean_Sttil` and `lrv_Sttil` per location; then
        return
        `log(pd) - log1pexp(beta - 4*gammatil/lrv_Sttil + 2*mean_Sttil/lrv_Sttil)`.
4.  Honour `return_prob`.

R reference: `log_prob_detect_r` updated; same semantics as exported.

C++ kernel: `log_prob_detect_cpp` updated to call `Sttil_mean_lrv`'s C++
entry point with the appropriate `compute_lrv` flag and to apply the
corresponding closed-form detection formula.

There is **no separate `detection_link_iv` function** in v2. The
`"ltsgr_iv"` detection formula is computed inline inside
`log_prob_detect` (and its R reference and C++ kernel), exactly as the
v1 `"ltsgr"` formula is computed inline today.

### 5.5 `loglik_math`

```         
loglik_math(params_math, env_dat, occ, model_structure,
            mask = NULL, negative = TRUE,
            validate = TRUE,
            num_threads = RcppParallel::defaultNumThreads())
```

Behaviour:

1.  If `validate`, run
    `validate_model_structure(model_structure,     env_dat)`. Then
    validate `params_math` and `mask` together: union of names must
    equal `mathscale_names(model_structure)` exactly, `mask` may take
    `Inf` only for permitted boundary slots (`sigltil*`, `sigrtil*`,
    `pd`).
2.  `params_math` may arrive unnamed from a C++ optimizer callback. Same
    recovery logic as v1: assign canonical free-parameter names if
    length matches.
3.  Combine `params_math` with `mask` into a full math-scale vector.
4.  Convert to `params_bio` via
    `math_to_bio(model_structure,     full_vec, validate = FALSE)`.
5.  Call `loglik_bio(..., validate = FALSE)` with the same
    `model_structure`. Apply `negative`.

R reference: routes through `math_to_bio_r` and `loglik_bio_r`.

C++ kernel: `loglik_math_cpp` reads `model_structure` and dispatches
internally. The XPtr factory `make_loglik_math_xptr_cpp` is updated
analogously: the closure captures `model_structure` along with the data,
mask, and threading config, and is constructed with `validate = FALSE`
semantics for the inner-loop body.

### 5.6 `mathscale_names`, and the retirement of `make_mask_names` / `num_par` / `num_env_var`

```         
mathscale_names(model_structure)
```

Returns a character vector of canonical math-scale names for the model
described by `model_structure`. **Required to work for any n \>= 1.**
For n=1 link `"ltsgr"`:
`c("mu1", ..., "mup", "sigltil1", ..., "sigltilp", "sigrtil1", ..., "sigrtilp", "ctil", "pd", "o_par1", ..., "o_par_q")`
— same names and order as `make_mask_names(p)` produces today, so
existing math-scale vectors interoperate unchanged. For n=1 link
`"ltsgr_iv"`: as above except `ctil` is replaced by `gammatil`, `beta`.
For n\>=2: per source2 §2.3, with the renaming of `gamma` to `gammatil`
(where applicable) and the boundary conventions of this plan.

The intent is that `mathscale_names` becomes in v2 exactly what it will
need to be in v3 when the n\>=2 case is implemented numerically; it
should not need to change at the v2-to-v3 transition. The n\>=2 logic
can be ported from source3.

`make_mask_names(p)`, `num_par(p)`, `num_env_var(n)` remain *exported*
thin aliases from Stage 2 through Stage 6, so the v1 export surface is
undisturbed and intermediate stages keep compiling without a flag-day
rewrite of every internal call site. `create_mask` and
`create_param_vector_masked` likewise remain exported through Stage 6.
**All five are deleted entirely in Stage 7** — both the exports and the
source files — alongside the rest of the math-scale-helper
consolidation. After Stage 7, the only canonical-names entry point is
`mathscale_names(model_structure)`, and the only full-vector-assembly
entry point is the internal helper `assemble_params_math_()`.

This is a deliberate consolidation: v1 carries five exported helpers
that are layered on top of one another but produce essentially the same
canonical-name machinery; v2 carries one (`mathscale_names`).

### 5.7 `math_to_bio` / `bio_to_math`

```         
math_to_bio(model_structure, params_math,
            validate = TRUE)                      # → params_bio

bio_to_math(model_structure, params_bio,
            validate = TRUE)                      # → named numeric (params_math)
```

`math_to_bio` is the function whose return shape changes most visibly:
v1 returned a flat list with names `mu`, `sigltil`, `sigrtil`, `ctil`,
`pd`, `o_mat`. v2 returns `params_bio` (a list of length 2). The
link-`"ltsgr"` n=1 case is the same numeric content in a different
shape: `list(list(mu, sigltil, sigrtil, [o_mat]), list(pd, ctil))`.

`bio_to_math` returns a math-scale named numeric vector — same shape as
v1, only the input is now `params_bio`.

R references: `math_to_bio_r` is updated. A new `bio_to_math_r` is added
to `R/internals.R` for parity testing (in v1 there was no R reference
for `bio_to_math`; this is one of the cleanups, see §7.1).

The `bio_to_math` documentation and error handling are simplified per
§7.1: `bio_to_math` now relies entirely on `expm::logm` to compute the
math-scale `o_par` entries from `o_mat`, and any `expm::logm` error is
propagated unmodified to the user. The docstring states this plainly and
stops attempting to enumerate the conditions under which `expm::logm`
may fail.

C++ kernels: `.math_to_bio_cpp` and `.build_canonical_param_vector_cpp`
take `model_structure` instead of `p`.

### 5.8 `start_parms`

```         
start_parms(env_dat, model_structure, mask = NULL, breadth = 1,
            num_starts = 100L, validate = TRUE)
```

Behaviour:

1.  If `validate`, validate `model_structure`.
2.  n=1 only; dispatch on `detection_link_type`.
3.  Build a `range_df` whose row-names are exactly
    `mathscale_names(model_structure)` (minus any `mask` names). For
    link `"ltsgr"` the existing data-driven heuristics apply unchanged.
    For link `"ltsgr_iv"` we need new heuristics for `gammatil` and
    `beta` ranges and a different rule for the detection-link parameter
    set (no `ctil`); heuristic details are listed in §10.2 (open
    question).
4.  Sobol' design via `pomp::sobol_design`, identical to v1.

`start_parms` stays pure R; no C++ port and no `_r` reference.

### 5.9 `optimize_likelihood`

```         
optimize_likelihood(env_dat, occ, model_structure, mask = NULL,
                    num_starts = 100L, breadth = 1,
                    parallel = FALSE,
                    num_threads = RcppParallel::defaultNumThreads(),
                    control = list(), verbose = FALSE,
                    validate = TRUE)
```

Behaviour: drop-in replacement for v1 with `model_structure` dispatch.
Validates once at the top (when `validate = TRUE`); the multi-start
inner loop and the XPtr factory are constructed with `validate = FALSE`.

### 5.10 `profile_likelihood`

```         
profile_likelihood(env_dat, occ, model_structure, mask = NULL,
                   params_math_optim, profile_parameter,
                   increment_left = 0.1, increment_right = increment_left,
                   num_steps_left = 20L, num_steps_right = num_steps_left,
                   alpha = 0.95,
                   num_threads = RcppParallel::defaultNumThreads(),
                   control = list(), verbose = FALSE,
                   validate = TRUE)
```

The argument previously named `optim_param_vector` is renamed
`params_math_optim`. Behaviour is identical for link `"ltsgr"` n=1; link
`"ltsgr_iv"` simply uses a different set of math-scale free parameter
names (no `ctil`, but `gammatil` and `beta`).

### 5.11 `dist_between_params`

```         
dist_between_params(model_structure, p1, p2,
                    mask = NULL, give_closest_rep = FALSE,
                    validate = TRUE)
```

The dispatcher reads `detection_link_type` from `model_structure`. The
Hungarian distance treats `mu`, the widths, and `o_mat` identically
across link types. The link-specific scalar parameters (`ctil` for
`"ltsgr"`; `gammatil` and `beta` for `"ltsgr_iv"`) and `pd` contribute
additively to the squared distance the same way `ctil` and `pd` did in
v1.

### 5.12 `validate_model_structure` and `validate_params_bio`

```         
validate_model_structure(model_structure, env_dat = NULL)
validate_params_bio(model_structure, params_bio, env_dat = NULL)
```

Both validators are exported. Both **must accept arbitrary n** (n=1 and
n\>=2) and apply all relevant structural checks. The intent is that they
are written in v2 to the v3-eventual specification and do not need to
change at the v2→v3 transition.

When `env_dat = NULL`, the validators run all checks that can be
performed without env_dat. They do not emit a `message()` about skipped
checks. The calling routine is responsible for passing `env_dat`
whenever it has env_dat in hand, so that the env_dat- dependent checks
(index-in-range, `max_lag <= dim(env_dat)[2] - 3`) do run when
applicable.

The check list is essentially that of source2 §5, with the overrides
listed in §4.4 of this document:

-   `detection_link_type` is a string (`"ltsgr"` or `"ltsgr_iv"`).
-   For n=1 the n\^2+4 element (`anchored_entry`) must be **absent**.
    Presence is a validation error. (Rationale: this aligns the user's
    mental model of "anchored_entry is unused for n=1" with the inputs
    the user actually constructs.)
-   For n=1 with p=1, `params_bio[[1]]$o_mat` must be **absent**.
    Presence is a validation error.
-   For n=1 link `"ltsgr_iv"`, `params_bio[[2]]` has names `pd`,
    `gammatil`, `beta`.

A lot of the validation logic can be ported from source3; do that rather
than re-deriving it from source2 §5.

### 5.13 `Sttil_mean_lrv` (replacement for `like_neg_ltsgr`)

```         
Sttil_mean_lrv(env_dat, mu, sigltil, sigrtil, o_mat,
               max_lag = 0L,
               compute_lrv = FALSE,
               num_threads = RcppParallel::defaultNumThreads())
```

This function replaces and extends `like_neg_ltsgr`. Like
`like_neg_ltsgr`, it is **not exported**, so its arguments do not need
to be in `params_bio` form and no `model_structure` argument is needed.
Like `like_neg_ltsgr`, **input validation is omitted**; the function
trusts callers.

Returns:

-   When `compute_lrv = FALSE` (the link `"ltsgr"` use case): a numeric
    vector of length `n_loc` giving, per location, $\mean{\tilde{S}_t}$
    — the time mean of $\tilde{S}_t$ as defined in source1 eq.
    (\ref{eq:St}), evaluated with the supplied widths.
-   When `compute_lrv = TRUE` (the link `"ltsgr_iv"` use case): a list
    (or two-column matrix) with two slots:
    -   `mean_Sttil` — per-location $\mean{\tilde{S}_t}$ as above;
    -   `lrv_Sttil` — per-location $\tilde{W}$, the long-run variance of
        $\{\tilde{S}_t\}$ estimated by the Bartlett-weighted truncated
        autocovariance sum of source1 eq. (\ref{eq:iv_comp}) with
        truncation lag `max_lag`.

Notes:

-   The return convention differs from `like_neg_ltsgr`:
    `like_neg_ltsgr` returned `0.5 * mean(\tilde{S}_t)`.
    `Sttil_mean_lrv` returns `mean(\tilde{S}_t)` (no factor of 1/2).
    Callers that previously used the half-quantity scale must multiply
    by `0.5` after the call. This is a deliberate cleanup so the
    function name is meaningful (`mean_Sttil` is what is returned).
-   Both an R version (`Sttil_mean_lrv` in `R/Sttil_mean_lrv.R`) and a
    C++ version (with a thin R wrapper) must be provided, returning
    exactly the same numeric values. They are tested against each other
    for parity at tolerance 1e-6, exactly as `like_neg_ltsgr` is today.
-   `Sttil_mean_lrv` supersedes `like_neg_ltsgr`. After Stage 3b,
    `like_neg_ltsgr` is removed. See §7 for the staging.

The `compute_lrv = TRUE` branch performs a single sweep of `env_dat` to
produce both the mean and the long-run variance, avoiding the
double-pass cost that a separate iv-only kernel would incur. The
underlying lower-level kernel (currently `like_ltsg`) may need to be
refactored or supplemented to expose the per-year `\tilde{S}_t` sequence
for the variance computation; this is an implementation detail of Stage
3a (R version first, C++ version second; both verified against
`like_neg_ltsgr` for the mean part, and against hand-computed short
series for the variance part).

### 5.14 `habitat_suitability` / `vsp`

```         
habitat_suitability(model_structure, params_bio, env_list,
                    output = "", overwrite = FALSE,
                    return_prob = TRUE, threads = 0L, wopt = list(),
                    validate = TRUE)

vsp(env_data, model_structure, params_bio,
    return_raster = FALSE, validate = TRUE)
```

The hot tile loop in `habitat_suitability` builds a per-tile `env_dat`
chunk and forwards `model_structure` and `params_bio` (with
`validate = FALSE`) to the `log_prob_detect` C++ kernel updated in Stage
5. `vsp` stays a thin wrapper around `log_prob_detect`.

### 5.15 `interpret_parameters`

```         
interpret_parameters(model_structure, params_bio, plot_indices,
                     plot_lims = NULL, env_dat = NULL, occ = NULL,
                     breadth = 1, validate = TRUE, ...)
```

The diagnostic plots interpret the *growth-environment* function, which
does not depend on the choice of detection link, so the inferred contour
shape is computed identically for both links. The plotting code reads
`params_bio[[1]]` for `mu`, `sigltil`, `sigrtil`, `o_mat`. The
link-specific parameters are not used for the plot.

------------------------------------------------------------------------

## 6. The "no-deltas" verification component

The invariant: for every legal v1 call, the equivalent v2 call (using
n=1 link `"ltsgr"`) returns the same value, up to documented tolerance.
The mechanism:

A single directory `xsdm-devel/plan/nodeltas/`. For each (A)-class
function `<func>`, there is a sub-folder:

```         
xsdm-devel/plan/nodeltas/<func>/
  ├── pre.R          # run on v1
  ├── post.R         # run on v2
  ├── snapshot.rds   # produced by pre.R, consumed by post.R
  └── notes.md       # human notes about case coverage and tolerance
```

-   `pre.R`:
    -   Sets a fixed seed.
    -   Calls a generator that produces a list of test cases (each a
        list of v1-style arguments). Cases must cover the small-n
        regime, the realistic regime (size of `example_1$env_array`,
        p=2), boundary models (Inf masks where applicable), the
        non-default `mask` case where applicable, and the non-default
        flag combinations.
    -   For each case, calls v1 `<func>` and stores the result.
    -   Calls a converter `to_v2_args(<v1-args>)` (in
        `plan/nodeltas/v1_to_v2.R`) to produce the equivalent v2
        arguments.
    -   Writes one rds bundle:
        `r     saveRDS(list(cases_v1 = ..., cases_v2 = ..., results_v1 = ...,                  metadata = list(R_version = ..., xsdm_version = "1.0.0",                                  generated_at = Sys.time())),             file = "snapshot.rds")`
-   `post.R`:
    -   Loads `snapshot.rds`, calls v2 `<func>` with `cases_v2[[i]]`,
        compares to `results_v1[[i]]` with
        `all.equal(..., tolerance =     1e-8)`, and reports PASS/FAIL.

This is one snapshot file per function. The snapshots are committed to
`version2` as part of each stage; they are versioned data and are not
regenerated after the initial commit unless the v1 reference itself
changes. The whole `plan/nodeltas/` tree is removed (or kept as
historical material under `plan/`, the user's call) once v2 ships; **no
permanent v1↔v2 regression test is retained**.

The shared converter file `plan/nodeltas/v1_to_v2.R` contains pure
functions: `to_v2_args_loglik_bio`, `to_v2_args_log_prob_detect`, etc.
They take v1-style argument lists and emit v2-style argument lists;
example for `loglik_bio`:

``` r
to_v2_args_loglik_bio <- function(v1_args) {
  p <- length(v1_args$mu)
  list(
    env_dat = v1_args$env_dat,
    occ = v1_args$occ,
    model_structure = list(
      seq_len(p),         # element 1: matrix-entry indices
      "Lyapunov",         # element 2: ltsgr_iv_method
      0L,                 # element 3: max_lag (unused for "ltsgr" n=1)
      "ltsgr"             # element 4: detection_link_type
    ),
    params_bio = list(
      c(list(mu = v1_args$mu, sigltil = v1_args$sigltil,
             sigrtil = v1_args$sigrtil),
        if (p > 1) list(o_mat = v1_args$o_mat) else list()),
      list(pd = v1_args$pd, ctil = v1_args$ctil)
    ),
    return_prob = v1_args$return_prob,
    sum_log_p = v1_args$sum_log_p,
    num_threads = v1_args$num_threads
  )
}
```

The converters are tested implicitly: if they're wrong, `post.R` fails.
The `_r_vs_cpp` parity tests live in `tests/testthat/` and are unrelated
to no-deltas.

------------------------------------------------------------------------

## 7. Staged plan

Stages are mostly ordered by dependency. Stage v1px sits at the very
beginning; the rest are renumbered relative to plan1 to make room.

### Stage v1px: pre-migration cleanup

**What.** A list of small cleanups to v1 code that should happen before
the v2 refactor begins. Some of these may be done on `main` itself and
shipped as part of the initial v1 CRAN release; others may land on
`version2` first. The decision is item-by-item.

Items in this stage:

#### 7.1.a Argument and identifier consistency

Audit identifier casing/spelling across the codebase and standardise
inconsistent forms. Known issues:

-   `parms` vs `params`: e.g. `bio_to_math` takes `parms_bio`, while
    every other place uses `params_*`. Standardise on `params_bio`,
    `params_math` everywhere.
-   `param_list` (in v1 `vsp`, `habitat_suitability`,
    `interpret_parameters`) vs the proposed `params_bio`. These all
    become `params_bio` in v2; for the cleanup stage, simply rename
    `param_list` → `params_bio` consistently in v1 (keeping the v1
    flat-list shape), so the v2 refactor doesn't touch the identifier in
    addition to the shape.
-   `numstarts` vs `num_starts` (in `get_start_parms_`). Standardise on
    `num_starts`.
-   `param_vector` vs `params_math` in v1 `loglik_math` / `math_to_bio`:
    rename `param_vector` → `params_math`, again purely cosmetic in v1.

A short script `plan/cleanup/v1px_renames.R` enumerates the renames and
the affected files, so the diff is mechanical.

#### 7.1.b `bio_to_math` documentation and error handling

The current `bio_to_math` has docs that key on whether
`det(o_mat) == -1`, suggesting that's the only problematic case. That's
not strictly true: `expm::logm` may also fail on other
special-orthogonal inputs that lie outside the principal-logarithm
domain. The cleanup is:

-   Simplify the docstring to state plainly that `bio_to_math` uses
    `expm::logm` to compute `o_par` from `o_mat`, and that any error
    from `expm::logm` is propagated to the caller unchanged.
-   Remove the special-case error messages for `det == -1`. Just let
    `expm::logm`'s error propagate.
-   Optionally: discuss whether to drop `bio_to_math` entirely. The
    decision is to **keep** it for completeness (it is the inverse of
    `math_to_bio` and is useful in tests and example construction), but
    pare it down to the minimum.

#### 7.1.c Other cleanup candidates (Claude-found)

The following minor issues have been observed in v1 source and are
candidates for cleanup. They are listed here for triage; subsequent
revisions of this plan should keep some, drop others, defer some to v2
stages.

1.  **`profile_likelihood`** has a commented-out mask validation block
    (around lines 246–253). Decide whether to delete the comment or
    fix-and-uncomment.
2.  **`vsp`** has two consecutive
    `if (!requireNamespace("terra",     quietly = TRUE))` checks; the
    second is dead code. Remove.
3.  **`dist_between_params`** has the opaque
    `checkmate::assert_integerish(sqrt(9 + 8 * (length(p1) +     length(mask))))`
    validation. Replace with a clearer `length`-based check that
    documents what it's actually checking.
4.  **`R/internals.R`** is large and mixes pure-R reference
    implementations with general utilities (logit, permutations,
    check_env_array). Consider splitting into `R/internals_references.R`
    and `R/internals_utils.R` for readability.
5.  **`tests/manual/`** scripts: review which are still relevant and
    which are stale.
6.  **`start_parms` / `get_range_df_quant_vec` deprecation**: the
    `quant_vec` argument has been soft-deprecated since v0.3. v1 can
    hard-remove it, since v1 itself is the first CRAN release and there
    are no external users to break.
7.  **Docstring style consistency**: `@returns` vs `@return`, casing of
    "Logical" vs "logical", the order of `@examples` blocks. Mostly
    cosmetic. Run `roxygen2` after; ensure NAMESPACE is deterministic.
8.  **Examples that reference packed rasters** (`bio01`, `bio12`) use
    `terra::unwrap` inline; ensure consistent `\donttest` / `\dontrun`
    usage.
9.  **`test-numpars.R`** mixes tests for `num_par` and `num_env_var`;
    consider splitting or renaming for clarity.
10. **`like_ltsg`** documentation is currently very minimal (it's
    flagged `@export` but documented like an internal). Either flesh out
    the docs or unexport it (it is a low-level building block and most
    users will never call it directly).
11. **`make_loglik_math_xptr` (the older R-callback XPtr factory in
    `src/loglik_math_xptr.cpp`)** is unused except by its own parity
    test (`tests/testthat/test-loglik_math_xptr_parity.R`), and the
    comment in `R/profile_likelihood.R:69` documents a known bug that
    makes it unsuitable for the production hot path. Delete the `.cpp`,
    the parity test, and the corresponding `RcppExports` entry; the
    pure-C++ `make_loglik_math_xptr_cpp` factory remains the canonical
    XPtr factory.

#### 7.1.d Other cleanups that may emerge

This subsection is a placeholder. Subsequent revisions of this plan will
finalize the list. Some items may slide into v2 stages if it's cleaner
to do them alongside the function-by-function refactor.

**Tests / checkpoints.** Each cleanup item lands as its own commit with
the existing test suite passing. No new tests are required; existing
tests must not regress.

**Branch / commit / merge.** Sub-branch `version2-stagev1px-<tag>` for
each substantial item; merge back to `version2`. Items that should also
reach `main` (the v1 CRAN release) are cherry-picked or re-applied on a
separate `main`-rooted branch and PR'd to `main` independently.

### Stage 0: planning artefacts

Land plan2.md (this file) in `plan/`. No code changes. Single commit on
`version2`.

### Stage 1: validators

Add `validate_model_structure(model_structure, env_dat = NULL)` and
`validate_params_bio(model_structure, params_bio, env_dat = NULL)` in
`R/validate_model_structure.R` and `R/validate_params_bio.R`. **Both
must handle arbitrary n.** Port logic from source3 where appropriate.

Tests:

-   `tests/testthat/test-validate_model_structure.R` covers the n=1
    happy paths (both link types), the n\>=2 happy paths (a couple of
    representative shapes), and one error per check listed in §4.4 and
    source2 §5.1.
-   `tests/testthat/test-validate_params_bio.R` likewise.

No-deltas: not applicable.

Branch: `version2-stage1-validators`; merge back to `version2`.

### Stage 2: `mathscale_names` and the canonical-names plumbing

Add `mathscale_names(model_structure)` in `R/mathscale_names.R`. **Must
work for any n \>= 1.** Implement n=1 link `"ltsgr"` / link `"ltsgr_iv"`
cases and the general n\>=2 case (porting from source3 where
appropriate). Re-implement `make_mask_names`, `num_par`, `num_env_var`
as thin wrappers around `mathscale_names`. They remain exported
(unchanged NAMESPACE) so the v1 export surface is undisturbed between
Stage 2 and Stage 7; this avoids touching call sites twice.
`create_mask` and `create_param_vector_masked` continue to work via the
aliased `make_mask_names` and the existing C++ kernel. Stage 7 deletes
all five at once.

Tests:

-   `tests/testthat/test-mathscale_names.R` covers n=1 link `"ltsgr"`
    (including agreement with the existing `make_mask_names(p)`), n=1
    link `"ltsgr_iv"`, and a few n\>=2 cases.
-   Existing tests for `make_mask_names`, `num_par`, `num_env_var`
    continue to pass through the alias.

No-deltas: trivial; covered by the existing tests.

Branch: `version2-stage2-mathscale-names`.

### Stage 3a: introduce `Sttil_mean_lrv` alongside `like_neg_ltsgr`

Add `Sttil_mean_lrv` (R version + C++ version with R wrapper) per §5.13,
**alongside** the existing `like_neg_ltsgr`. Do not yet remove or modify
`like_neg_ltsgr`.

Tests:

-   `tests/testthat/test-Sttil_mean_lrv.R` (new) covers the
    `compute_lrv = FALSE` path (mean of `S_t`) by hand-computed short
    series and (for the same widths) by parity with
    `2 * like_neg_ltsgr(...)` (since `like_neg_ltsgr` returns half the
    mean).
-   Same file covers the `compute_lrv = TRUE` path (`mean_Sttil` and `lrv_Sttil`)
    with hand-computed short series at several values of `max_lag`.
-   R-vs-C++ parity at tolerance 1e-6.

No-deltas: not applicable (new function).

Branch: `version2-stage3a-Sttil_mean_lrv`.

### Stage 3b: replace `like_neg_ltsgr`

Replace every internal call site of `like_neg_ltsgr` with
`Sttil_mean_lrv`, adjusting the factor of 0.5 at the call site where
necessary. Remove `like_neg_ltsgr` from `R/`, `src/`, `NAMESPACE`,
documentation, and tests. The corresponding pure-R reference
`like_neg_ltsgr_r` is also removed.

Tests:

-   `test-like_neg_ltsgr*.R` files are removed.
-   The full test suite continues to pass on `Sttil_mean_lrv`.

No-deltas: per §6, the no-deltas snapshots for `log_prob_detect`,
`loglik_bio`, etc. (which transitively depend on `like_neg_ltsgr`) are
produced *now*, before stage 5 and beyond, so they will naturally
exercise this replacement. Stage 3b itself doesn't need a dedicated
no-deltas check.

Branch: `version2-stage3b-remove-like_neg_ltsgr`.

### Stage 4: `math_to_bio` / `bio_to_math` reshape

Change `math_to_bio` to take `model_structure, params_math` and return a
`params_bio` list. `bio_to_math` likewise (input is now `params_bio`,
output unchanged in shape). Update `.math_to_bio_cpp` and
`.build_canonical_param_vector_cpp`. Update R references. Add
`bio_to_math_r` to `R/internals.R`.

`bio_to_math`'s simplified docstring (per §7.1.b) lands in this stage if
it didn't land during Stage v1px.

No-deltas: pre-run for `math_to_bio` and `bio_to_math` before this stage
starts. The shape difference is handled in the converter and assertion:
the snapshot stores the v1-shape output of `math_to_bio`, and `post.R`
calls a small helper `flatten_params_bio_to_v1_list(v2_out)` (in
`plan/nodeltas/v1_to_v2.R`) before comparing.

Tests:

-   `test-math_to_bio.R`, `test-math_to_bio_cpp.R`,
    `test-math_to_bio_r_vs_cpp.R` updated to the new signature.
-   `test-bio_to_math.R` updated; new parity test
    `test-bio_to_math_r_vs_cpp.R` added.
-   `plan/nodeltas/math_to_bio/post.R` PASS.
-   `plan/nodeltas/bio_to_math/post.R` PASS.

Branch: `version2-stage4-math-bio`.

### Stage 5: `log_prob_detect`

Change the signature; add the `"ltsgr_iv"` branch (which calls
`Sttil_mean_lrv(..., compute_lrv = TRUE)` and applies the closed
formula); update R reference; update C++ kernel. The link-type dispatch
is the only logic addition relative to v1 — both branches are simple
closed-form expressions over `Sttil_mean_lrv`'s outputs, exactly the
spirit in which v1 implements the link `"ltsgr"` formula.

No-deltas: pre-run before this stage. Snapshot covers the link-`"ltsgr"`
n=1 case at p=1 and p=2, several seeds, both `return_prob` flags.

Tests:

-   `test-log_prob_detect.R`, `test-log_prob_detect_cpp.R`,
    `test-log_prob_detect_r_vs_cpp.R` updated. **Per AP2, link
    `"ltsgr_iv"` tests live in the same files as the link `"ltsgr"`
    tests, in clearly-marked sections; no separate `*_link1.R` files.**
-   `plan/nodeltas/log_prob_detect/post.R` PASS.

Branch: `version2-stage5-log-prob-detect`.

### Stage 6: `loglik_bio`

Change the signature; thread the link dispatch through; update R
reference; update C++ kernel.

No-deltas: pre-run before this stage. Snapshot covers all four flag
combinations of `(return_prob, sum_log_p)` and both p=1 and p=2.

Tests:

-   `test-loglik_bio.R`, `test-loglik_bio_cpp.R`,
    `test-loglik_bio_r_vs_cpp.R` updated; link `"ltsgr_iv"` cases in the
    same files.
-   `plan/nodeltas/loglik_bio/post.R` PASS.

Branch: `version2-stage6-loglik-bio`.

### Stage 7: `loglik_math`, math-scale helper consolidation, XPtr factory

This stage does two things at once: (i) the `model_structure`-aware
refactor of `loglik_math` and the XPtr factory, and (ii) the full
consolidation of the math-scale-helper layer described in §5.6. After
this stage, the v1 helpers `make_mask_names`, `num_par`, `num_env_var`,
`create_mask`, and `create_param_vector_masked` no longer exist in the
codebase.

Work items:

-   Refactor `loglik_math` to take `model_structure` (R wrapper, R
    reference, and C++ kernel `loglik_math_cpp`). The XPtr factory
    `make_loglik_math_xptr_cpp` is updated analogously: it captures
    `model_structure` along with the data, mask, and threading config,
    and constructs the closure with `validate = FALSE` semantics for the
    inner-loop body.
-   Introduce a single internal helper
    `assemble_params_math_(model_structure, params_math, mask)` in
    `R/internals.R`. It validates that the names of `params_math` and
    `mask` together equal `mathscale_names(model_structure)` with no
    overlap and no missing slots, then delegates the merge to the C++
    kernel `.build_canonical_param_vector_cpp` (whose signature is
    updated to take `model_structure` instead of `p`). This helper
    contains all the validation logic that previously lived in
    `create_param_vector_masked`; it is the only place that logic lives.
-   Update every internal call site (`optimize_likelihood`,
    `loglik_math_r`, `dist_between_params`, `dist_between_params_r`)
    that previously called `create_param_vector_masked` to call
    `assemble_params_math_()` instead.
-   Update every internal call site that used
    `names(make_mask_names(p))` to use
    `mathscale_names(model_structure)` directly. This affects
    `R/dist_between_params.R`, `R/internals.R` (math_to_bio_r and
    dist_between_params_r), `R/loglik_math.R`, `R/math_to_bio.R`, and
    `R/bio_to_math.R`. The one `bio_to_math.R` site that uses
    `make_mask_names(p)` as a NA-vector scaffold is rewritten to either
    call `setNames(rep(NA_real_,     length(nm)), nm)` once or, more
    cleanly, to assemble the math-scale vector slot-by-slot without a NA
    scaffold (since all slots are overwritten anyway).
-   Delete `R/make_mask_names.R`, `R/num_par.R`, `R/num_env_var.R`,
    `R/create_mask.R`, `R/create_param_vector_masked.R`. Delete the
    corresponding `tests/testthat/test-make_mask_names.R`,
    `test-numpars.R`, `test-create_mask.R`,
    `test-create_param_vector_masked.R`,
    `test-create_param_vector_masked_r_vs_cpp.R`. (Their behaviour is
    now exercised through the public-API tests of `loglik_math`,
    `optimize_likelihood`, and `dist_between_params`, plus a small new
    `test-assemble_params_math.R` covering the validation paths.)
-   Update NAMESPACE: remove the five exports.
-   Update any roxygen `[make_mask_names()]` / `[num_par()]` /
    `[create_mask()]` / `[create_param_vector_masked()]` references in
    other docstrings to point at `mathscale_names()` (or remove the
    cross-reference if it no longer makes sense).

No-deltas: pre-run for `loglik_math`. Snapshot covers `mask = NULL`,
mask with one fixed `mu`, mask with `Inf` for `sigl`/`sigr`/`pd`. No
no-deltas snapshots are produced for the deleted helpers; their job is
now subsumed inside `loglik_math` and `optimize_likelihood`, whose own
snapshots transitively exercise the new helper.

Tests:

-   The corresponding `test-*` files for `loglik_math` and the XPtr
    factory updated.
-   New `tests/testthat/test-assemble_params_math.R` covers the
    validation paths (overlap rejection, missing-slot rejection, Inf in
    non-permitted slots, etc.) directly on the internal helper via
    `xsdm:::assemble_params_math_`.
-   `plan/nodeltas/loglik_math/post.R` PASS.
-   `devtools::check()` is expected to pass at the end of this stage.

Branch: `version2-stage7-loglik-math`.

### Stage 8: `start_parms`, `optimize_likelihood`

Update both for the new arguments. The `start_parms` heuristic for link
`"ltsgr_iv"` is new — see §10.2 for the open question on what ranges to
use for `gammatil` and `beta`.

No-deltas: pre-run. `start_parms` snapshots are deterministic given the
seed used inside `pomp::sobol_design`. `optimize_likelihood` is
non-deterministic in the sense that floating-point noise across
optimizer iterations can change the last digit of the optimum; the
no-deltas check therefore uses `tolerance = 1e-4` for `loglik` and
`1e-2` for the parameter-vector components, with a documented caveat in
`plan/nodeltas/<func>/notes.md`.

Tests:

-   `test-start_parms.R`, `test-optimize_likelihood.R`,
    `test-get_start_parms.R`, `test-get_range_df.R`,
    `test-num_starts_validation.R`, `test-optimize_helpers.R` updated.
-   `plan/nodeltas/start_parms/post.R` PASS.
-   `plan/nodeltas/optimize_likelihood/post.R` PASS (relaxed tolerance).

Branch: `version2-stage8-starts-opt`.

### Stage 9: `profile_likelihood`

Update for the new arguments. Internal logic is unchanged.

No-deltas: pre-run with one small case (1 step left, 1 step right) and
one realistic-size case. Same loose tolerances as Stage 8.

Tests:

-   `test-profile_likelihood.R` updated.
-   `plan/nodeltas/profile_likelihood/post.R` PASS.

Branch: `version2-stage9-profile`.

### Stage 10: `dist_between_params`

Update for the new arguments. The link-`"ltsgr_iv"` branch is small
(different scalar parameters).

No-deltas: pre-run. Tolerance 1e-8 (deterministic).

Tests:

-   `test-dist_between_params.R`, `test-solve_lsap_cpp.R` updated.
-   `plan/nodeltas/dist_between_params/post.R` PASS.

Branch: `version2-stage10-dist`.

### Stage 11: `habitat_suitability` and `vsp`

Update both for the new arguments. The hot-tile loop in
`habitat_suitability` forwards `model_structure` and `params_bio` (with
`validate = FALSE`) to the C++ kernel updated in Stage 5.

No-deltas: pre-run with one small raster case to keep the snapshot
small.

Tests:

-   `test-habitat_suitability.R`, `test-vsp.R` updated.
-   `plan/nodeltas/habitat_suitability/post.R` PASS.
-   `plan/nodeltas/vsp/post.R` PASS.

Branch: `version2-stage11-rasters`.

### Stage 12: `interpret_parameters`

Update signature. Internal logic unchanged in shape; reads fields out of
`params_bio[[1]]` instead of a flat list.

No-deltas: function is mostly side-effects, so the snapshot stores only
the numeric `z` matrix (or 1D `y` vector) computed inside `compute_y`
(refactor `compute_y` into a small internal `compute_contour_grid_` so
it is callable without drawing).

Tests:

-   `test-interpret_parameters.R` updated.
-   `plan/nodeltas/interpret_parameters/post.R` PASS.

Branch: `version2-stage12-interpret`.

### Stage 13: examples, vignettes, docs, NEWS

Update every example block (most use `example_1$par_vec` or similar and
must be reshaped). Update `vignettes/`. Bump `DESCRIPTION` to
`Version: 2.0.0`. Add NEWS.md entry summarising the breaking-API change
and the new detection link.

Tests: `devtools::check(vignettes = TRUE)` passes. `devtools::test()`
passes.

Branch: `version2-stage13-docs`.

### Stage 14: performance regression

A dedicated stage to address performance concerns identified during
implementation. Items to look at, all of which are deferred from the
function-by-function stages:

-   Validators called inside the optimizer hot path: confirm via
    profiling that the `validate = FALSE` propagation (§5.2) brings the
    per-loglik validation cost to negligible.
-   The `Sttil_mean_lrv(compute_lrv = TRUE)` path: profile a realistic
    optimization (≈100 starts × ≈2000 iterations each on
    `example_1`-size data) and confirm the iv computation is not a
    bottleneck. For `max_lag = 4` and `ts_length = 39` and
    `n_loc = 2000`, the iv cost is \~300k extra ops per loglik; this
    should be negligible compared to the matrix sweep already done in
    the link-`"ltsgr"` path, but verify.
-   The XPtr factory in `loglik_math_xptr_cpp` captures more state (the
    parsed `model_structure`); confirm no per-call R round-trip.
-   General optimizer-hot-path benchmarking on a representative fitting
    workload, with link `"ltsgr"` for parity to v1 and link `"ltsgr_iv"`
    to characterise the new path.

This stage produces a short report `plan/perf/stage15_report.md` and any
code adjustments needed to hit reasonable performance targets. Targets
are not numeric in this plan; they will be set after the first benchmark
run.

Branch: `version2-stage14-perf`.

### Stage 15: integration regression and final review

From the now-fully-updated `version2`, run a single final regression
that calls every public function in v2 with the v1 test data and a
reasonable v2 wrapping. This is essentially the union of all
`plan/nodeltas/*/post.R` runs, executed end-to-end via a
`plan/nodeltas/run_all.R` driver. Document any deltas that exceed
documented tolerances. Produce `plan/nodeltas/final_report.md`.

Branch: `version2-stage15-final`.

### Dependency map

```         
v1px ─→ 0 ─→ 1 ─→ 2 ─→ 3a ─→ 3b ─→ 4 ─→ 5 ─→ 6 ─→ 7
                                                      ↓
              ┌───────────────────────────────────────┤
              ↓                                        ↓
              8                                        9
              ↓                                        ↓
              └─→ 10 ─→ 11 ─→ 12 ─→ 13 ─→ 14 ─→ 15
```

Stage 1 (validators) and Stage 3a (`Sttil_mean_lrv`) could in principle
be done in parallel since neither depends on the other. We do them
sequentially to keep the linear narrative simple, but a fast-track
parallel approach is acceptable.

------------------------------------------------------------------------

## 8. Branch / commit / merge conventions

-   **Branch.** `version2` is the trunk for all v2 work (already
    exists). Each stage gets a sub-branch named
    `version2-stage<N>-<short-tag>`, branched off `version2`.
    Sub-branches are short-lived. Force-pushes are allowed on
    sub-branches up until the merge; never on `version2` and never on
    `main`.
-   **Commits.** Small and atomic, with informative subject lines
    starting with the stage tag (e.g.
    `stage5: log_prob_detect     signature & R wrapper`,
    `stage5: log_prob_detect link "ltsgr_iv"     C++ branch`,
    `stage5: log_prob_detect parity tests for "ltsgr_iv"`).
-   **Merging back to `version2`:** standard `git merge --no-ff` so each
    stage shows up as a merge commit.
-   **`main`:** items from Stage v1px that should reach the v1 CRAN
    release are cherry-picked or re-applied on a `main`-rooted branch
    and PR'd to `main` independently. Otherwise `main` is untouched
    until v2 is ready to release. Final merge to `main` is a single PR
    that goes through `devtools::check()`, `plan/nodeltas/run_all.R`,
    and final review.
-   **Tags.** Tag `v1.0.0` on `main` before any v2 work merges back (if
    not already tagged). Tag `v2.0.0` on the v2-merge commit.
-   **No-deltas snapshots are committed** to `version2` as part of each
    stage; they are versioned data and not regenerated after the initial
    commit.

------------------------------------------------------------------------

## 9. Testing summary

-   **Unit tests** (`tests/testthat/`): every (A)-class function has its
    v2 unit tests; every (B)-class function has new unit tests; every
    `*_r_vs_cpp` parity test is updated. Per AP2, link `"ltsgr_iv"`
    tests live in the same `test-<func>.R` file as the link `"ltsgr"`
    tests, in clearly-labelled sections — there are no separate
    `*_link1.R` files.
-   **No-deltas tests** (`plan/nodeltas/*/post.R`): one per (A)-class
    function; all pass at their documented tolerances.
-   **Vignette compilation** (`devtools::check(vignettes = TRUE)`):
    green from Stage 13 onward.
-   **R CMD check** (`devtools::check()`): green from Stage 13 onward
    (Stage 13 is where examples and vignettes are reshaped to match
    every (A)-class function's new signature; before that, mid-stage
    documentation references can still trip the check).
-   **Manual visual inspection**: only `interpret_parameters` requires
    this. A short visual-comparison script lives in
    `tests/manual/interpret_parameters_v1_v2_visual.R`.
-   **Performance benchmarks**: a small benchmarking harness lives under
    `plan/perf/`; targets are set during Stage 14.

------------------------------------------------------------------------

## 10. Open questions and things to keep developing

These items still need a decision before the corresponding stage begins.

### 10.1 (resolved) Link `"ltsgr_iv"` parameter naming

Resolved: the bio-scale and math-scale names are `pd`, `gammatil`,
`beta`. On the math scale, `pd` is on the logit scale, `gammatil` and
`beta` are on the log scale. The names are identical across scales (no
`_log` / `_logit` suffixes). This convention is wired into §4.3 and used
throughout the plan.

### 10.2 (open) `start_parms` ranges for link `"ltsgr_iv"`

Data-driven heuristics are needed for the `gammatil` and `beta` ranges.
A defensible default:

-   `beta` math-scale range `log(c(1e-2, 1, 1e2)) = c(-4.6, 0, 4.6)`
    (centred at `beta = 1`, breadth ±2 orders of magnitude).
-   `gammatil` math-scale range derived from the empirical scale of
    `mean_Sttil` and `lrv_Sttil` at the centred starting point: pick
    `gammatil` such that
    `4 * gammatil / lrv_Sttil ≈ mean_Sttil / lrv_Sttil`, then take ±
    a factor.

These are guesses; before Stage 8 we should run a sensitivity analysis
on a virtual species to confirm reasonable optimization behaviour.

### 10.3 (resolved) Permanent v1↔v2 regression test

Resolved: no permanent snapshot is retained after v2 ships. The
no-deltas mechanism in `plan/nodeltas/` is sufficient during the
migration; once v2 ships, the snapshots are removed (or kept as
historical material, the user's call).

### 10.4 (resolved) Validator behaviour when `env_dat` is missing

Resolved: validators take `env_dat` with a default of `NULL`. When
`NULL`, all checks possible without env_dat are run; no message is
emitted. When given, the additional env_dat-dependent checks are also
run. Calling routines pass env_dat whenever they have it.

Related: for n=1, the n\^2+4 element (`anchored_entry`) must be absent.
Presence at n=1 is an error, not a silent skip. (Wired into §4.4 and
§5.12.)

Related: validators must work for arbitrary n (not just n=1) in v2. The
intent is that the validators are written in v2 to the v3-eventual
specification and don't need to change at the v2→v3 transition. Source3
contains code that can serve as a starting point. (Wired into §5.12.)

### 10.5 (resolved) Detection-link indicator: name and value

Resolved: the n\^2+3 element of `model_structure` is named
`detection_link_type` (when names are present) and its value is the
string `"ltsgr"` or `"ltsgr_iv"`. The n\^2+1 element is named
`ltsgr_iv_method` (note: with underscore between `ltsgr` and `iv`).

### 10.6 (resolved) `validate = TRUE` flag

Resolved: every (A)-class function takes `validate = TRUE` (default).
Functions that call other (A)-class functions internally do so with
`validate = FALSE` after their own validation. Wired into §5.2 and
referenced in every (A)-class signature in §5.

### 10.7 (resolved) `like_neg_ltsgr` replacement

Resolved: replaced by `Sttil_mean_lrv`. Stage 3a introduces
`Sttil_mean_lrv` alongside `like_neg_ltsgr`; Stage 3b removes
`like_neg_ltsgr`. There is no separate `iv_estimator` or
`detection_link_iv` function in v2; the iv detection probability is
applied inline inside `log_prob_detect` (using `Sttil_mean_lrv`'s two
outputs). Wired into §5.13 and §7.

### 10.8 (open) `o_mat` slot in n=1 `params_bio[[1]]` when p=1

Resolved: `o_mat` must be **absent** in `params_bio[[1]]` when p=1. The
validator rejects a list that contains an `o_mat` slot at p=1. (Wired
into §4.2 and §5.12.)

### 10.9 (resolved) Boundary models for link `"ltsgr_iv"`

Resolved: no new boundary model is introduced for link `"ltsgr_iv"`. The
only boundary models are the existing ones: one or more of
`sigltil*`/`sigrtil*` equal to `Inf` on the bio scale (equivalently
`Inf` on the math scale via `log`), and `pd = 1` on the bio scale
(equivalently `Inf` on the math scale via `logit`). For `gammatil` and
`beta`, `Inf` is rejected in `mask` validation.

### 10.10 (open) Final list of Stage v1px cleanup items

The list in §7.1.c is a triage candidate set. Future revisions of this
plan finalize which items go to v1 (cherry-pick to `main`) versus to
`version2` only versus deferred to a v2 stage versus dropped.

### 10.11 (deferred) Should `model_structure` be an S3 class?

Promote `model_structure` to an `S3` class `xsdm_model_structure` with
constructor and `print` method? Postponed to a future version; not
blocking for v2.0. Listed in §11.

### 10.12 (open) C++ implementation of `Sttil_mean_lrv` `compute_lrv = TRUE`

Whether to implement the `compute_lrv = TRUE` path in C++ as a single
sweep that produces both `mean_Sttil` and `lrv_Sttil` (preferred), or as two
sweeps reusing the existing `like_ltsg` kernel for the mean and a
separate kernel for the variance. The former is recommended; the latter
is acceptable as a fallback. Decide during Stage 3a based on what the
lower-level kernel can cleanly expose.

------------------------------------------------------------------------

## 11. Ideas for possible future versions

Items deferred from the v2 plan that are worth revisiting later:

-   **Promote `model_structure` to an S3 class**
    (`xsdm_model_structure`) with a constructor, a `print` method, and
    possibly subset accessors. This would catch shape errors at
    construction time and make `print(model_structure)` informative.
    Cost: a small style decision that ripples through every function.
    Defer to v2.x or v3.
-   **AIC/BIC comparison across detection links and other structural
    differences**. With v2's `model_structure`, fitting with `"ltsgr"`
    and again with `"ltsgr_iv"` and comparing by AIC/BIC is
    straightforward at the user level. A higher-level helper that
    automates this comparison (and similar comparisons across other
    structural choices once n\>=2 lands) is a candidate for v2.x or
    later.
-   **Permanent compact v1↔v2 regression** (e.g. one tiny snapshot on
    `example_1`). Considered and rejected for v2 (see §10.3); may be
    revisited if practice shows the no-deltas tolerance gets
    accidentally relaxed across patches.
-   **Hard-deprecate `make_mask_names`, `num_par`, `num_env_var`,
    `create_mask`, `create_param_vector_masked`** fully (they are
    deleted in v2 Stage 7 — both the exports and the internal aliases —
    so nothing further is required at the v2 release; listed here for
    historical completeness).
-   **Refactor `like_ltsg` to expose the per-year `S_t` grid** as a
    first-class output of a lower-level kernel; this is an option for
    the Stage 3a implementation of `Sttil_mean_lrv` C++ path, and may be
    revisited later if the v2 implementation chose the fallback
    two-sweep approach.

------------------------------------------------------------------------

## 12. Summary checklist

A condensed tick-list of what "v2 is done" looks like:

-   [ ] Stage v1px cleanups landed (some on `main`, some on `version2`).
-   [ ] All exported v1 functions have v2 signatures per §5.
-   [ ] No exported v1 signature survives.
-   [ ] `validate_model_structure`, `validate_params_bio`, and
    `mathscale_names` are exported and **work for arbitrary n**.
-   [ ] Every (A)-class function takes `validate = TRUE` (default) and
    propagates `validate = FALSE` to inner (A)-class calls.
-   [ ] `Sttil_mean_lrv` replaces `like_neg_ltsgr` (Stages 3a, 3b).
-   [ ] `make_mask_names`, `num_par`, `num_env_var`, `create_mask`,
    `create_param_vector_masked` are all removed (Stage 7); the only
    remaining canonical-names entry point is `mathscale_names`, and the
    only full-vector-assembly entry point is the internal
    `assemble_params_math_()`.
-   [ ] R reference + C++ kernel duality preserved for every hot-path
    function.
-   [ ] All `*_r_vs_cpp` parity tests pass at tolerance 1e-6.
-   [ ] Per-link tests live in the same `test-<func>.R` file (no
    `*_link1.R` files).
-   [ ] All `plan/nodeltas/*/post.R` scripts pass at their documented
    tolerances.
-   [ ] Stage 14 performance benchmarking complete; no obvious
    regressions vs v1 on link `"ltsgr"` workloads.
-   [ ] `devtools::check()` passes.
-   [ ] `devtools::check(vignettes = TRUE)` passes.
-   [ ] DESCRIPTION bumped to `Version: 2.0.0`.
-   [ ] NEWS.md describes the breaking change and the new detection
    link.
-   [ ] Vignettes and examples reshaped to use `model_structure` /
    `params_bio`.
-   [ ] All Stage-1-through-Stage-15 sub-branches merged back to
    `version2`.
-   [ ] `version2` merged to `main` in a single reviewed PR.
-   [ ] Tag `v2.0.0` lands on `main`.

------------------------------------------------------------------------

*End of plan2.md.*
