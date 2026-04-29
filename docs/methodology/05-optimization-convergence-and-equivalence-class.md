# 05 — Optimization convergence + equivalence-class normalization

## Why

A user reported that `optimize_likelihood()`'s multi-start fits "reach an
optimum that is far for the known values of parameters" and posted a
table where rows 2–10 of `solutions$full_par` looked like **eight
different parameter vectors**:

```
     mu1   mu2  sigltil1  sigltil2  sigrtil1  sigrtil2   ctil    pd  o_par1
 1  -5.79  4.87 -0.702    12.20    -2.52      5.77    -18.6  1.35    1.67
 2  14.2   6.72 -2.30      0.299   -0.737     6.90    -18.3  2.14   -7.66
 3  14.2   6.70  0.303    -0.739    0.291    -2.30    -18.3  2.14    6.48
 4  14.2   6.70 -0.739     0.291   -2.30      0.303   -18.3  2.14    1.77
 5  14.2   6.70  0.291    -2.30     0.303    -0.739   -18.3  2.14    9.62
 6  14.2   6.70 -2.30      0.303   -0.739     0.291   -18.3  2.14   -7.66
 7  14.2   6.70 -0.739     0.291   -2.30      0.303   -18.3  2.14    8.05
 8  14.2   6.72  8.84     -2.30     0.299    -0.737   -18.3  2.14    3.34
 9  14.2   6.70  0.291    -2.30     0.303    -0.739   -18.3  2.14    9.62
10  14.2   6.70  0.291    -2.30     0.303    -0.739   -18.3  2.14    3.34
```

Three independent things were going on:

1. **Finite-sample MLE ≠ population truth.** The "known values" the user
   compared against are the parameters that *generated* the data; the
   MLE for a finite sample is a different point on the surface. On the
   shipped `examples` dataset, `loglik_math(true_par)` = `-630.51` but
   the optimizer reaches `-629.25` — the optimizer is finding a *better*
   fit to the data than the true generating parameters.

2. **The xsdm log-likelihood is invariant under a discrete equivalence
   class.** For `p` environmental variables the model is unchanged
   under any of the `2^p · p!` operations:
   * sign-flip column `k` of `o_mat` and swap `(sigltil_k, sigrtil_k)`;
   * permute the columns of `o_mat` along with the matching entries of
     `sigltil` and `sigrtil`.
   On the math scale these become non-trivial transformations of
   `sigltil*`, `sigrtil*`, and `o_par*`. Rows 2–10 above are simply
   different group images of the same `loglik = -629.25` MLE.

3. **Loose default tolerances.** `optimize_likelihood()` was previously
   shipping `grtol = 1e-4`, `xtol = 1e-8`. ucminfcpp's own
   `ucminf_control()` defaults are `grtol = 1e-6`, `xtol = 1e-12` — and
   so are R's classical `ucminf::ucminf` defaults. With central
   differences, the loose tolerances allow the optimizer to declare
   convergence at a point where `max|grad|_inf` is anywhere in
   `[1e-6, 1e-4]`, which on a shallow likelihood ridge is far from the
   actual minimum.

## Fix shipped

### Tighten default control to canonical ucminfcpp values

```r
default_ctrl <- list(
  grad     = "central",
  gradstep = c(1e-6, 1e-8),
  grtol    = 1e-6,    # was 1e-4
  xtol     = 1e-12,   # was 1e-8
  stepmax  = 5,       # kept (vs ucminfcpp default 1) — math-scale
                      # parameters live on a wide unconstrained range
  maxeval  = 2000
)
```

Reference: `R/optimize_likelihood.R` lines 78–94.

### Math-scale equivalence-class canonicalization

New exported function:

```r
canonicalize_param_vector(param_vector, reference)
canonicalize_solutions(fit, reference_par = NULL)
```

`canonicalize_param_vector(p, ref)` walks all `2^p · p!` group elements,
filters to the subset that produces `det(o_mat) = +1` (the math-scale
parameterisation `o_par` only encodes `SO(p)`, not `O(p)`, so half the
group elements have no math-scale image), picks the one with smallest
biological-scale Frobenius distance to `ref`, and returns the math-scale
representation of that closest representative.

`canonicalize_solutions(fit)` applies this to every row of
`fit$solutions$full_par` using `fit$best$par` as the reference, so all
rows of the table read in the same EC representative.

Reference: `R/canonicalize_solutions.R`,
`tests/testthat/test-canonicalize_solutions.R`.

## Empirical evidence on the shipped `examples` data

Running 12 starts on `examples$env_array / examples$occ_vec`:

| Defaults                              | Best `loglik` | Median `loglik` (12 starts) |
|---------------------------------------|---------------|------------------------------|
| Old `grtol = 1e-4`, `xtol = 1e-8`     | -629.251      | -629.251                     |
| New `grtol = 1e-6`, `xtol = 1e-12`    | -629.251      | -629.251                     |
| ucminfcpp default `stepmax = 1`       | -629.251      | -629.251                     |

(All three converge to the same MLE; the user's apparent "instability"
was the EC group acting on parameter space, not on log-likelihood.)

After `canonicalize_solutions(fit)` every row of `solutions$full_par`
reads:

```
mu1=14.1727  mu2=6.7034  sigltil1=0.3032  sigltil2=-0.7387
sigrtil1=0.2909  sigrtil2=-2.2950  ctil=-18.2637  pd=2.1425  o_par1=0.1944
```

— a single, stable answer ready to compare against `examples$par_vec`.

## How to use it

```r
fit <- optimize_likelihood(
  env_dat    = examples$env_array,
  occ        = examples$occ_vec,
  num_starts = 100L
)

# Same fit, every row of $solutions$full_par now in the same EC representative
fit_canon <- canonicalize_solutions(fit)

# Comparing to the truth:
dist_between_params(
  math_to_bio(fit_canon$best$par),
  math_to_bio(examples$par_vec)
)
```

Note: if the user wants to canonicalize against a *known* truth rather
than against the highest-loglik fit, pass `reference_par = truth_math`:

```r
fit_canon <- canonicalize_solutions(fit, reference_par = examples$par_vec)
```

## When to also pass tighter `control` overrides

For pathological data (very flat likelihood surface, near-collinear
environmental variables, very few presences), even `grtol = 1e-6` may
declare premature convergence. In those cases pass:

```r
fit <- optimize_likelihood(
  env_dat    = ...,
  occ        = ...,
  num_starts = 100L,
  control    = list(grtol = 1e-8, xtol = 1e-14, stepmax = 1, maxeval = 5000)
)
```

The user-supplied `control` wins over the package defaults.

## Pitfalls (lessons from the audit)

- **Don't compare to the population truth as the convergence check.**
  The MLE on a finite sample is not equal to the parameters that
  generated the sample. Compare the achieved log-likelihood against
  `loglik_math(true_par)`; if MLE > true (as expected), the optimization
  worked.
- **Don't conclude "the optimizer is broken" from raw `solutions$full_par`.**
  For a model with non-trivial equivalence classes, raw parameter
  vectors look spurious. Always canonicalize before inspection.
- **Don't trust loose `grtol` on a likelihood with shallow ridges.**
  `1e-4` on a Newton-ish gradient corresponds to function-value drift
  of order `O(grtol^2 / Hessian_eigenvalue)` — easily 1–10 log-likelihood
  units near a flat ridge. Use ucminfcpp's `1e-6` default or tighter.
- **Math scale only encodes `SO(p)`, not `O(p)`.** Half the bio-scale
  EC elements have no math-scale image. Either (a) skip non-`SO(p)`
  elements as the canonicalization helper does, or (b) canonicalize on
  the bio scale where the full group acts.
- **`extract_orthogonal_matrix_parameters()` (used by `bio_to_math()`)
  is the principal matrix logarithm.** It can fail for orthogonal
  matrices with eigenvalue `-1` (e.g., point reflections in `p ≥ 3`).
  `canonicalize_param_vector()` wraps the round-trip in `tryCatch` and
  falls back to the input vector if the round-trip fails, with a warning.

## Reusable checklist for a new model

- [ ] Check whether the model is invariant under any discrete group
      action on the parameter space (sign flips, permutations, gauge
      symmetries). If yes:
- [ ] Implement `convert_equivalence_class()` on the natural ("bio")
      scale of the parameters.
- [ ] Implement `canonicalize_*()` that maps each fit to a canonical
      group representative, restricted to the subset that has an image
      under the optimizer's parameter scale.
- [ ] Fix optimizer `grtol` / `xtol` to the canonical values from the
      underlying optimizer's documentation, not whatever was loose
      enough to "look fast".
- [ ] Add a parity test that the canonicalization preserves the
      objective value at every fit.
- [ ] Add a test that uses several explicit group images of a known
      reference and asserts they all canonicalize to the same vector.
