# ucminfcpp migration plan for xsdm optimizers

## Scope completed in this issue

- `optimize_likelihood()` and `profile_likelihood()` now optimize through
  `ucminfcpp::ucminf_xptr()`.
- A compiled objective-pointer factory (`make_loglik_math_xptr`) was added in
  C++ and wired into the optimization helpers.
- Package metadata was updated to depend on `ucminfcpp` and C++17.

## Parameter and objective mapping

- **Free parameters (`par`)**: unchanged, still the math-scale vector of free
  parameters.
- **Fixed parameters (`mask`)**: passed into the XPtr objective closure and
  merged by `loglik_math()`.
- **Data (`env_dat`, `occ`)**: captured by the XPtr objective closure.
- **Objective value**: `loglik_math(..., negative = TRUE)` (minimization target).
- **Gradient for XPtr path**: finite differences computed inside the compiled
  objective using `control$grad` and `control$gradstep`.

## Risks and limitations

- The XPtr objective currently evaluates `loglik_math()` from C++, so heavy
  model internals in R are still present; speedups depend on workload.
- Since gradients are now evaluated inside the objective closure, numerical
  sensitivity to `gradstep` remains important.
- Requires a C++17-capable toolchain and installation of `ucminfcpp`.

## Follow-up recommendations

1. Move more of `loglik_math` / `loglik_bio` internals into compiled C++ to
   reduce R callback overhead further.
2. Add benchmark tests comparing old `ucminf` and current `ucminfcpp::ucminf_xptr`
   paths on representative datasets.
3. Add a user-facing fallback option if `ucminfcpp` is unavailable in
   constrained environments.
