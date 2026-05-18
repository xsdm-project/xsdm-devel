# Changelog

## xsdm 1.0.0

First public release.

### Highlights

- **End-to-end C++ likelihood chain.** The optimization hot path
  (`loglik_math`, `loglik_bio`, `log_prob_detect`, `math_to_bio`,
  `create_param_vector_masked`, `like_neg_ltsgr`) is now driven entirely
  by C++ kernels via thin R wrappers. The pure-R reference
  implementations are preserved as non-exported `*_r` functions in
  `R/internals.R` and used by the test suite to assert numerical parity
  at `tol = 1e-6`.
- **[`optimize_likelihood()`](https://xsdm-project.github.io/xsdm-devel/reference/optimize_likelihood.md)
  and
  [`profile_likelihood()`](https://xsdm-project.github.io/xsdm-devel/reference/profile_likelihood.md)**
  drive
  [`ucminfcpp::ucminfcpp_xptr()`](https://CRAN.R-project.org/package=ucminfcpp)
  with a pure-C++ XPtr factory (`make_loglik_math_xptr_cpp`), so there
  are no R callbacks in the inner loop.
- **Parameter masking** is fully integrated: any subset of the canonical
  parameter vector can be fixed during optimization, and the remaining
  free parameters are optimized end-to-end on the math scale.
- **Profile-likelihood name propagation fix.**
  [`profile_likelihood()`](https://xsdm-project.github.io/xsdm-devel/reference/profile_likelihood.md)
  now preserves free-parameter names across the `ucminf` boundary, so
  profile rows are reliably labelled with the parameter being held
  fixed.

### API changes

- New canonical exported names (with thin C++ wrappers):

  | exported                     | non-exported reference                |
  |------------------------------|---------------------------------------|
  | `loglik_math`                | `xsdm:::loglik_math_r`                |
  | `loglik_bio`                 | `xsdm:::loglik_bio_r`                 |
  | `log_prob_detect`            | `xsdm:::log_prob_detect_r`            |
  | `math_to_bio`                | `xsdm:::math_to_bio_r`                |
  | `create_param_vector_masked` | `xsdm:::create_param_vector_masked_r` |
  | `like_neg_ltsgr`             | `xsdm:::like_neg_ltsgr_r`             |
  | `dist_between_params`        | `xsdm:::dist_between_params_r`        |

- `like_neg_ltsgr_cpp` is no longer exported; it remains as an
  unexported back-compat alias (`xsdm:::like_neg_ltsgr_cpp`) that
  forwards to `like_neg_ltsgr`. New code should call `like_neg_ltsgr`
  directly.

- [`dist_between_params()`](https://xsdm-project.github.io/xsdm-devel/reference/dist_between_params.md)
  is now backed by a pure-C++ implementation in
  `src/dist_between_params.cpp` that builds the pairing cost matrix and
  solves the linear sum assignment problem via a clean-room Hungarian /
  Kuhn–Munkres routine (O(n^3) potentials variant from Kuhn 1955,
  Munkres 1957, Jonker & Volgenant 1987). No code is taken from the
  `clue` package, avoiding the GPL-2 / AGPL-3 licence mismatch. The
  pure-R reference `xsdm:::dist_between_params_r` continues to call
  [`clue::solve_LSAP()`](https://rdrr.io/pkg/clue/man/solve_LSAP.html)
  and is used by the parity tests in
  `tests/testthat/test-dist_between_params_r_vs_cpp.R`. The legacy
  brute-force reference `distance_between_params()` is also preserved
  (non-exported) as `xsdm:::distance_between_params_r`.

### Validation & robustness

- [`optimize_likelihood()`](https://xsdm-project.github.io/xsdm-devel/reference/optimize_likelihood.md),
  [`start_parms()`](https://xsdm-project.github.io/xsdm-devel/reference/start_parms.md),
  and
  [`get_start_parms_()`](https://xsdm-project.github.io/xsdm-devel/reference/get_start_parms_.md)
  now reject `num_starts < 3` early with an informative `checkmate`
  error, rather than crashing inside
  [`sobol::sobol_design()`](https://alrobles.github.io/sobol/reference/sobol_design.html)
  (which segfaults at `nseq = 0` and returns a malformed vector at
  `nseq = 1`).
- [`start_parms()`](https://xsdm-project.github.io/xsdm-devel/reference/start_parms.md)
  runtime:
  [`sobol::sobol_design`](https://alrobles.github.io/sobol/reference/sobol_design.html)
  output is defensively coerced to a `data.frame`, protecting against
  future behaviour changes.

### Build & packaging

- `terra` and `tibble` moved from `Suggests:` to `Imports:`. Both are
  used unconditionally in mainline package code
  ([`env_data_array()`](https://xsdm-project.github.io/xsdm-devel/reference/env_data_array.md),
  [`habitat_suitability()`](https://xsdm-project.github.io/xsdm-devel/reference/habitat_suitability.md),
  [`start_parms()`](https://xsdm-project.github.io/xsdm-devel/reference/start_parms.md),
  [`vsp()`](https://xsdm-project.github.io/xsdm-devel/reference/vsp.md)).
- `src/*.hpp` headers renamed to `src/*.h` to satisfy `R CMD check`’s
  “Subdirectory ‘src’ contains” warning.
- New `R-CMD-check.yaml` GitHub workflow runs `R CMD check --as-cran` on
  Ubuntu (release / devel / oldrel-1 / oldrel-2), Windows, and macOS.
- Test suite: 1006 tests, full suite ~9s, all green; parity tests live
  under `tests/testthat/test-*_r_vs_cpp.R`.

### Internals

- All non-exported helpers (`check_env_array`, `logit`, `permutations`,
  `auto_plot_lims_`, etc.) consolidated in `R/internals.R`.
- [`profile_likelihood()`](https://xsdm-project.github.io/xsdm-devel/reference/profile_likelihood.md)
  cleanly stops when the requested parameter is outside the masked
  design.
