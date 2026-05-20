# Basic (non-adaptive) tool for profiling the likelihood

Basic (non-adaptive) tool for profiling the likelihood

## Usage

``` r
profile_likelihood(
  profile_parameter = "mu1",
  increment_left = 0.1,
  increment_right = increment_left,
  num_steps_left = 20L,
  num_steps_right = num_steps_left,
  alpha = 0.95,
  optim_param_vector,
  env_dat,
  occ,
  mask = NULL,
  num_threads = RcppParallel::defaultNumThreads(),
  control = list(),
  verbose = FALSE
)
```

## Arguments

- profile_parameter:

  Character. Name of the parameter to profile. Profiles are done on the
  math scale.

- increment_left:

  Numeric. Step size (math scale) when moving to the left, from the
  start point of the parameter point estimate, to construct the profile.

- increment_right:

  Numeric. Step size (math scale) when moving to the right, from the
  start point of the paraneter point estimate, to construct the profile.

- num_steps_left:

  Integer. Maximum number of steps to take to the left.

- num_steps_right:

  Integer. Maximum number of steps to take to the right.

- alpha:

  Numeric value between 0 and 1. Confidence level used for the
  likelihood ratio (LR) threshold: threshold = MLE_loglik -
  qchisq(alpha, 1)/2.

- optim_param_vector:

  Named numeric. MLE parameters on math scale.

- env_dat:

  3D array (locations x time x variables).

- occ:

  Logical, either 0 or 1, vector (length = number of locations).

- mask:

  Named numeric or NULL. Parameters kept fixed (math scale).

- num_threads:

  Integer. Threads used internally by log-likelihood.

- control:

  Named list. Control passed to `ucminfcpp::ucminf_xptr(control = ...)`.
  User-specified entries override defaults:

  - `grad = "central"`

  - `gradstep = c(1e-6, 1e-8)`

  - `grtol = 1e-5`

  - `xtol = 1e-12`

  - `stepmax = 5`

  - `maxeval = 2000`

  If you want optimizer iteration trace, set `control$trace > 0`.

- verbose:

  Logical. If `TRUE`, prints compact progress messages; otherwise
  silent.

## Value

A list with:

- `profile`: data.frame with columns `param`, `value_math`, `loglik`,
  `convergence`, and a list-column `full_par`. (No side/step columns.)

- `found_better`: logical; TRUE if any profiled point exceeds the MLE
  log-likelihood.

- `threshold`: numeric; LR threshold used.

- `parameters`: data.frame with the parameters found in each step of the
  profiling

## Examples

``` r
## Minimal profiling example (fast): 1 step left + 1 step right
res <- profile_likelihood(
  profile_parameter = "mu1",
 increment_left = 0.2,
 increment_right = 0.2,
 num_steps_left = 1L, # one iteration on the left
 num_steps_right = 1L, # one iteration on the right
 alpha = 0.95,
 optim_param_vector = example_1$optim_par_vec,
 env_dat = example_1$env_array,
 occ = example_1$occ_df$presence,
 num_threads = 1L, # keep it fast and deterministic
 control = list(maxeval = 20),
 verbose = FALSE
)
# Check the structure of the output:
res$profile
#>   param value_math    loglik convergence
#> 2   mu1   8.730929 -1018.311           3
#> 1   mu1   8.930929 -1009.447          NA
#> 3   mu1   9.130929 -1010.026           3
res$threshold
#> [1] -1011.367
res$found_better
#> [1] FALSE
## Full math-scale parameter vectors used at each evaluated point:
res$parameter_df
#> NULL
```
