# Helper. Profile one side of a likelihood profile (internal)

Fixes one parameter on the math scale and re-optimizes the remaining
parameters along a single direction (left/right) until the LR threshold
is reached or a step cap is hit.

## Usage

``` r
profile_one_side_(
  direction,
  increment,
  max_steps,
  profile_parameter,
  optim_param_vector,
  env_dat,
  occ,
  mask,
  num_threads,
  optim_ll,
  thresh,
  base_control,
  start_full = optim_param_vector,
  invh_lt = NULL,
  verbose = FALSE
)
```

## Arguments

- direction:

  Integer. -1 (left) or +1 (right).

- increment:

  Numeric. Step size on math scale for this side.

- max_steps:

  Integer. Maximum iterations for this side.

- profile_parameter:

  Character. Name of the parameter to profile. Profiles are done on the
  math scale.

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

- base_control:

  Named list. Control passed to `ucminfcpp::ucminf_xptr(control = ...)`.
  User-specified entries should be merged in the caller (see
  `profile_likelihood`).

- start_full:

  Named numeric. Full warm-start parameter vector. Defaults to
  `optim_param_vector`.

- invh_lt:

  Optional numeric. Lower triangle of the inverse Hessian for
  warm-start.

- verbose:

  Logical. If `TRUE`, prints compact progress messages; otherwise
  silent.

## Value

A list with elements `ll`, `vals`, `fulls`, `conv`, `last_full`,
`last_invh`, `steps`, and `crossed`.

## See also

[`profile_likelihood`](https://xsdm-project.github.io/xsdm-devel/reference/profile_likelihood.md)
