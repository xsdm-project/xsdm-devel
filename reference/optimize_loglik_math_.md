# Internal helper: run ucminfcpp for one starting vector

Internal helper: run ucminfcpp for one starting vector

## Usage

``` r
optimize_loglik_math_(
  param_vector,
  env_dat,
  occ,
  mask,
  num_threads,
  base_control,
  invh_lt = NULL,
  optimizer_fun = ucminfcpp::ucminf_xptr
)
```

## Value

A list with `par`, `value`, `convergence`, and optionally
`invhessian.lt`.
