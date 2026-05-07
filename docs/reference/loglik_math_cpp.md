# Pure-C++ log-likelihood for the xsdm model (math-scale parameters)

Computes the log-likelihood directly in C++ without any R callback in
the inner loop. Semantically equivalent to the R function loglik_math.

## Usage

``` r
loglik_math_cpp(
  param_vector,
  env_dat,
  occ,
  mask = NULL,
  negative = TRUE,
  num_threads = 0L
)
```

## Arguments

- param_vector:

  Named numeric vector of math-scale parameters. When \`mask\` is NULL,
  must contain every canonical name for the dimension p implied by
  \`env_dat\`. When \`mask\` is supplied, contains only the free
  (non-masked) parameters.

- env_dat:

  3D numeric array with dimensions (n_loc, ts_length, p). No missing
  values allowed.

- occ:

  Integer or logical vector of length n_loc, 0/1 or FALSE/TRUE.

- mask:

  Optional named numeric vector of fixed parameters.

- negative:

  Logical; if TRUE (default) returns the negative log-likelihood (the
  value to be minimized).

- num_threads:

  Integer; 0 leaves the RcppParallel default.

## Value

A scalar double.
