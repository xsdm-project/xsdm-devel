# Log-likelihood function for the xsdm model, parameters on the biological scale.

Computes the log-likelihood for the xsdm model given environmental data,
a vector of occurrences and pseudo-absences, and model parameters on the
biological scale.

## Usage

``` r
loglik_bio(
  env_dat,
  occ,
  mu,
  sigltil,
  sigrtil,
  o_mat,
  ctil,
  pd,
  return_prob = FALSE,
  sum_log_p = TRUE,
  num_threads = RcppParallel::defaultNumThreads()
)
```

## Arguments

- env_dat:

  The environmental data array, dimensions `n_loc x n_time x p` (number
  of locations x time-series length x number of environmental
  variables). Must be a 3-dimensional array with no missing values.

- occ:

  Presence/pseudo-absence binary vector. Same length as dimension 1 of
  `env_dat`.

- mu:

  Vector of optimal environmental values. Length `p=dim(env_dat)[3]`.
  Unconstrained real numbers.

- sigltil:

  Vector specifying width of the growth-environment function. Length
  `p=dim(env_dat)[3]`. Positive real numbers, `Inf` entries also
  allowed.

- sigrtil:

  Vector specifying width of the growth-environment function. Length
  `p=dim(env_dat)[3]`. Positive real numbers, `Inf` entries also
  allowed.

- o_mat:

  An orthogonal matrix, dimensions `p` by `p`.

- ctil:

  Scalar. Relates to the center of the detection-link function.

- pd:

  Maximum probability of detection of the species. Parameter between 0
  and 1.

- return_prob:

  Logical (default FALSE). Flag to return likelihood instead of
  log-likelihood.

- sum_log_p:

  Logical (default TRUE). If FALSE, returns the individual
  log-likelihoods (or likelihoods, if `return_prob` is TRUE) associated
  with the individual locations, instead of their sum (product, if
  `return_prob` is TRUE).

- num_threads:

  Number of threads for parallel computation. Defaults to
  [`RcppParallel::defaultNumThreads()`](https://rdrr.io/pkg/RcppParallel/man/setThreadOptions.html).

## Value

A single value, the log-likelihood (or the likelihood, if `return_prob`
is TRUE); or a vector of location specific values of `sum_log_p` is
FALSE.

## Details

This is a thin R wrapper around the C++ implementation `loglik_bio_cpp`;
the optimizer hot path (`sum_log_p = TRUE`, `return_prob = FALSE`) is
pure C++. The non-default flag combinations (`sum_log_p = FALSE` or
`return_prob = TRUE`) are computed by delegating to the C++-backed
`log_prob_detect` and reducing in R. A pure-R reference implementation,
`loglik_bio_r`, is kept internal to the package and is used only by the
parity tests in `tests/testthat/test-loglik_bio_r_vs_cpp.R`.

## Examples

``` r
ll <- loglik_bio(
  env_dat = example_1$env_array,
  occ = example_1$occ_vec,
  mu = example_1$true_par_list$mu,
  sigltil = example_1$true_par_list$sigltil,
  sigrtil = example_1$true_par_list$sigrtil,
  o_mat = example_1$true_par_list$o_mat,
  ctil = example_1$true_par_list$ctil,
  pd = example_1$true_par_list$pd
)
ll
#> [1] -15829.08
```
