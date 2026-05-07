# Pure-C++ log-likelihood for the xsdm model (biological-scale parameters)

Computes the log-likelihood directly in C++ without any R callback.
Equivalent to loglik_bio(..., sum_log_p = TRUE, return_prob = FALSE).

## Usage

``` r
loglik_bio_cpp(
  env_dat_vec,
  env_dat_dims,
  occ,
  mu,
  sigltil,
  sigrtil,
  o_mat,
  ctil,
  pd,
  num_threads = 0L
)
```

## Arguments

- env_dat_vec:

  Flat numeric vector containing env_dat in column-major order (as
  produced by as.vector(env_dat)).

- env_dat_dims:

  Integer vector of length 3: c(n_loc, ts_length, p).

- occ:

  Integer vector of length n_loc, 0 or 1.

- mu:

  Numeric vector, length p.

- sigltil:

  Positive numeric vector, length p.

- sigrtil:

  Positive numeric vector, length p.

- o_mat:

  A p x p orthogonal matrix (column-major).

- ctil:

  Scalar.

- pd:

  Scalar in (0, 1\].

- num_threads:

  Number of threads for the inner xtensor kernel (0 = RcppParallel
  default).

## Value

Scalar log-likelihood.
