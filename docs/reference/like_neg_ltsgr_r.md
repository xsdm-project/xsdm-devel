# Long-term stochastic growth rate worker function for the xsdm model, R version

Computes the negative of the long-term stochastic growth rate, plus
log(lambda_max), for the xsdm model, for each location. This is the R
version of a worker function, see also the accompanying C version, which
should produce identical results but faster.

## Usage

``` r
like_neg_ltsgr_r(env_dat, mu, sigltil, sigrtil, o_mat)
```

## Arguments

- env_dat:

  The environmental data array, dimensions (number of locations) x (time
  series length) x (number of environmental variables). Must not contain
  missing values.

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

## Value

A vector of length equal to the number of locations, as described above.

## Details

Being an internal function, there is no error checking. Note that
env_dat must be a 3d array (not a matrix or a vector) even if one of its
dimensions is 1. And `o_mat` must be a matrix even when `p` is 1 (in
that case it's a 1 x 1 matrix, but not a scalar).
