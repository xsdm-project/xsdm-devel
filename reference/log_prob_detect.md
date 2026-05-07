# Probability of detection of the species in each location

Computes the probability of detection of the species in each location
for the xsdm model, given environmental data and model parameters.

## Usage

``` r
log_prob_detect(
  env_dat,
  mu,
  sigltil,
  sigrtil,
  o_mat,
  ctil,
  pd,
  return_prob = FALSE,
  num_threads = RcppParallel::defaultNumThreads()
)
```

## Arguments

- env_dat:

  The environmental data array, dimensions `n_loc x n_time x p` (number
  of locations x time-series length x number of environmental
  variables). Must be a 3-dimensional array with no missing values.

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

  Logical (default FALSE). Flag to return probabilities of detection
  instead their logs.

- num_threads:

  Number of threads for parallel computation. Defaults to
  [`RcppParallel::defaultNumThreads()`](https://rdrr.io/pkg/RcppParallel/man/setThreadOptions.html).

## Value

A vector of length equal to the number of locations, containing the
probabilities of detection (or their logs) of the species in each
location.

## Details

This is a thin R wrapper around the C++ implementation
`log_prob_detect_cpp`; the optimizer hot path is pure C++. A pure-R
reference implementation, `log_prob_detect_r`, is kept internal to the
package and is used only by the parity tests in
`tests/testthat/test-log_prob_detect_r_vs_cpp.R`.

## Examples

``` r
mu <- c(-1, 5.046939)
sigltil <- c(1.036834, 1.556083)
sigrtil <- c(1.538972, 1.458738)
ctil <- -2
pd <- 0.9
o_mat <- matrix(c(-0.4443546, 0.8958510, -0.8958510, -0.4443546), ncol = 2)
```
