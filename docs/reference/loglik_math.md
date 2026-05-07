# Log-likelihood function for the xsdm model, parameters on the math scale.

Computes the log-likelihood for the xsdm model given environmental data,
a vector of occurrences and pseudo-absences, and model parameters on the
math scale. This is the function that one optimizes to fit xsdm with
data.

## Usage

``` r
loglik_math(
  param_vector,
  env_dat,
  occ,
  mask = NULL,
  num_threads = RcppParallel::defaultNumThreads(),
  negative = TRUE
)
```

## Arguments

- param_vector:

  A \*\*named numeric vector\*\* of math-scale parameters. When
  `mask = NULL`, the names must exactly match the canonical schema
  returned by
  [`make_mask_names`](https://xsdm-project.github.io/xsdm-devel/reference/make_mask_names.md)`(p)`
  where `p = dim(env_dat)[3]`, and the length must equal
  [`num_par`](https://xsdm-project.github.io/xsdm-devel/reference/num_par.md)`(p)`.
  When `mask` is supplied, `param_vector` should contain only the names
  \*\*not\*\* present in `mask`, in the canonical order. In both cases
  the vector is combined with `mask` via
  [`create_param_vector_masked`](https://xsdm-project.github.io/xsdm-devel/reference/create_param_vector_masked.md)
  and then mapped to biological-scale parameters via
  [`math_to_bio`](https://xsdm-project.github.io/xsdm-devel/reference/math_to_bio.md).
  Must not contain missing values. See Details for the full naming and
  ordering conventions.

- env_dat:

  The environmental data array, dimensions `n_loc x n_time x p` (number
  of locations x time-series length x number of environmental
  variables). Must be a 3-dimensional array with no missing values.

- occ:

  Presence/pseudo-absence binary vector. Same length as dimension 1 of
  `env_dat`.

- mask:

  For optionally keeping some parameters at fixed values during
  optimization. Either NULL or a named numeric vector. The NULL case
  means all parameters are in param_vector, corresponding to the case
  where all parameters will be adjustable by the optimizer when this
  function is passed to it as the objective function. In the non-NULL
  case, names of entries of `mask` must correspond precisely to
  parameter names (see Details), and then those values are used. In that
  case, `param_vector` is construed to contain the values of the other
  parameters, in order (see Details). So, in particular, the length of
  `param_vector` plus the length of `mask` must equal the total number
  of parameters for the model, which is determined by dim(env_dat)\[3\].
  The most common case for most applications will be `mask=NULL`.
  Entries of `mask` are interpreted on the math scale.

- num_threads:

  Number of threads for parallel computation. Defaults to
  [`RcppParallel::defaultNumThreads()`](https://rdrr.io/pkg/RcppParallel/man/setThreadOptions.html).

- negative:

  Logical. If TRUE returns the negative of the log-likelihood instead of
  the log-likelihood itself. Facilitates optimization with some
  optimizers.

## Value

A single value, the log-likelihood (or the negative log-likelihood, if
`negative` is TRUE).

## Details

Optimizing the likelihood and profiling requires conventions for
transforming parameters from unconstrained spaces to the constrained
space of possible parameters which can be accepted by `loglik_bio`. This
function and `math_to_bio` implement those conventions, and also allow
for optimizations while keeping one or more parameters fixed, including
potentially at boundary values. Typically `loglik_math` is the function
one optimizes numerically in order to fit xsdm or a boundary model with
data, or to profile a fitted model. For what follows, denote
`dim(env_dat)[3]` by `p`.

We start by explaining the case `mask=NULL`, for which all model
parameters are in `param_vector`. The parameters of `param_vector` are
assumed to appear in the following order:

1.  Parameters for `mu`, of which there are `p`;

2.  Parameters which are `exp`-transformed to get the entries of
    `sigltil`, of which there are `p`;

3.  Parameters which are `exp`-transformed to get the entries of
    `sigrtil`, of which there are `p`;

4.  The parameter `ctil`;

5.  A parameter which is `expit` transformed to get `pd`;

6.  Parameters which are inserted via column-major order into the lower-
    triangle of a skew-symmetric matrix which is then transformed by the
    matrix exponential to get `o_mat`, of which there are (p^2-p)/2.

Thus, when `mask` is `NULL`, `param_vector` must be an unconstrained
numeric vector of length `3*p+2+(p^2-p)/2` with no missing values.

The argument `mask` is used in the event one wants to fix certain
parameters and optimize over the remaining parameters. This argument
must be a named numeric vector with unique names being some but not all
of the `3*p+2+(p^2-p)/2` following: `mu1`, `mu2`, ..., `mup`,
`sigltil1`, `sigltil2`, ..., `sigltilp`, `sigrtil1`, `sigrtil2`, ...,
`sigrtilp`, `ctil`, `pd`, and `o_mati` for `i` ranging from 1 to
`(p^2-p)/2`. These names must be used exactly. See the function
`make_mask_names`, which facilitates the construction of a correctly
formatted `mask` argument. Entries of `mask` are on the math scale;
`bio_to_math` can convert biological-scale constraints to math scale.

The missing entries of `mask` are filled in using the entries of
`param_vector`, in the order specified above, and then the
transformations described above (implemented by `math_to_bio`) are
applied to get biological-scale parameters which are passed to
`loglik_bio` to get the log likelihood.

Entries of `mask` corresponding to `sigltil` or `sigrtil` can be `Inf`.
Likewise, the entry of `mask` corresponding to `pd` can be `Inf` (on the
math scale, corresponding to a biological-scale value of 1). This
functionality is used to fit boundary models. Entries of `param_vector`
must be finite.

## Implementation

This is a thin R wrapper around the C++ implementation
`loglik_math_cpp`; the optimizer hot path is pure C++. A pure-R
reference, `loglik_math_r`, is kept internal to the package and is used
only by the parity tests in
`tests/testthat/test-loglik_math_r_vs_cpp.R`.

## Examples

``` r
# Testing the function with the example data
loglik_math(
  param_vector = example_1$par_vec,
  env_dat = example_1$env_array,
  occ = example_1$occ_vec
)
#> [1] 630.5103
# Mute one parameter to use the mask
par_vec <- example_1$par_vec[-2]
mask_parameters_a <- c(mu2 = 6.5)
loglik_math(
  param_vector = par_vec,
  env_dat = example_1$env_array,
  occ = example_1$occ_vec,
  mask = mask_parameters_a
)
#> [1] 789.1959
# Return the negative
loglik_math(
  param_vector = example_1$par_vec,
  env_dat = example_1$env_array,
  occ = example_1$occ_vec,
  negative = TRUE
)
#> [1] 630.5103
```
