# Compute likelihood for LTSG model

Compute likelihood for LTSG model

Compute likelihood for LTSG model

## Arguments

- mu:

  Numeric vector of means (length equal to number of rows in \`env_m\`)

- env_m:

  Numeric matrix of environmental data. Must be column-major with time
  varying fastest: column `j` corresponds to `(location, time)` via
  index `j*q + i`. This matches the memory layout expected by the
  underlying C++ implementation.

- dl_mat:

  Diagonal matrix (as NumericMatrix)

- drl_mat:

  Diagonal matrix (as NumericMatrix)

- ortho_m:

  Numeric matrix (orthogonal basis)

- q:

  Integer, number of rows for reshaping

- r:

  Integer, number of columns for reshaping.

## Value

A numeric vector of length \`r\` with computed sums.

A numeric vector of length \`r\` with computed sums.

## Details

This function calculates a likelihood-like measure using orthogonal
matrices, environmental data, and diagonal matrices, leveraging parallel
computation.

This function calculates a likelihood-like measure using orthogonal
matrices, environmental data, and diagonal matrices, leveraging parallel
computation.

## Examples

``` r
mu <- c(1, 2)
ortho_m <- matrix(1:4, nrow = 2)
env_m <- matrix(1:4, nrow = 2)
dl_mat <- diag(2)
drl_mat <- diag(2)
like_ltsg(mu, env_m, dl_mat, drl_mat, ortho_m, q = 1, r = 2)
#> [1]   0 416
```
