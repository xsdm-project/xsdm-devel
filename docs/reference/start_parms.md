# Starting parameters for the optimization

Generate starting parameters for the optimization of the xsdm
log-likelihood. Starting points are constructed from the environmental
conditions at observed presences (where `occ == 1`) using a Latin
hypercube design for the parameters based on the Sobol' low-discrepancy
sequence.

## Usage

``` r
start_parms(env_dat, mask = NULL, breadth = 1, num_starts = 100)
```

## Arguments

- env_dat:

  The environmental array for only the observed occurrences

- mask:

  Either NULL or a named numeric vector. Names must be as specified by
  calling `make_mask_names`. The NULL case means planned optimizations
  will be over all model parameters, so start parameter sets should
  include all model parameters. The non-NULL case means some parameters
  will not be specified in the output of this function because the
  optimizations which are planned will fix those parameters anyway. The
  most common case for most applications will be `mask=NULL`.

- breadth:

  Scalar in `[0, 1]` controlling how wide the search ranges are around
  their (fixed, data-driven) center. `breadth = 1` (the default)
  reproduces the pre-v0.3 behaviour that corresponded to
  `quant_vec = c(0.1, 0.5, 0.9)`; `breadth = 0` collapses every range to
  essentially a single point, equivalent to
  `quant_vec = c(0.5 - 1e-6, 0.5, 0.5 + 1e-6)`. Values in between
  interpolate linearly.

- num_starts:

  The number of samples of the hypercube

## Value

A data frame with samples for each parameter to optimize

## Details

The bounds and center of the search range for a mu parameter are based
on the quantiles in quant_vec applied to all observations of that
environmental variable, over space and time. The relationship between
quant_vec and the width of the ranges selected for the other parameters
varies, but generally wider ranges in quant_vec produce wider ranges for
start parameters.

## Examples

``` r
env_dat <- example_1$env_array[example_1$occ_vec == 1, , ]
start_parms(env_dat)
#> # A tibble: 100 × 9
#>      mu1   mu2 sigltil1 sigltil2 sigrtil1 sigrtil2   ctil     pd o_par1
#>    <dbl> <dbl>    <dbl>    <dbl>    <dbl>    <dbl>  <dbl>  <dbl>  <dbl>
#>  1  9.58  4.00   0.364    0.359   0.172     0.620  -1.26   0.858 -4.86 
#>  2 11.8   5.29   1.06    -0.334   0.865    -0.0736 -0.624 -1.34   4.57 
#>  3 12.9   3.36   0.711    0.0127  1.21      0.273  -0.944  1.96  -0.147
#>  4 10.7   4.65   0.0178  -0.680   0.518    -0.420  -1.58  -0.240  9.28 
#>  5 11.2   3.04   0.884    0.533   1.04     -0.247  -1.42  -1.89  -7.22 
#>  6 13.4   4.33   0.191   -0.161   0.345     0.446  -0.784  0.309  2.21 
#>  7 12.3   3.68  -0.155    0.186  -0.00144   0.0997 -1.10  -0.790 -2.50 
#>  8 10.1   4.97   0.538   -0.507   0.692     0.793  -1.74   1.41   6.92 
#>  9 10.4   3.52  -0.0688  -0.0740  0.778     0.706  -1.18  -1.06  -8.39 
#> 10 12.6   4.81   0.624    0.619   0.0852    0.0130 -0.544  1.13   1.03 
#> # ℹ 90 more rows
start_parms(env_dat, mask = c(mu2 = 5, pd = 1))
#> # A tibble: 100 × 7
#>      mu1 sigltil1 sigltil2 sigrtil1 sigrtil2   ctil o_par1
#>    <dbl>    <dbl>    <dbl>    <dbl>    <dbl>  <dbl>  <dbl>
#>  1  9.58   0.299   -0.139    0.995   -0.290  -0.664 -1.91 
#>  2 11.8    0.993    0.554    0.302    0.403  -1.30   7.51 
#>  3 12.9   -0.0472   0.208    0.648    0.750  -0.984  2.80 
#>  4 10.7    0.646   -0.486   -0.0448   0.0564 -1.62  -6.63 
#>  5 11.2   -0.220    0.381    1.17     0.576  -1.46  -4.27 
#>  6 13.4    0.473   -0.312    0.475   -0.117  -0.824  5.15 
#>  7 12.3    0.126   -0.659    0.822   -0.463  -1.14   0.442
#>  8 10.1    0.819    0.0343   0.129    0.230  -0.504 -8.98 
#>  9 10.4    0.0395  -0.572    0.562    0.316  -0.584 -0.736
#> 10 12.6    0.733    0.121    1.25    -0.377  -1.22   8.69 
#> # ℹ 90 more rows
```
