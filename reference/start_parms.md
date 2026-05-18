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
#>  1  9.58  5.09   0.429    0.316    0.800    0.0564 -1.40  -0.790  4.57 
#>  2 11.8   3.80   1.12    -0.377    0.107    0.750  -0.764  1.41  -4.86 
#>  3 12.9   4.45   0.0828  -0.0306   1.15    -0.290  -1.72  -1.89   9.28 
#>  4 10.7   3.16   0.776   -0.724    0.453    0.403  -1.08   0.309 -0.147
#>  5 11.2   4.12   0.256    0.143    0.973   -0.463  -1.56  -1.34   6.92 
#>  6 13.4   5.41   0.949   -0.550    0.280    0.230  -0.924  0.858 -2.50 
#>  7 12.3   3.48  -0.0905   0.489    0.627   -0.117  -1.24  -0.240  2.21 
#>  8 10.1   4.77   0.603   -0.204   -0.0664   0.576  -0.604  1.96  -7.22 
#>  9 10.4   3.32   1.04    -0.464    0.713   -0.0303 -1.48  -1.06   1.03 
#> 10 12.6   4.61   0.343    0.229    0.0202   0.663  -0.844  1.13  -8.39 
#> # ℹ 90 more rows
start_parms(env_dat, mask = c(mu2 = 5, pd = 1))
#> # A tibble: 100 × 7
#>      mu1 sigltil1 sigltil2 sigrtil1 sigrtil2   ctil o_par1
#>    <dbl>    <dbl>    <dbl>    <dbl>    <dbl>  <dbl>  <dbl>
#>  1  9.58   0.884   -0.0740   0.952   0.338   -1.18  -3.98 
#>  2 11.8    0.191    0.619    0.258  -0.355   -0.544  5.45 
#>  3 12.9    0.538   -0.421    0.605   0.685   -1.50  -8.69 
#>  4 10.7   -0.155    0.273   -0.0881 -0.00862 -0.864  0.736
#>  5 11.2    0.364   -0.247    0.778   0.511   -1.66  -6.33 
#>  6 13.4    1.06     0.446    0.0852 -0.182   -1.02   3.09 
#>  7 12.3    0.0178  -0.594    1.12    0.165   -1.34  -1.62 
#>  8 10.1    0.711    0.0993   0.432  -0.528   -0.704  7.80 
#>  9 10.4   -0.0688   0.533    0.172   0.251   -1.26  -5.15 
#> 10 12.6    0.624   -0.161    0.865  -0.442   -0.624  4.27 
#> # ℹ 90 more rows
```
