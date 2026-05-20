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
env_dat <- example_1$env_array[example_1$occ_df$presence == 1, , ]
start_parms(env_dat)
#> # A tibble: 100 × 9
#>      mu1   mu2 sigltil1 sigltil2 sigrtil1 sigrtil2   ctil     pd o_par1
#>    <dbl> <dbl>    <dbl>    <dbl>    <dbl>    <dbl>  <dbl>  <dbl>  <dbl>
#>  1  8.28  4.70  -0.0275   0.243   0.0966    0.163  -1.26  -0.790  4.57 
#>  2  9.51  3.40   0.666   -0.450  -0.597     0.856  -0.776  1.41  -4.86 
#>  3 10.1   4.05  -0.374   -0.104   0.443    -0.184  -1.50  -1.89   9.28 
#>  4  8.90  2.75   0.319   -0.797  -0.250     0.510  -1.02   0.309 -0.147
#>  5  9.20  3.73  -0.201    0.0697  0.270    -0.357  -1.38  -1.34   6.92 
#>  6 10.4   5.03   0.492   -0.623  -0.423     0.336  -0.896  0.858 -2.50 
#>  7  9.82  3.08  -0.547    0.416  -0.0767   -0.0103 -1.14  -0.240  2.21 
#>  8  8.59  4.38   0.146   -0.277  -0.770     0.683  -0.656  1.96  -7.22 
#>  9  8.74  2.92   0.579   -0.537   0.00998   0.0764 -1.32  -1.06   1.03 
#> 10  9.98  4.22  -0.114    0.156  -0.683     0.770  -0.836  1.13  -8.39 
#> # ℹ 90 more rows
start_parms(env_dat, mask = c(mu2 = 5, pd = 1))
#> # A tibble: 100 × 7
#>      mu1 sigltil1 sigltil2 sigrtil1 sigrtil2   ctil o_par1
#>    <dbl>    <dbl>    <dbl>    <dbl>    <dbl>  <dbl>  <dbl>
#>  1  8.28   0.427   -0.147    0.248    0.445  -1.09  -3.98 
#>  2  9.51  -0.266    0.546   -0.445   -0.249  -0.611  5.45 
#>  3 10.1    0.0809  -0.493   -0.0983   0.791  -1.33  -8.69 
#>  4  8.90  -0.612    0.200   -0.791    0.0980 -0.851  0.736
#>  5  9.20  -0.0924  -0.320    0.0750   0.618  -1.45  -6.33 
#>  6 10.4    0.601    0.373   -0.618   -0.0753 -0.971  3.09 
#>  7  9.82  -0.439   -0.667    0.422    0.271  -1.21  -1.62 
#>  8  8.59   0.254    0.0264  -0.272   -0.422  -0.731  7.80 
#>  9  8.74  -0.526    0.460   -0.532    0.358  -1.15  -5.15 
#> 10  9.98   0.167   -0.234    0.162   -0.335  -0.671  4.27 
#> # ℹ 90 more rows
```
