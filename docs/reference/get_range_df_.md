# Helper. Reasonable ranges for initial conditions for optimizations seeking to maximize the likelihood of the xsdm model

Given environmental data, constructs reasonable ranges for each xsdm
model parameter within which to select initial conditions for the
optimizer.

## Usage

``` r
get_range_df_(env_dat, breadth = 1, quant_vec = NULL)
```

## Arguments

- env_dat:

  A 3D numeric array of environmental time series data with dimensions
  `(locations) x (time) x (environmental variables)`. Missing values are
  not allowed. Assumed to pertain only to locations where the species
  was observed.

- breadth:

  Scalar in `[0, 1]` controlling how wide the search ranges are around
  their (fixed, data-driven) center. `breadth = 1` (the default)
  reproduces the pre-v0.3 behaviour that corresponded to
  `quant_vec = c(0.1, 0.5, 0.9)`; `breadth = 0` collapses every range to
  essentially a single point, equivalent to
  `quant_vec = c(0.5 - 1e-6, 0.5, 0.5 + 1e-6)`. Values in between
  interpolate linearly.

- quant_vec:

  Deprecated. If supplied, a deprecation warning is emitted and
  `breadth` is set to `(quant_vec[3] - quant_vec[1]) / 0.8`. This
  argument will be removed in a future release.

## Value

A data.frame with three columns (lower bound, center, upper bound)
giving the search range for each parameter.

## Details

The center of every range is fixed and data-driven: it is the empirical
median for `mu`-type parameters, `0` for orthogonal-matrix angles, `0.5`
for `pd`, and the median of the evaluated log-likelihood surface for
`ctil`. Only the half-width around that center is user-controllable, via
`breadth`.

Internally, `breadth` is mapped to two monotonic quantities:

- `half_width = 1e-6 + breadth * (0.4 - 1e-6)` — the half-width in
  probability space used for quantiles of `mu`, `pd`, and `ctil`.
  Quantile arguments are clamped so that `logit` stays finite.

- `fact = 1 + breadth` — the multiplicative factor on the log scale for
  `sigltil` and `sigrtil` ranges. `breadth = 1` gives `fact = 2` (the
  legacy value); `breadth = 0` gives `fact = 1`, i.e. a degenerate
  single-point range.

## Examples

``` r
set.seed(1)
env <- array(rnorm(10 * 5 * 2), dim = c(10, 5, 2))
xsdm:::get_range_df_(env)                 # default breadth = 1
#>               lower        center      upper
#> mu1      -0.8219844  0.1291041538  0.9594551
#> mu2      -1.1388880  0.1137973312  1.4362768
#> sigltil1 -0.7837805 -0.0906332806  0.6025139
#> sigltil2 -0.7809272 -0.0877800581  0.6053671
#> sigrtil1 -1.0179876 -0.3248403828  0.3683068
#> sigrtil2 -0.6927734  0.0003737865  0.6935210
#> ctil     -1.5860885 -0.8513780741 -0.5249049
#> pd       -2.1972246  0.0000000000  2.1972246
#> o_par1   -9.4247780  0.0000000000  9.4247780
xsdm:::get_range_df_(env, breadth = 0.3) 
#>                lower        center       upper
#> mu1      -0.07583513  0.1291041538  0.44434919
#> mu2      -0.19940643  0.1137973312  0.35210235
#> sigltil1 -0.35299755 -0.0906332806  0.17173098
#> sigltil2 -0.35014432 -0.0877800581  0.17458421
#> sigrtil1 -0.58720465 -0.3248403828 -0.06247612
#> sigrtil2 -0.26199048  0.0003737865  0.26273805
#> ctil     -0.96471193 -0.8513780741 -0.69016869
#> pd       -0.48955120  0.0000000000  0.48955120
#> o_par1   -9.42477796  0.0000000000  9.42477796
```
