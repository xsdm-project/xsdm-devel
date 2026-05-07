# Numerically stable \`log(1 - exp(-a))\`

Computes \\\log(1 - \exp(-a))\\ accurately for non-negative \`a\`, using
two different formulas depending on whether \`a\` is above or below
\`log(2)\`.

## Usage

``` r
log1mexp(a, cutoff = log(2))
```

## Arguments

- a:

  Numeric vector of non-negative values. \`NA\` values are preserved;
  negative values emit a warning and return \`NaN\`.

- cutoff:

  Positive numeric scalar. Threshold between the two formulas;
  \`log(2)\` is near-optimal.

## Value

A numeric vector the same length as \`a\` with \\\log(1 - \exp(-a))\\.

## References

Mächler, M. (2012). \*Accurately Computing log(1 − exp(− \|a\|)).\* CRAN
package \`copula\` vignette.

## See also

[`log1pexp`](https://xsdm-project.github.io/xsdm-devel/reference/log1pexp.md),
[`log1p`](https://rdrr.io/r/base/Log.html),
[`expm1`](https://rdrr.io/r/base/Log.html)

## Examples

``` r
a <- 2^seq(-20, 5, length.out = 10)
cbind(a, log(1 - exp(-a)), log1mexp(a))
#>                  a                            
#>  [1,] 9.536743e-07 -1.386294e+01 -1.386294e+01
#>  [2,] 6.540253e-06 -1.193754e+01 -1.193754e+01
#>  [3,] 4.485274e-05 -1.001215e+01 -1.001215e+01
#>  [4,] 3.075979e-04 -8.086871e+00 -8.086871e+00
#>  [5,] 2.109492e-03 -6.162363e+00 -6.162363e+00
#>  [6,] 1.446679e-02 -4.243124e+00 -4.243124e+00
#>  [7,] 9.921257e-02 -2.359687e+00 -2.359687e+00
#>  [8,] 6.803950e-01 -7.060641e-01 -7.060641e-01
#>  [9,] 4.666116e+00 -9.453283e-03 -9.453283e-03
#> [10,] 3.200000e+01 -1.265654e-14 -1.266417e-14
```
