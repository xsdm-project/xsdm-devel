# Numerically stable \`log(1 + exp(x))\`

Computes \\\log(1 + \exp(x))\\ accurately for any real \`x\`, avoiding
overflow as \`x -\> +Inf\` and catastrophic cancellation as \`x -\>
-Inf\`.

## Usage

``` r
log1pexp(x, c0 = -37, c1 = 18, c2 = 33.3)
```

## Arguments

- x:

  Numeric vector. \`NA\` values are preserved.

- c0, c1, c2:

  Numeric scalars defining the switch points between four asymptotically
  optimal formulas. Defaults (-37, 18, 33.3) are from Mächler (2012) and
  should not normally be changed.

## Value

A numeric vector the same length as \`x\` with \\\log(1 + \exp(x))\\.

## References

Mächler, M. (2012). \*Accurately Computing log(1 − exp(− \|a\|)).\* CRAN
package \`copula\` vignette.

## See also

[`log1mexp`](https://xsdm-project.github.io/xsdm-devel/reference/log1mexp.md),
[`log1p`](https://rdrr.io/r/base/Log.html),
[`expm1`](https://rdrr.io/r/base/Log.html)

## Examples

``` r
x <- seq(-40, 40, by = 10)
cbind(x, log1p(exp(x)), log1pexp(x))
#>         x                          
#>  [1,] -40 4.248354e-18 4.248354e-18
#>  [2,] -30 9.357623e-14 9.357623e-14
#>  [3,] -20 2.061154e-09 2.061154e-09
#>  [4,] -10 4.539890e-05 4.539890e-05
#>  [5,]   0 6.931472e-01 6.931472e-01
#>  [6,]  10 1.000005e+01 1.000005e+01
#>  [7,]  20 2.000000e+01 2.000000e+01
#>  [8,]  30 3.000000e+01 3.000000e+01
#>  [9,]  40 4.000000e+01 4.000000e+01
```
