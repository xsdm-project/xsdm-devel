# Functions to take the expit of numerical vectors. expit exp(x)/(1 + exp(x))

Functions to take the expit of numerical vectors. expit exp(x)/(1 +
exp(x))

## Usage

``` r
expit(x)
```

## Arguments

- x:

  A numeric value

## Value

A real vector corresponding to the expits of x

## Examples

``` r
expit(0)
#> [1] 0.5
expit(0.5)
#> [1] 0.6224593
expit(-1)
#> [1] 0.2689414
```
