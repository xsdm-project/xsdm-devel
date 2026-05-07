# Get the number of environmental variables given the number of parameters

Inverts `num_par(p)` to recover `p` from `n`, the number of parameters
of the main xsdm model. Uses the closed-form solution of the quadratic:
\\2n = p^2 + 5p + 4\\, i.e. \\p = (-5 + \sqrt{9 + 8n})/2\\. Errors if
\\n\\ is not a valid value of `num_par(p)` for some integer \\p \ge 1\\.

## Usage

``` r
num_env_var(n)
```

## Arguments

- n:

  Integerish scalar: total number of parameters.

## Value

A single integer `p`, the number of environmental variables.

## Examples

``` r
num_env_var(5) # -> 1  (since num_par(1) = 5)
#> [1] 1
num_env_var(9) # -> 2  (since num_par(2) = 9)
#> [1] 2
num_env_var(14) # -> 3  (since num_par(3) = 14)
#> [1] 3
# round-trip check:
p <- 4
stopifnot(num_env_var(num_par(p)) == p)
```
