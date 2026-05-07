# Function to facilitate the creation of the argument `mask` to the function `loglik_math`

The argument `mask` to the function `loglik_math` is required to follow
some very specific conventions in order to reduce the risk of errors
coming from mismatched arguments. This function facilitates the creation
of such vectors.

## Usage

``` r
make_mask_names(p)
```

## Arguments

- p:

  A positive integer representing the number of environmental variables
  to be used in the xsdm model, i.e., `dim(env_dat)[3]` for the
  `env_dat` argument of the `loglik_math` function

## Value

A named numeric vector full of NAs, with the names generated according
to the conventions in Details of the function `loglik_math`. See also
Details below.

## Details

The output has length `3*p+(p^2-p)/2+2`. The names of the entries are
`mu1`, `mu2`, ..., `mup`, `sigltil1`, `sigltil2`, ..., `sigltilp`,
`sigrtil1`, `sigrtil2`, ..., `sigrtilp`, `o_pari` for `i` ranging from 1
to `(p^2-p)/2`, `ctil`, and `pd`. All entries are `NA`.

## Examples

``` r
make_mask_names(2)
#>      mu1      mu2 sigltil1 sigltil2 sigrtil1 sigrtil2     ctil       pd 
#>       NA       NA       NA       NA       NA       NA       NA       NA 
#>   o_par1 
#>       NA 
```
