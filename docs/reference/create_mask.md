# Create a parameter mask aligned with \`make_mask_names()\`

Constructs a named numeric vector whose names and length match the
canonical schema returned by \`make_mask_names()\`. Optionally fills
selected entries from a user-supplied named vector \`mask\`. The output
is intended for use within downstream functions (e.g., \`loglik_math\`)
that need parameters in a fixed order and with standard names: \`mu1\`
... \`mup\`, \`sigltil1\`, ..., \`sigltilp\`, \`sigrtil1\`, ...,
\`sigrtilp\`,\`ctil\`, \`pd\`, and `o_mati` for `i` ranging from 1 to
`(p^2-p)/2`.

## Usage

``` r
create_mask(mask = NULL, p = 1)
```

## Arguments

- mask:

  Named numeric vector (default \`NULL\`). Names must be a subset of
  those produced by \`make_mask_names(p)\`. Values are inserted into the
  corresponding positions; all other entries of the output are
  \`NA_real\_\`.

- p:

  A positive integer representing the number of environmental variables
  to be used in the xsdm model, i.e., `dim(env_dat)[3]` for the
  `env_dat` argument of the `loglik_math` function.

## Value

A named numeric vector of length \`num_par(p)\` with names in the
canonical order given above. Entries are initialized to \`NA_real\_\`
except for those provided in \`mask\`.

## See also

\[make_mask_names()\], \[num_par()\]

## Examples

``` r
# Empty mask for p = 2 (all NA values)
create_mask(p = 2)
#>      mu1      mu2 sigltil1 sigltil2 sigrtil1 sigrtil2     ctil       pd 
#>       NA       NA       NA       NA       NA       NA       NA       NA 
#>   o_par1 
#>       NA 

# Partially filled mask; unspecified entries remain NA
create_mask(mask = c(mu1 = 11, sigltil1 = Inf, pd = 1, ctil = -2), p = 2)
#>      mu1      mu2 sigltil1 sigltil2 sigrtil1 sigrtil2     ctil       pd 
#>       11       NA      Inf       NA       NA       NA       -2        1 
#>   o_par1 
#>       NA 

# p = 1 has no o_par entries
create_mask(mask = c(mu1 = 7, pd = 0.5), p = 1)
#>      mu1 sigltil1 sigrtil1     ctil       pd 
#>      7.0       NA       NA       NA      0.5 
```
