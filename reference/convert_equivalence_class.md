# Converts a set of parameters to other representatives of the same equivalence class

Model parameters in the biological scale are only determined up to an
equivalence class. This function converts a set of parameters to another
set of equivalent parameters.

## Usage

``` r
convert_equivalence_class(p, flip, perm)
```

## Arguments

- p:

  Named list with entries mu, sigltil, sigrtil, ctil, pd, and o_mat

- flip:

  Vector of binaries corresponding to which columns of o_mat are to have
  their sign change (which a concomitant switch of the corresponding
  entries of sigltil and sigrtil). Length must equal the number of
  columns of o_mat.

- perm:

  Permutation to be applied to the columns of o_mat.

## Value

List with entries o_mat, sigltil, and sigrtil

## Examples

``` r
convert_equivalence_class(
  p = example_1$optim_par_list,
  flip = c(1, 0),
  perm = c(1, 2)
)
#> $o_mat
#>            [,1]       [,2]
#> [1,] -0.9758770 -0.2183212
#> [2,] -0.2183212  0.9758770
#> 
#> $sigltil
#> [1] 0.5354781 0.1498244
#> 
#> $sigrtil
#> [1] 0.2624513 0.8570480
#> 
```
