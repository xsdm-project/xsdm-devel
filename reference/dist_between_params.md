# Distance in parameter space between two sets of parameters

Computes a distance in parameter space between two parameter sets of the
xsdm model. The `o_mat`, `sigltil`, and `sigrtil` parameters are only
determined up to an equivalence class; this function returns the minimum
distance over all equivalence-class representatives of `p1`, using the
Hungarian (Kuhn–Munkres) linear sum assignment algorithm to avoid
enumerating every permutation and sign flip. Distance is measured in
sum-squared-errors on the biological scale, except that `sigltil` and
`sigrtil` are inverted before comparison (that is the scale on which
distance is most meaningful for those parameters).

## Usage

``` r
dist_between_params(p1, p2, mask = NULL, give_closest_rep = FALSE)
```

## Arguments

- p1:

  First set of parameters. May be math-scale (a named numeric vector
  whose names complement `mask`) or biological-scale (a named list with
  entries `mu`, `sigltil`, `sigrtil`, `ctil`, `pd`, `o_mat`).

- p2:

  Second set of parameters; same format options as `p1`.

- mask:

  Same format as the `mask` argument to
  [`loglik_math`](https://xsdm-project.github.io/xsdm-devel/reference/loglik_math.md)
  and
  [`start_parms`](https://xsdm-project.github.io/xsdm-devel/reference/start_parms.md).
  Ignored if both `p1` and `p2` are on the biological scale. Otherwise
  the names of `mask` must exactly complement the names of whichever of
  `p1` / `p2` is on the math scale.

- give_closest_rep:

  If `TRUE`, also returns the member of the equivalence class of `p1`
  that attains the minimum distance (biological scale). Default `FALSE`.

## Value

If `give_closest_rep` is `FALSE`, a single number: the distance.
Otherwise a list with entries `distance` and `representative`.

## Details

All numerics besides the linear sum assignment are computed in R. The
assignment problem itself is solved by an in-package C++ implementation
of the classical O(n^3) Hungarian algorithm, exposed (unexported) as
`xsdm:::.solve_lsap_cpp`. An R-level alternative is
[`clue::solve_LSAP`](https://rdrr.io/pkg/clue/man/solve_LSAP.html); the
two are compared in `tests/testthat/test-solve_lsap_cpp.R`.

## References

H. W. Kuhn (1955). The Hungarian Method for the Assignment Problem.
*Naval Research Logistics Quarterly* 2(1-2), 83–97.

J. Munkres (1957). Algorithms for the Assignment and Transportation
Problems. *Journal of the SIAM* 5(1), 32–38.

R. Jonker and A. Volgenant (1987). A Shortest Augmenting Path Algorithm
for Dense and Sparse Linear Assignment Problems. *Computing* 38,
325–340.

K. Hornik (2005). A CLUE for CLUster Ensembles. *Journal of Statistical
Software* 14(12). (See also the clue package on CRAN for an alternative
R-level LSAP implementation.)

## Examples

``` r
# Using lists on the biological scale
par_list <- math_to_bio(example_1$optim_par_vec)
par_list_equivalent <- math_to_bio(example_1$optim_par_vec_equivalent)
dist_between_params(
  p1 = par_list,
  p2 = par_list_equivalent
)
#> [1] 5.150436e-15

# Using vectors on the math scale
dist_between_params(
  p1 = example_1$optim_par_vec,
  p2 = example_1$optim_par_vec_equivalent
)
#> [1] 5.150436e-15
```
