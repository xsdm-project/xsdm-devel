# Extract a math-scale real-parameter vector corresponding to a special orthogonal matrix

Computes the principal matrix logarithm of a special orthogonal matrix
\`o_mat\`, then returns the strictly lower-triangular entries of the
resulting skew-symmetric matrix. This is a partial inverse of
\`build_orthogonal_matrix\`, up to the periodicity of the exponential
map.

## Usage

``` r
extract_orthogonal_matrix_parameters(o_mat)
```

## Arguments

- o_mat:

  A \\k \times k\\ special orthogonal matrix (`t(o_mat) %*% o_mat = I`
  and `det(o_mat) = 1`).

## Value

A numeric vector of length \\k(k-1)/2\\ containing the strictly
lower-triangular entries of the skew-symmetric generator. For \\k=1\\,
returns `NULL` (the identity matrix).

## Details

The matrix exponential \\\exp: \mathfrak{so}(k) \to SO(k)\\ is
surjective but not injective: different skew-symmetric matrices can
exponentiate to the same orthogonal matrix. This function uses the
\*\*principal matrix logarithm\*\* as implemented in \`expm::logm\`.
Consequently, it may fail (or produce complex results) for matrices that
have eigenvalues equal to \\-1\\ (i.e., rotations by \\\pi\\). Such
matrices lie on the cut locus of the exponential map and do not possess
a unique real logarithm. If you encounter this, consider perturbing the
matrix slightly away from the problematic rotation.

## Examples

``` r
o_par2 <- 0.25
O2 <- build_orthogonal_matrix(o_par2)
extract_orthogonal_matrix_parameters(O2)
#> [1] 0.25
```
