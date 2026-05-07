# Build an orthogonal matrix from a real-parameter vector

Constructs a \\k \times k\\ orthogonal matrix \\O\\ by exponentiating a
skew-symmetric matrix \\S\\ built by assigning its lower-triangular
entries from the input vector. Specifically, the function sets
\\S\_{ij}\\ (for \\i\>j\\) from \`entries\`, mirrors it to enforce \\S =
L - L^\top\\, and returns \`expm::expm(S)\`, which is guaranteed
orthogonal because \\\exp(S)\\ is orthogonal whenever \\S\\ is real
skew-symmetric.

## Usage

``` r
build_orthogonal_matrix(entries)
```

## Arguments

- entries:

  A numeric vector (possibly \`NULL\`). If \`NULL\`, returns the \`1 x
  1\` identity. Otherwise, its length must be \\n = k(k-1)/2\\ for some
  integer \\k \ge 2\\, supplying the strictly lower-triangular entries
  of a skew-symmetric generator.

## Value

A \`k x k\` orthogonal matrix. If \`entries\` is \`NULL\`, returns
\`matrix(1, 1, 1)\` (identity).

## Details

The dimension \`k\` is inferred from \`length(entries)\` via the
relation \\n = k(k-1)/2\\, so \`length(entries)\` must equal a
triangular number.

\- Dimension inference uses \\k = \frac{1 + \sqrt{1 + 8n}}{2}\\ where
\\n = \text{length(entries)}\\. If \\k\\ is not an integer, the input is
invalid and an error is thrown. - Orthogonality follows from the fact
that \\S^\top = -S\\ implies \\\exp(S)^\top \exp(S) = I\\. - Note the
function actually returns a special orthogonal matrix, i.e., the
determinant is +1.

## Examples

``` r
# 1x1 identity (NULL input)
build_orthogonal_matrix(NULL)
#>      [,1]
#> [1,]    1

# 2x2 orthogonal matrix from one parameter
O2 <- build_orthogonal_matrix(0.0)
all.equal(t(O2) %*% O2, diag(2)) # should be TRUE
#> [1] TRUE

# 3x3 example: length(entries) = 3 (= 3*2/2), so k = 3
O3 <- build_orthogonal_matrix(c(0.1, -0.2, 0.3))
all.equal(t(O3) %*% O3, diag(3), tolerance = 1e-10)
#> [1] TRUE
```
