# Converts parameters from the biological scale to the math (unconstrained) scale

Given a list with biological-scale parameters (\`mu\`, \`sigltil\`,
\`sigrtil\`, \`ctil\`, \`pd\`, \`o_mat\`), returns a named numeric
vector on the math scale, in the canonical order produced by
\`make_mask_names(p)\`: \`mu1..mup\`, \`sigltil1..p\`, \`sigrtil1..p\`,
\`o_par1..q\`, \`ctil\`, \`pd\`, where \`q = p\*(p-1)/2\` and \`p =
length(mu) = nrow(o_mat)\`.

## Usage

``` r
bio_to_math(parms_bio)
```

## Arguments

- parms_bio:

  A named list with entries: \`mu\`, \`sigltil\`, \`sigrtil\`, \`ctil\`,
  \`pd\`, \`o_mat\`.

## Value

A named numeric vector on the math scale, ordered per
\`make_mask_names(p)\`.

## Details

Transformations: - \`mu\` : identity - \`sigltil\` : \`log()\` -
\`sigrtil\` : \`log()\` - \`ctil\` : identity - \`pd\` : \`logit()\` -
\`o_mat\` : lower-triangular parameters recovered via
\`extract_orthogonal_matrix_parameters()\` (see Details)

The \`o_mat\` entries are mapped to a vector via the principal matrix
logarithm, i.e., one of the (skew-symmetric) matrices S such that
\`o_mat = expm(S)\`. The math-scale parameters \`o_par\` are then the
strictly lower-triangular elements of \`S\`. For \`p = 1\`, there are no
\`o_par\` entries. Note, however, that the principal logarithm is not
defined for all special orthogonal matrices (even though all such
matrices are in the image of the matrix exponential), so this function
may fail for some valid \`o_mat\` inputs.

## See also

\[math_to_bio()\], \[make_mask_names()\], \[build_orthogonal_matrix()\]

## Examples

``` r
## --- p = 1 (no o_par entries) ---
mu1 <- 10
sigltil1 <- 1.2
sigrtil1 <- 0.8
bio_parameters <- list(
  mu      = c(mu1),
  sigltil = c(sigltil1),
  sigrtil = c(sigrtil1),
  ctil    = 0.3,
  pd      = 0.85,
  o_mat   = matrix(1, 1, 1) # 1x1 orthogonal
)
math1 <- bio_to_math(bio_parameters)
# Canonical names
names(math1)
#> [1] "mu1"      "sigltil1" "sigrtil1" "ctil"     "pd"      
# Back to biological scale
math_parameters <- math_to_bio(math1)
all.equal(math_parameters$mu, bio_parameters$mu)
#> [1] TRUE
all.equal(math_parameters$sigltil, bio_parameters$sigltil)
#> [1] TRUE
all.equal(math_parameters$sigrtil, bio_parameters$sigrtil)
#> [1] TRUE
all.equal(math_parameters$ctil, bio_parameters$ctil)
#> [1] TRUE
all.equal(math_parameters$pd, bio_parameters$pd)
#> [1] TRUE

## --- p = 2 (includes one o_par) ---
mu2 <- c(11, 5)
sigltil2 <- c(1.1, 1.5)
sigrtil2 <- c(1.4, 1.3)
ctil2 <- -0.2
pd2 <- 0.9
o_par2 <- 0.25
O2 <- build_orthogonal_matrix(o_par2)
bio_parameters_2d <- list(
  mu      = mu2,
  sigltil = sigltil2,
  sigrtil = sigrtil2,
  ctil    = ctil2,
  pd      = pd2,
  o_mat   = O2
)
math_parameters_2d <- bio_to_math(bio_parameters_2d)
# check canonical name order produced by make_mask_names(2)
identical(names(math_parameters_2d), names(make_mask_names(2)))
#> [1] TRUE
```
