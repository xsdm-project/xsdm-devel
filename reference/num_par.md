# Get the number of parameters of the main xsdm model given the number of environmental variables to be considered

Get the number of parameters of the main xsdm model given the number of
environmental variables to be considered

## Usage

``` r
num_par(p)
```

## Arguments

- p:

  A positive integer representing the number of environmental variables
  to be used in the xsdm model, i.e., `dim(env_dat)[3]` for the
  `env_dat` argument of the `loglik_math` function

## Value

An integer with the number of parameters

## Details

For instance, in the \`p=1\` case, the xsdm model parameters are \`mu\`,
\`sigltil\`, and \`sigrtil\` (which are scalars in the \`p=1\` case);
\`ctil\`, and \`pd\` (which are scalars for any value \`p\`). That makes
5 parameters, so this function returns 5. In the \`p=2\` case, the
parameters are \`mu\`, \`sigltil\`, and \`sigrtil\` (each of which is
now a length-2 vector); \`ctil\`, and \`pd\` (again scalars); and the
single parameter pertaining to \`o_mat\`; for a total of 9.

## Examples

``` r
num_par(2)
#> [1] 9
```
