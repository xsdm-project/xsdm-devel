# Generate a virtual species probability map with presence/absence sampling

Creates a virtual species probability-of-detection map based on
environmental time-series data and a set of species-specific parameters,
then samples presence/absence points based on a user-defined probability
threshold.

## Usage

``` r
vsp(param_list, env_data, size_presence, size_absence, threshold = 0.5)
```

## Arguments

- param_list:

  A named list of biological‑scale parameters required by
  \`log_prob_detect()\`. Must include \`mu\`, \`sigltil\`, \`sigrtil\`,
  \`ctil\`, \`pd\`, and \`o_mat\`. Values like \`sigltil\`/\`sigrtil\`
  can be \`Inf\`.

- env_data:

  A named list of time‑series raster objects (e.g., from the \`terra\`
  package). Each element must be a \`SpatRaster\` with the same geometry
  and number of layers.

- size_presence:

  Integer. Number of sample points to draw from cells where the
  detection probability \*\*exceeds\*\* \`threshold\`.

- size_absence:

  Integer. Number of sample points to draw from cells where the
  detection probability is \*\*less than or equal to\*\* \`threshold\`.

- threshold:

  Numeric in \`\[0, 1\]\`. Probability cutoff used to distinguish
  presence vs. absence sampling areas. Default \`0.5\`.

## Value

A tibble with columns \`lon\`, \`lat\`, \`presence\` (0/1), where each
row corresponds to a sampled point. The presence/absence is drawn from a
binomial distribution using the habitat suitability value as the success
probability.

## Details

Internally the function:

1.  Computes a habitat suitability raster using
    \`habitat_suitability()\`.

2.  Splits the raster into two layers based on \`threshold\`: cells with
    prob \> threshold (presence pool) and ≤ threshold (absence pool).

3.  Samples \`size_presence\` and \`size_absence\` points from each pool
    (without replacement), with probabilities proportional to the
    suitability value.

4.  Generates a binomial outcome for each sampled point using its
    suitability as the probability of success.

## See also

\[habitat_suitability()\], \[log_prob_detect()\],
\[terra::spatSample()\]

## Examples

``` r
# \donttest{
data("example_1", package = "xsdm")
bio1_ts  <- terra::unwrap(example_1$bio01) / 100
bio12_ts <- terra::unwrap(example_1$bio12) / 100
env_data <- list(bio1 = bio1_ts, bio12 = bio12_ts)

vsp(
  param_list    = example_1$par_list,
  env_data      = env_data,
  size_presence = 100,
  size_absence  = 100,
  threshold     = 0.7
)
#> # A tibble: 200 × 3
#>         lon      lat presence
#>       <dbl>    <dbl>    <int>
#>  1  -873723 1417222.        1
#>  2 -1083723 1427222.        1
#>  3 -1113723 1392222.        1
#>  4 -1083723 1252222.        0
#>  5 -1138723 1397222.        1
#>  6  -748723 1627222.        0
#>  7  -758723 1552222.        1
#>  8  -853723 1532222.        1
#>  9 -1158723 1287222.        1
#> 10  -793723 1532222.        1
#> # ℹ 190 more rows
# }
```
