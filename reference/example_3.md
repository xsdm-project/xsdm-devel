# Consolidated example data for the xsdm. This is environmental data array and an occurrence presence absence vector Blarina carolinensis A named list containing all example datasets used in the package's documentation and examples.

Consolidated example data for the xsdm. This is environmental data array
and an occurrence presence absence vector Blarina carolinensis A named
list containing all example datasets used in the package's documentation
and examples.

## Usage

``` r
example_3
```

## Format

A list of 2 objects:

- `env_array`:

  A 3‑D numeric array with dimensions \`1156 (locations) × 39 (time) × 2
  (variables)\`. Contains the environmental data (bio1 and bio12)
  extracted from the rasters for all locations.

- `occ_vec`:

  An integer vector of length 1156 Binary presence/absence (0/1) for the
  same locations as \`env_array\`.

## Source

Berti et al., 2025
([doi:10.1101/2024.10.30.621023](https://doi.org/10.1101/2024.10.30.621023)
)

## Examples

``` r
# Access the list
names(example_3)
#> [1] "env_array" "occ_vec"  
```
