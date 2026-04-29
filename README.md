
<!-- README.md is generated from README.Rmd. Please edit that file -->

# xsdm

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/xsdm-project/xsdm-devel/badge.svg)](https://app.codecov.io/gh/xsdm-project/xsdm-devel)
[![R-CMD-check](https://github.com/xsdm-project/xsdm-devel/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/xsdm-project/xsdm-devel/actions/workflows/R-CMD-check.yaml)
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

**xsdm** is an R package that integrates concepts of *stochastic
demography* into species distribution modelling (SDM). Instead of
treating environmental conditions as a static snapshot, xsdm uses
**multi-year environmental time-series** together with species
presence/absence records to:

- Reconstruct a species’ **fundamental ecological niche** via
  maximum-likelihood estimation.
- Account for **inter-annual climate variability** when estimating niche
  breadth and position.
- Project the species’ **potential geographic range** under current or
  future climate scenarios.

The statistical underpinning is described in:

> Berti, E., Robles Fernández, A.L., Rosenbaum, B., Peterson, T.A.,
> Soberón, J., & Reuman, D.C. (2025). *The impacts of climate
> variability on the niche concept and distributions of species*.
> bioRxiv. <https://doi.org/10.1101/2024.10.30.621023>

------------------------------------------------------------------------

## Installation

### CRAN (coming soon)

``` r
install.packages("xsdm")   # not yet on CRAN — see development install below
```

### Development version

Install the latest version directly from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("xsdm-project/xsdm-devel")
```

xsdm requires **R ≥ 4.1** and links to C++ via Rcpp/RcppParallel, so a
working compiler (e.g. Rtools on Windows, Xcode CLT on macOS) is needed.

------------------------------------------------------------------------

## Quick start

The package ships with a built-in dataset (`example_1`) containing all
objects needed to run a complete workflow without loading external
files.

``` r
library(xsdm)

# The example_1 dataset is a named list bundling all example objects.
names(example_1)
#> [1] "par_vec"                  "bio01"
#> [3] "bio12"                    "env_array"
#> [5] "occ_df"                   "occ_vec"
#> [7] "optim_par_list"           "optim_par_vec"
#> [9] "optim_par_vec_equivalent" "par_list"
#> [11] "par_vec_vsp"             "par_table"

# Fit the model from the pre-built environmental array
result <- optimize_likelihood(
  env_dat    = example_1$env_array,
  occ        = example_1$occ_vec,
  num_starts = 10L
)

# Inspect the best solution
result$best$loglik
result$best$convergence

# Convert the best math-scale parameters to biological scale
best_bio <- math_to_bio(result$best$par)
str(best_bio)
```

------------------------------------------------------------------------

## Full end-to-end workflow

### 1 · Load and prepare environmental time series

xsdm expects **bioclimatic time-series rasters** — one `SpatRaster` per
variable with each layer representing one year (or time step). The
built-in data cover southern New Mexico, USA, over 39 years (1980–2018)
using [CHELSA v2.1](https://www.chelsa-climate.org/) bio1 (mean annual
temperature) and bio12 (annual precipitation).

``` r
library(xsdm)
library(terra)

# The rasters are stored as packed SpatRasters to minimise package size.
# Unpack them before use.
bio1_ts  <- terra::unwrap(example_1$bio01)
bio12_ts <- terra::unwrap(example_1$bio12)

# Scale to interpretable units: CHELSA bio1 is in 0.1 °C → convert to °C,
# and bio12 is in kg/m² → scale to match temperature magnitude.
bio1_ts  <- bio1_ts  / 100
bio12_ts <- bio12_ts / 100
```

### 2 · Build the environmental data array

`env_data_array()` extracts and stacks environmental values at the
locations given in a presence/absence data frame, returning a 3-D array
`(locations × time × variables)`.

``` r
# Named list of raster time series (names become the variable labels)
env_data <- list(bio1 = bio1_ts, bio12 = bio12_ts)

# example_1$occ_df has columns: name, x (lon), y (lat), presence (0/1)
env_dat <- env_data_array(env_data, occ = example_1$occ_df)
dim(env_dat)
#> [1] 1000   39    2    (locations × time-steps × variables)

occ <- example_1$occ_df$presence
```

### 3 · Fit the model

`optimize_likelihood()` runs multiple optimizations from Latin-hypercube
starting points (Sobolʼ design) and returns every solution sorted by
decreasing log-likelihood. Internally, optimization uses
`ucminfcpp::ucminf_xptr()` with a compiled objective pointer.

``` r
result <- optimize_likelihood(
  env_dat    = env_dat,
  occ        = occ,
  num_starts = 20L,       # increase for a real analysis
  parallel   = TRUE,       # set TRUE to use future/furrr parallelism
  verbose    = TRUE
)

# Solutions data frame: one row per starting point
head(result$solutions[, c("start_id", "loglik", "convergence")])

# Best solution
result$best$loglik
```

### 4 · Interpret the fitted parameters

Convert the best math-scale parameter vector to the biologically
interpretable scale and plot the inferred log growth–environment
function:

``` r
best_bio <- math_to_bio(result$best$par)

# key biological parameters:
# mu       – optimal environmental values
# sigltil  – left-side niche widths  (Inf = no left boundary)
# sigrtil  – right-side niche widths (Inf = no right boundary)
# pd       – maximum probability of detection (0–1)
# o_mat    – rotation matrix (correlations between environmental axes)
str(best_bio)

# Visualise the niche contours (two-panel: presences vs non-detections)
interpret_parameters(
  best_bio,
  plot_indices = c(1, 2),
  env_dat      = env_dat,
  occ          = occ
)
```

### 5 · Project the range (virtual species probability map)

`vsp()` uses the fitted parameter list and the full-grid rasters to
return per-cell detection probabilities:

``` r
# Use the full raster extent (not just occurrence points)
env_data_full <- list(bio1 = bio1_ts, bio12 = bio12_ts)

# Returns a tibble with x, y, probs by default
prob_tbl  <- vsp(env_data_full, best_bio)
head(prob_tbl)

# Or as a SpatRaster for mapping
prob_rast <- vsp(env_data_full, best_bio, return_raster = TRUE)
terra::plot(prob_rast, main = "Detection probability")
```

------------------------------------------------------------------------

## Key functions

| Function | Purpose |
|----|----|
| `env_data_array()` | Build a `(locations × time × variables)` array from raster time series and an occurrence table |
| `optimize_likelihood()` | Multi-start MLE fitting; returns solutions sorted by log-likelihood |
| `loglik_math()` | Evaluate the log-likelihood at any math-scale parameter vector |
| `math_to_bio()` | Convert math-scale vector → biological-scale parameter list |
| `bio_to_math()` | Convert biological-scale parameter list → math-scale vector |
| `start_parms()` | Generate Latin-hypercube starting points from presence-only data |
| `profile_likelihood()` | Profile one parameter while re-optimising over the rest |
| `vsp()` | Produce a spatial probability-of-detection map |
| `interpret_parameters()` | Diagnostic plots of the niche shape |
| `dist_between_params()` | Distance between two parameter sets (Hungarian algorithm, equivalence-class aware) |

------------------------------------------------------------------------

## Uncertainty quantification: profile likelihood

`profile_likelihood()` fixes a target parameter at a grid of values,
re-optimises all others at each grid point, and returns the profile
log-likelihood together with a likelihood-ratio confidence threshold.

``` r
prof <- profile_likelihood(
  profile_parameter  = "mu1",       # parameter to profile (math scale)
  increment_left     = 0.2,
  increment_right    = 0.2,
  num_steps_left     = 20L,
  num_steps_right    = 20L,
  alpha              = 0.95,        # 95 % LR confidence level
  optim_param_vector = result$best$par,
  env_dat            = env_dat,
  occ                = occ,
  verbose            = FALSE
)

# Profile data frame: param, value_math, loglik, convergence
head(prof$profile[, c("param", "value_math", "loglik")])
prof$threshold    # LR cut-off for the 95 % CI
prof$found_better # TRUE if profiling found a better point than the MLE
```

------------------------------------------------------------------------

## Parameter comparison

Because the xsdm model has a partial non-identifiability (parameters are
defined only up to an equivalence class), ordinary Euclidean distance is
not appropriate. `dist_between_params()` uses the Hungarian algorithm to
find the equivalence-class representative of `p1` closest to `p2`:

``` r
# Both representations belong to the same equivalence class → distance ≈ 0
dist_between_params(
  p1 = example_1$optim_par_vec,
  p2 = example_1$optim_par_vec_equivalent
)
```

------------------------------------------------------------------------

## For developers and contributors

### Install from source

``` r
# Clone the repository and install locally
install.packages("devtools")
devtools::install_github("xsdm-project/xsdm-devel")
```

### Run the test suite

``` r
devtools::test()
```

Tests are written with **testthat** (≥ 3.0) and cover all exported
functions. Code coverage is tracked via
[Codecov](https://app.codecov.io/gh/xsdm-project/xsdm-devel).

### Manual visual tests

Some diagnostics for `interpret_parameters()` require human inspection.
See `tests/manual/README.md` for instructions on running those scripts.

### Reporting issues

Please file bugs and feature requests on the [GitHub issue
tracker](https://github.com/xsdm-project/xsdm-devel/issues). When
reporting a bug, include the output of `sessionInfo()` and a minimal
reproducible example.

### Extending xsdm

The internal log-likelihood engine is implemented in C++ via **Rcpp**
and **RcppParallel**. For statistically or computationally sophisticated
users, the raw likelihood function (`loglik_math()`) and
parameter-transformation utilities (`math_to_bio()`, `bio_to_math()`,
`make_mask_names()`) are fully exported, enabling custom optimization
workflows.

------------------------------------------------------------------------

## Citation

If you use xsdm in published work, please cite:

    Berti, E., Robles Fernández, A.L., Rosenbaum, B., Peterson, T.A.,
    Soberón, J., & Reuman, D.C. (2025). The impacts of climate variability
    on the niche concept and distributions of species. bioRxiv.
    https://doi.org/10.1101/2024.10.30.621023

BibTeX:

``` bibtex
@article{bertiXSDM1,
  title  = {The impacts of climate variability on the niche concept and
            distributions of species},
  author = {Berti, E. and Fern\'{a}ndez, ALR and Rosenbaum, B and
            Peterson, TA and Sober\'{o}n, J and Reuman, DC},
  journal = {bioRxiv},
  doi    = {10.1101/2024.10.30.621023},
  year   = {2025},
  url    = {https://doi.org/10.1101/2024.10.30.621023}
}
```

------------------------------------------------------------------------

## License

xsdm is released under the [GNU General Public License v3 or
later](LICENSE.md).
