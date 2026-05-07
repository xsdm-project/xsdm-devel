# End-to-end MLE workflow with xsdm

## Overview

`xsdm` reconstructs a species’ fundamental ecological niche by
maximising a likelihood that combines:

- an **environmental time-series** (locations × time × variables),
- **presence / absence** records on the same locations.

This vignette walks through the complete workflow on the bundled
`example_1` fixture: build a starting design, run multi-start MLE,
inspect results, and project them spatially.

## Setup

``` r

library(xsdm)
#> Package 'xsdm' version 1.0.0
#> Type 'citation("xsdm")' for citing this R package in publications.

# The example_1 object ships with the package and provides a small,
# fully-formed environmental array + occurrence vector.
str(example_1, max.level = 1)
#> List of 12
#>  $ par_vec                 : Named num [1:9] 14.2 6.8 0.336 -0.693 0.3 ...
#>   ..- attr(*, "names")= chr [1:9] "mu1" "mu2" "sigltil1" "sigltil2" ...
#>  $ bio01                   :
#> Loading required namespace: terra
#> Formal class 'PackedSpatRaster' [package "terra"] with 3 slots
#>  $ bio12                   :Formal class 'PackedSpatRaster' [package "terra"] with 3 slots
#>  $ env_array               : num [1:2000, 1:39, 1:2] 12.97 15.34 13.95 12.88 9.77 ...
#>   ..- attr(*, "dimnames")=List of 3
#>  $ occ_df                  :Classes 'tbl_df', 'tbl' and 'data.frame':    4000 obs. of  4 variables:
#>  $ occ_vec                 : int [1:2000] 1 0 1 1 1 1 0 1 1 1 ...
#>  $ optim_par_list          :List of 6
#>  $ optim_par_vec_equivalent: Named num [1:9] 14.111 6.562 -2.252 0.296 -0.771 ...
#>   ..- attr(*, "names")= chr [1:9] "mu1" "mu2" "sigltil1" "sigltil2" ...
#>  $ optim_par_vec           : Named num [1:9] 14.1105 6.5619 -0.7706 0.0787 -2.2524 ...
#>   ..- attr(*, "names")= chr [1:9] "mu1" "mu2" "sigltil1" "sigltil2" ...
#>  $ par_list                :List of 6
#>  $ par_vec_vsp             : Named num [1:9] 14.2 6.8 0.336 -0.693 0.3 ...
#>   ..- attr(*, "names")= chr [1:9] "mu1" "mu2" "sigltil1" "sigltil2" ...
#>  $ par_table               :Classes 'tbl_df', 'tbl' and 'data.frame':    100 obs. of  9 variables:
```

## 1. Multi-start maximum-likelihood estimation

[`optimize_likelihood()`](https://xsdm-project.github.io/xsdm-devel/reference/optimize_likelihood.md)
is the user-facing entry point. It builds a Sobol’ multi-start design of
`num_starts` parameter vectors over a data-driven range, runs each start
through the C++ engine, and returns both the per-start solutions and the
best one.

``` r

fit <- optimize_likelihood(
  env_dat     = example_1$env_array,
  occ         = example_1$occ_vec,
  num_starts  = 8L,        # use 50-200 for real analyses
  parallel    = FALSE,     # set TRUE on multi-core machines
  control     = list(maxeval = 200),
  verbose     = FALSE
)
#>  ■■■■■                             12% |  ETA: 16s
#>  ■■■■■■■■■                         25% |  ETA: 15s
#>  ■■■■■■■■■■■■■■■■                  50% |  ETA: 10s
#>  ■■■■■■■■■■■■■■■■■■■■              62% |  ETA:  7s
#>  ■■■■■■■■■■■■■■■■■■■■■■■           75% |  ETA:  5s
#>  ■■■■■■■■■■■■■■■■■■■■■■■■■■■       88% |  ETA:  2s

# Best solution and its log-likelihood
fit$best$loglik
#> [1] -629.2514
fit$best$par
#>         mu1         mu2    sigltil1    sigltil2    sigrtil1    sigrtil2 
#>  14.1727105   6.7033915   0.3031849  -0.7387138   0.2909225  -2.2950287 
#>        ctil          pd      o_par1 
#> -18.2637143   2.1425111  -6.0887805
```

## 2. Habitat-suitability projection

``` r

# Compute habitat suitability from the best-fit parameters.
# Skipped at vignette-build time because terra rasters are packed.
best_bio <- math_to_bio(fit$best$par)
bio01 <- terra::unwrap(example_1$bio01) / 100
bio12 <- terra::unwrap(example_1$bio12) / 100
suit <- habitat_suitability(
  param_list  = best_bio,
  env_list    = list(bio1 = bio01, bio12 = bio12),
  return_prob = TRUE
)
```

## 3. Profile likelihood

``` r

# Profile the detection probability `pd` around the MLE.
prof <- profile_likelihood(
  profile_parameter  = "pd",
  optim_param_vector = fit$best$par,
  env_dat            = example_1$env_array,
  occ                = example_1$occ_vec,
  increment_left     = 0.2,
  increment_right    = 0.2,
  num_steps_left     = 20L,
  num_steps_right    = 20L
)
```

## See also

- `?xsdm-package` — package overview.
- [`?optimize_likelihood`](https://xsdm-project.github.io/xsdm-devel/reference/optimize_likelihood.md)
  — full argument reference.
- [`?profile_likelihood`](https://xsdm-project.github.io/xsdm-devel/reference/profile_likelihood.md)
  — uncertainty estimation.
- [`?make_mask_names`](https://xsdm-project.github.io/xsdm-devel/reference/make_mask_names.md)
  — fixing specific parameters during optimization.
