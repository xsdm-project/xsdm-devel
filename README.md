
<!-- README.md is generated from README.Rmd. Please edit that file -->

# xsdm

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/xsdm-project/xsdm-devel/badge.svg)](https://app.codecov.io/gh/xsdm-project/xsdm-devel)
–\> <!-- badges: end -->

The goal of xsdmMle is to fit

## Installation

You can install the development version of xsdmMle from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("alrobles/xsdmMle")
```

## Example

This is a basic example to fit the distribution of a species:

``` r
library(xsdmMle)
## 01 We read the rasters for our example

bio1_ts <- terra::unwrap(example_1_bio01)
bio12_ts <- terra::unwrap(example_1_bio12)

# 02 We transform the raster of precipitation to a similar scale of temperature

bio1_ts <- bio1_ts / 100
bio12_ts <- bio12_ts / 100

# 03 We create a list of raster time series
envData <- list(bio1 = bio1_ts, bio12 = bio12_ts)

# 04 we create a data array providing the list of environmental data and the occurence
# data frame (sp_virtual_example)
envdat <- env_data_array(envData, occ = example_1_occurrence_df)
```

Following this we fit using the function optimize_likelihood:

``` r
# 01 We store the presence absence vector
occ <- example_1_occurrence_df$presence

# 03 Then we provide the environmental data array and the occ vector in the optim_mll.
# We can set parameters here like the number of initial points, or the flag if
# we want to run in parallel
optim_df <- optimize_likelihood(envdat, occ, num_starts = 5, parallel = TRUE)
```

\`\`\`

Then we got a data frame with the maximum likelihood estimation for each
set of initial starting points in the optimization. We have this on the
biological scale of parameters (Is the unconstrained space of parameters
where the optimizer found the maximum of the function).
