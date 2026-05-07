# Tiled habitat-suitability map from environmental raster stacks

Evaluates the log detection probability (or its exponential, the
probability of detection) for every cell of a list of multi-layer
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
objects, processing the inputs in memory-bounded blocks so that
arbitrarily large grids can be handled without loading the entire
dataset into R memory. Each block is forwarded to
[`log_prob_detect_cpp`](https://xsdm-project.github.io/xsdm-devel/reference/log_prob_detect_cpp.md),
the xtensor-backed C++ kernel that consolidates the
`like_neg_ltsgr() -> like_ltsg()` call chain.

## Usage

``` r
habitat_suitability(
  param_list,
  env_list,
  output = "",
  overwrite = FALSE,
  return_prob = TRUE,
  threads = 0L,
  wopt = list()
)
```

## Arguments

- param_list:

  A named list of biological-scale parameters. Must contain `mu`,
  `sigltil`, `sigrtil`, `o_mat`, `ctil` and `pd`. See
  [`log_prob_detect`](https://xsdm-project.github.io/xsdm-devel/reference/log_prob_detect.md)
  for details of each element.

- env_list:

  A list of
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  objects, one per environmental variable. Each raster must have the
  same number of layers (time steps) and identical spatial geometry
  (extent, resolution, CRS). Minimum length 1.

- output:

  Character scalar. File path for the output GeoTIFF. The empty string
  `""` (default) creates an in-memory `SpatRaster`.

- overwrite:

  Logical scalar. If `TRUE`, an existing file at `output` is
  overwritten. Default `FALSE`.

- return_prob:

  Logical scalar. If `TRUE` (default), the output cell values are
  probabilities of detection (range \\\[0, 1\]\\). If `FALSE`, the cell
  values are log-probabilities (range \\(-\infty, 0\]\\).

- threads:

  Integer scalar. Number of parallel threads forwarded to
  [`log_prob_detect_cpp`](https://xsdm-project.github.io/xsdm-devel/reference/log_prob_detect_cpp.md).
  Use `0` (default) to let RcppParallel pick the number of threads
  automatically.

- wopt:

  List. Additional write options forwarded to
  [`writeStart`](https://rspatial.github.io/terra/reference/readwrite.html).
  Default [`list()`](https://rdrr.io/r/base/list.html).

## Value

A `SpatRaster` with one layer named either `"habitat_suitability"` (when
`return_prob = TRUE`) or `"log_prob_detect"` (when
`return_prob = FALSE`). The raster is returned invisibly when
`output != ""`.

## Details

Internally the function uses terra's streaming block-loop API:

1.  [`readStart`](https://rspatial.github.io/terra/reference/readwrite.html)
    is called on every raster in `env_list`.

2.  [`writeStart`](https://rspatial.github.io/terra/reference/readwrite.html)
    is called on the output raster, which returns a block schedule
    chosen by terra's memory manager.

3.  For each block,
    [`readValues`](https://rspatial.github.io/terra/reference/readwrite.html)
    reads a horizontal strip from every input raster into a matrix; the
    strips are packed into a flat column-major vector and passed to
    [`log_prob_detect_cpp`](https://xsdm-project.github.io/xsdm-devel/reference/log_prob_detect_cpp.md).
    Cells that are NA in any variable or time step are masked out and
    re-inserted as NA in the output.

4.  [`writeValues`](https://rspatial.github.io/terra/reference/readwrite.html)
    writes the per-cell results.

5.  [`readStop`](https://rspatial.github.io/terra/reference/readwrite.html)
    and
    [`writeStop`](https://rspatial.github.io/terra/reference/readwrite.html)
    are called via [`on.exit`](https://rdrr.io/r/base/on.exit.html) to
    ensure file handles are released even if an error occurs.

At most one block of pixels is held in R memory at any time, making the
function suitable for continental or global rasters.

## See also

[`log_prob_detect_cpp`](https://xsdm-project.github.io/xsdm-devel/reference/log_prob_detect_cpp.md),
[`log_prob_detect`](https://xsdm-project.github.io/xsdm-devel/reference/log_prob_detect.md),
[`vsp`](https://xsdm-project.github.io/xsdm-devel/reference/vsp.md),
[`writeStart`](https://rspatial.github.io/terra/reference/readwrite.html)

## Examples

``` r
# \donttest{
data("example_1", package = "xsdm")
bio01 <- terra::unwrap(example_1$bio01) / 100
bio12 <- terra::unwrap(example_1$bio12) / 100
env_list <- list(bio01 = bio01, bio12 = bio12)
suit <- habitat_suitability(
  param_list  = example_1$par_list,
  env_list    = env_list,
  return_prob = TRUE
)
suit
#> class       : SpatRaster
#> size        : 128, 123, 1  (nrow, ncol, nlyr)
#> resolution  : 5000, 5000  (x, y)
#> extent      : -1231223, -616223, 989721.9, 1629722  (xmin, xmax, ymin, ymax)
#> coord. ref. : +proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs
#> source(s)   : memory
#> name        : habitat_suitability
#> min value   :                   0
#> max value   :            0.864304
# }
```
