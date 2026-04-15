# rxbioclim

<!-- badges: start -->
[![R-CMD-check](https://github.com/alrobles/rxbioclim/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/alrobles/rxbioclim/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/alrobles/rxbioclim/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/alrobles/rxbioclim/actions/workflows/test-coverage.yaml)
<!-- badges: end -->

An R package for computing the 19 standard bioclimatic variables (BIO01–BIO19) from monthly climate data, following the [WorldClim](https://www.worldclim.org/data/bioclim.html) specification. This is an R implementation of the [xbioclim](https://github.com/alrobles/xbioclim) C++ library.

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("alrobles/rxbioclim")
```

## Usage

### Single-pixel (vector) interface

```r
library(rxbioclim)

# Monthly climate data (12 values, one per month)
tas    <- c(5, 7, 10, 14, 18, 22, 25, 24, 20, 15, 10, 6)
tasmax <- c(8, 10, 14, 18, 23, 28, 32, 31, 26, 19, 13, 9)
tasmin <- c(1, 3, 6, 10, 13, 17, 20, 19, 15, 10, 6, 2)
pr     <- c(60, 55, 50, 40, 30, 15, 5, 10, 25, 45, 55, 65)

# Compute all 19 bioclimatic variables at once
result <- bioclim(tas, tasmax, tasmin, pr)
print(result)

# Or compute individual variables
bio01(tas)         # Mean Annual Temperature
bio12(pr)          # Annual Precipitation
bio04(tas)         # Temperature Seasonality
bio15(pr)          # Precipitation Seasonality
```

### Raster (SpatRaster) interface

For large rasters, `bioclim_raster()` uses terra's block-loop architecture to
process data one block at a time, keeping memory use bounded regardless of
raster size. Multi-core processing within each block is supported via the
`ncores` argument.

```r
library(rxbioclim)
library(terra)

# Each SpatRaster must have exactly 12 layers (one per month)
# tas    <- rast("path/to/monthly_tas.tif")
# tasmax <- rast("path/to/monthly_tasmax.tif")
# tasmin <- rast("path/to/monthly_tasmin.tif")
# pr     <- rast("path/to/monthly_pr.tif")

# Sequential (memory-efficient block processing)
bio <- bioclim_raster(tas, tasmax, tasmin, pr)

# Write directly to file to avoid loading the full result into RAM
bio <- bioclim_raster(tas, tasmax, tasmin, pr,
                       filename = "bioclim_output.tif",
                       overwrite = TRUE)

# Multi-core: process cells within each block in parallel
bio <- bioclim_raster(tas, tasmax, tasmin, pr, ncores = 4L)

nlyr(bio)    # 19
names(bio)   # "bio01" ... "bio19"
```

### Native GDAL engine (`bioclim_engine`)

For maximum control and minimal file sizes, `bioclim_engine()` reads climate
data directly via GDAL, writes each bioclimatic variable to a **separate
single-band GeoTIFF** inside an output directory, and lets you select which
of the 19 variables to compute.

```r
library(rxbioclim)

# Compute all 19 variables — one file each in a directory
result <- bioclim_engine(
  "tas.tif", "tasmax.tif", "tasmin.tif", "pr.tif",
  output = "bioclim_output/",
  overwrite = TRUE
)
list.files("bioclim_output/")
# "bio01.tif" "bio02.tif" ... "bio19.tif"

# Compute only BIO01 (mean annual temp) and BIO12 (annual precip)
result <- bioclim_engine(
  "tas.tif", "tasmax.tif", "tasmin.tif", "pr.tif",
  output    = "bioclim_subset/",
  variables = c(1L, 12L),
  overwrite = TRUE
)
names(result)  # "bio01" "bio12"
```

## Bioclimatic Variables

| Variable | Description |
|----------|-------------|
| BIO01 | Mean Annual Temperature |
| BIO02 | Mean Diurnal Range |
| BIO03 | Isothermality (100 × BIO02 / BIO07) |
| BIO04 | Temperature Seasonality (100 × population SD) |
| BIO05 | Max Temperature of Warmest Month |
| BIO06 | Min Temperature of Coldest Month |
| BIO07 | Temperature Annual Range (BIO05 − BIO06) |
| BIO08 | Mean Temperature of Wettest Quarter |
| BIO09 | Mean Temperature of Driest Quarter |
| BIO10 | Mean Temperature of Warmest Quarter |
| BIO11 | Mean Temperature of Coldest Quarter |
| BIO12 | Annual Precipitation |
| BIO13 | Precipitation of Wettest Month |
| BIO14 | Precipitation of Driest Month |
| BIO15 | Precipitation Seasonality (CV) |
| BIO16 | Precipitation of Wettest Quarter |
| BIO17 | Precipitation of Driest Quarter |
| BIO18 | Precipitation of Warmest Quarter |
| BIO19 | Precipitation of Coldest Quarter |

## Documentation

Comprehensive vignettes are available after installing the package:

```r
vignette("getting-started",  package = "rxbioclim")  # Introduction & real-world examples
vignette("terra-comparison", package = "rxbioclim")  # rxbioclim vs terra
vignette("benchmarking",     package = "rxbioclim")  # Block-based performance
vignette("architecture",     package = "rxbioclim")  # Design & internals
```

## License

MIT
