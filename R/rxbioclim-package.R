#' rxbioclim: Bioclimatic Variables from Monthly Climate Data
#'
#' Computes the 19 standard bioclimatic variables (BIO01-BIO19) from monthly
#' climate data following the WorldClim specification. This is an R
#' implementation of the xbioclim C++ library.
#'
#' @details
#' The package provides two levels of interface:
#'
#' - **Individual variable functions** (`bio01()` … `bio19()`): compute a
#'   single bioclimatic variable from the required monthly inputs.
#' - **Unified wrapper** (`bioclim()`): compute all 19 variables in one call,
#'   returning a named numeric vector of length 19.
#'
#' All functions accept numeric vectors of length 12 (one value per calendar
#' month) and validate their inputs before computing.
#'
#' ## Vignettes
#'
#' - `vignette("getting-started", package = "rxbioclim")` – introduction,
#'   real-world examples, and multi-pixel processing patterns.
#' - `vignette("terra-comparison", package = "rxbioclim")` – side-by-side
#'   comparison with the **terra** package.
#' - `vignette("benchmarking", package = "rxbioclim")` – block-based
#'   processing strategies and performance benchmarks.
#' - `vignette("architecture", package = "rxbioclim")` – internal design,
#'   formula reference, and extension guide for advanced users.
#'
#' @seealso
#' - [bioclim()] for the unified all-variables function.
#' - [bioclim-variables] for documentation of each individual `bioXX()` function.
#' - The upstream C++ library: <https://github.com/alrobles/xbioclim>
#' - WorldClim bioclimatic variable definitions:
#'   <https://www.worldclim.org/data/bioclim.html>
#'
#' @docType package
#' @name rxbioclim-package
"_PACKAGE"
