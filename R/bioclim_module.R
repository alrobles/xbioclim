# Load the "bioclim_mod" Rcpp module after the package DLL is available.
.onLoad <- function(libname, pkgname) {
  Rcpp::loadModule("bioclim_mod", TRUE)
}

# Suppress R CMD check NOTE: 'no visible binding for global variable ClimateBlock'.
# ClimateBlock is injected into the package namespace by loadModule() in .onLoad().
utils::globalVariables("ClimateBlock")

#' C++ ClimateBlock class for batch bioclim computation
#'
#' An Rcpp module class exposing the C++ \code{ClimateBlock} implementation of
#' the xbioclim library. Accepts four \code{n_pixels x 12} matrices of monthly
#' climate data and computes all 19 bioclimatic variables for each pixel via a
#' compiled C++ back-end.
#'
#' \code{ClimateBlock} is loaded into the package namespace when the package
#' is attached (via \code{\link[Rcpp]{loadModule}} in \code{.onLoad}).
#'
#' @section Constructor:
#' \code{new(ClimateBlock, tas, tasmax, tasmin, pr)}
#'
#' \describe{
#'   \item{\code{tas}}{Numeric matrix of dimensions \code{n_pixels x 12}:
#'     monthly mean temperature for each pixel.}
#'   \item{\code{tasmax}}{Numeric matrix of dimensions \code{n_pixels x 12}:
#'     monthly maximum temperature for each pixel.}
#'   \item{\code{tasmin}}{Numeric matrix of dimensions \code{n_pixels x 12}:
#'     monthly minimum temperature for each pixel.}
#'   \item{\code{pr}}{Numeric matrix of dimensions \code{n_pixels x 12}:
#'     monthly precipitation for each pixel.}
#' }
#'
#' All four matrices must have exactly 12 columns and the same number of rows.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{n_pixels()}}{Returns the number of pixels (integer).}
#'   \item{\code{compute()}}{Computes the 19 bioclimatic variables for all
#'     pixels and returns an \code{n_pixels x 19} numeric matrix with columns
#'     named \code{bio01} through \code{bio19}.}
#' }
#'
#' @name ClimateBlock
#' @export ClimateBlock
#' @examples
#' tas    <- matrix(rep(1:12, 3), nrow = 3, byrow = TRUE)
#' tasmax <- tas + 1
#' tasmin <- tas - 1
#' pr     <- matrix(rep(1:12, 3), nrow = 3, byrow = TRUE)
#' block  <- new(ClimateBlock, tas, tasmax, tasmin, pr)
#' block$n_pixels()   # 3
#' result <- block$compute()   # 3 x 19 matrix
NULL

#' Create a C++ ClimateBlock and compute bioclimatic variables
#'
#' A convenience wrapper around the C++ \code{ClimateBlock} class exposed by
#' the \code{bioclim_mod} Rcpp module.  Accepts the same four monthly climate
#' vectors as \code{\link{bioclim}} but delegates the computation to the
#' compiled C++ back-end via Rcpp.
#'
#' @param tas    Numeric vector of length 12: monthly mean temperature.
#' @param tasmax Numeric vector of length 12: monthly maximum temperature.
#' @param tasmin Numeric vector of length 12: monthly minimum temperature.
#' @param pr     Numeric vector of length 12: monthly precipitation.
#'
#' @return A named numeric vector of length 19 (\code{bio01} through
#'   \code{bio19}), identical in meaning to the output of
#'   \code{\link{bioclim}}.
#'
#' @export
#' @examples
#' tas    <- c(5, 7, 10, 14, 18, 22, 25, 24, 20, 15, 10, 6)
#' tasmax <- c(8, 10, 14, 18, 23, 28, 32, 31, 26, 19, 13, 9)
#' tasmin <- c(1,  3,  6, 10, 13, 17, 20, 19, 15, 10,  6, 2)
#' pr     <- c(60, 55, 50, 40, 30, 15,  5, 10, 25, 45, 55, 65)
#' bioclim_block(tas, tasmax, tasmin, pr)
bioclim_block <- function(tas, tasmax, tasmin, pr) {
  # validate_monthly() is defined in R/primitives.R
  validate_monthly(tas,    "tas")
  validate_monthly(tasmax, "tasmax")
  validate_monthly(tasmin, "tasmin")
  validate_monthly(pr,     "pr")

  # Wrap each 12-element vector as a 1 × 12 matrix for ClimateBlock
  to_mat <- function(x) matrix(as.double(x), nrow = 1L, ncol = 12L)

  block  <- methods::new(ClimateBlock,
                         to_mat(tas), to_mat(tasmax),
                         to_mat(tasmin), to_mat(pr))
  result <- block$compute()

  # Return a named vector (same format as bioclim())
  result[1L, ]
}
