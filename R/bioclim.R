#' Compute Bioclimatic Variables from Monthly Climate Data
#'
#' @description
#' Computes the 19 standard bioclimatic variables (BIO01-BIO19) from monthly
#' climate data following the WorldClim specification. This is an R
#' implementation of the xbioclim C++ library.
#'
#' @details
#' The 19 bioclimatic variables are:
#' \itemize{
#'   \item BIO01: Mean Annual Temperature
#'   \item BIO02: Mean Diurnal Range (mean of monthly (tasmax - tasmin))
#'   \item BIO03: Isothermality (100 * BIO02 / BIO07)
#'   \item BIO04: Temperature Seasonality (100 * population SD of monthly tas)
#'   \item BIO05: Max Temperature of Warmest Month
#'   \item BIO06: Min Temperature of Coldest Month
#'   \item BIO07: Temperature Annual Range (BIO05 - BIO06)
#'   \item BIO08: Mean Temperature of Wettest Quarter
#'   \item BIO09: Mean Temperature of Driest Quarter
#'   \item BIO10: Mean Temperature of Warmest Quarter
#'   \item BIO11: Mean Temperature of Coldest Quarter
#'   \item BIO12: Annual Precipitation
#'   \item BIO13: Precipitation of Wettest Month
#'   \item BIO14: Precipitation of Driest Month
#'   \item BIO15: Precipitation Seasonality (CV)
#'   \item BIO16: Precipitation of Wettest Quarter
#'   \item BIO17: Precipitation of Driest Quarter
#'   \item BIO18: Precipitation of Warmest Quarter
#'   \item BIO19: Precipitation of Coldest Quarter
#' }
#'
#' All functions accept either plain numeric vectors of length 12 (single
#' pixel) or a [BioclimData-class] object holding a raster block (multiple
#' pixels).  When a [BioclimData-class] object is supplied, computation is
#' delegated to the compiled C++ backend (xbioclim) and a numeric vector
#' (one value per pixel) is returned.
#'
#' @param tas    Numeric vector of length 12 **or** [BioclimData-class]:
#'   monthly mean temperature.
#' @param tasmax Numeric vector of length 12: monthly maximum temperature.
#' @param tasmin Numeric vector of length 12: monthly minimum temperature.
#' @param pr     Numeric vector of length 12: monthly precipitation.
#' @param ...    Additional arguments. For [BioclimData-class] inputs the
#'   argument `na.rm` is accepted: if `TRUE`, missing months are omitted and
#'   each BIO is computed from the available months (a quarter needs at least
#'   one valid month). The default `na.rm = FALSE` makes a pixel all-`NA` if
#'   any input month is `NA`.
#'
#' @return
#' * For `bio01`–`bio19` with plain numeric-vector inputs: a single numeric
#'   value.
#' * For `bioclim()` with plain numeric-vector inputs: a named numeric vector
#'   of length 19 (names `bio01`–`bio19`).
#' * For [BioclimData-class] inputs to `bio01`–`bio19`: a numeric vector with
#'   one value per pixel.
#' * For [BioclimData-class] inputs to `bioclim()`: a numeric matrix with one
#'   row per pixel and 19 named columns (`bio01`–`bio19`).
#'
#' @name bioclim-variables
#' @aliases bio01 bio02 bio03 bio04 bio05 bio06 bio07 bio08 bio09 bio10 bio11 bio12 bio13 bio14 bio15 bio16 bio17 bio18 bio19 bioclim
#' @aliases bio01,ANY-method bio01,BioclimData-method
#' @aliases bio02,ANY,ANY-method bio02,BioclimData,missing-method
#' @aliases bio03,ANY,ANY-method bio03,BioclimData,missing-method
#' @aliases bio04,ANY-method bio04,BioclimData-method
#' @aliases bio05,ANY-method bio05,BioclimData-method
#' @aliases bio06,ANY-method bio06,BioclimData-method
#' @aliases bio07,ANY,ANY-method bio07,BioclimData,missing-method
#' @aliases bio08,ANY,ANY-method bio08,BioclimData,missing-method
#' @aliases bio09,ANY,ANY-method bio09,BioclimData,missing-method
#' @aliases bio10,ANY-method bio10,BioclimData-method
#' @aliases bio11,ANY-method bio11,BioclimData-method
#' @aliases bio12,ANY-method bio12,BioclimData-method
#' @aliases bio13,ANY-method bio13,BioclimData-method
#' @aliases bio14,ANY-method bio14,BioclimData-method
#' @aliases bio15,ANY-method bio15,BioclimData-method
#' @aliases bio16,ANY-method bio16,BioclimData-method
#' @aliases bio17,ANY-method bio17,BioclimData-method
#' @aliases bio18,ANY,ANY-method bio18,BioclimData,missing-method
#' @aliases bio19,ANY,ANY-method bio19,BioclimData,missing-method
#' @aliases bioclim,ANY,ANY,ANY,ANY-method bioclim,BioclimData,missing,missing,missing-method
NULL

