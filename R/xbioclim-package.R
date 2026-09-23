#' xbioclim: Bioclimatic Variables from Monthly Climate Data
#'
#' Computes the 19 standard bioclimatic variables (BIO01-BIO19) from monthly
#' climate data following the WorldClim specification. This is an R
#' implementation of the xbioclim C++ library, with a compiled C++ back-end
#' exposed through Rcpp Modules.
#'
#' @section ERA5-Land monthly aggregation:
#' In addition to the BIO01–BIO19 computation, xbioclim provides helpers to
#' aggregate ERA5-Land hourly `t2m`/`tp` into the CHELSA-compatible monthly
#' variables (`tas`, `tasmax`, `tasmin`, `pr`) used by the bioclim functions:
#' * [era5_to_monthly()] – single-pass hourly → monthly aggregation.
#' * [era5_t2m_to_monthly()] / [era5_tp_to_monthly()] – variable-specific helpers.
#' * [era5_bioclim()] / [era5_bioclim_years()] – end-to-end ERA5-Land → BIO pipeline.
#'
#' @section Error and warning propagation:
#' xbioclim mirrors the `SpatMessages` pattern used by the terra package.
#' C++ routines record errors and warnings into an internal message store
#' rather than throwing directly.  R-side wrappers around C++ calls should
#' invoke [check_messages()] after each call to convert any stored messages
#' into native R conditions.  Users can also inspect the store
#' programmatically:
#' * [bioclim_errors()] / [bioclim_warnings()] – retrieve stored messages.
#' * [has_error()] / [has_warning()] – test whether messages exist.
#' * [clear_messages()] – reset the store.
#'
#' @return No return value; package-level documentation.
#'
#' @docType package
#' @name xbioclim-package
#' @useDynLib xbioclim, .registration = TRUE
#' @import methods
#' @importFrom Rcpp evalCpp
"_PACKAGE"
