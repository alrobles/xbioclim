#' xbioclim: Bioclimatic Variables from Monthly Climate Data
#'
#' Computes the 19 standard bioclimatic variables (BIO01-BIO19) from monthly
#' climate data following the WorldClim specification. This is an R
#' implementation of the xbioclim C++ library, with a compiled C++ back-end
#' exposed through Rcpp Modules.
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
#' @docType package
#' @name xbioclim-package
#' @useDynLib xbioclim, .registration = TRUE
#' @import methods
#' @importFrom Rcpp evalCpp
"_PACKAGE"
