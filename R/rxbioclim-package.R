#' rxbioclim: Bioclimatic Variables from Monthly Climate Data
#'
#' Computes the 19 standard bioclimatic variables (BIO01-BIO19) from monthly
#' climate data following the WorldClim specification. This is an R
#' implementation of the xbioclim C++ library.
#'
#' @section Error and warning propagation:
#' rxbioclim mirrors the `SpatMessages` pattern used by the terra package.
#' C++ routines record errors and warnings into an internal message store
#' rather than throwing directly.  After each C++ call the R-side wrapper
#' invokes [check_messages()] to convert any stored messages into native R
#' conditions.  Users can also inspect the store programmatically:
#' * [bioclim_errors()] / [bioclim_warnings()] – retrieve stored messages.
#' * [has_error()] / [has_warning()] – test whether messages exist.
#' * [clear_messages()] – reset the store.
#'
#' @docType package
#' @name rxbioclim-package
"_PACKAGE"
