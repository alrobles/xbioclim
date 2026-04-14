# R/engine.R — thin R wrapper for the BioclimEngine XPtr-based interface.
#
# The low-level XPtr functions (engine_create, engine_open, engine_set_output,
# engine_set_mask, engine_set_threads, engine_compute) are generated in
# RcppExports.R by Rcpp::compileAttributes().  This file adds:
#   * has_gdal()  — user-facing GDAL availability check

#' Check whether the package was built with GDAL support
#'
#' Returns \code{TRUE} when rxbioclim was compiled with GDAL and the native
#' \code{BioclimEngine} tiled pipeline is available, \code{FALSE} otherwise.
#'
#' Internally the function calls \code{\link{gdal_can_open}} with a dummy path
#' and inspects whether the resulting error message indicates an absent GDAL
#' build.
#'
#' @return Logical scalar: \code{TRUE} if GDAL is available, \code{FALSE}
#'   otherwise.
#' @export
#' @examples
#' has_gdal()
has_gdal <- function() {
  tryCatch({
    gdal_can_open("/dev/null")
    TRUE  # gdal_can_open() returned without throwing — GDAL is present
  }, error = function(e) {
    if (grepl("built without GDAL", conditionMessage(e), fixed = TRUE)) {
      FALSE  # stop() was called by the no-GDAL stub
    } else {
      TRUE   # GDAL is present but the path is invalid — that's fine
    }
  })
}
