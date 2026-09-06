# R/engine.R — thin R wrapper for the BioclimEngine XPtr-based interface.
#
# The low-level XPtr functions (engine_create, engine_open, engine_set_output,
# engine_set_mask, engine_set_threads, engine_compute) are generated in
# RcppExports.R by Rcpp::compileAttributes().  This file adds:
#   * has_gdal()  — user-facing GDAL availability check
#   * has_cuda()  — user-facing CUDA GPU availability check
#   * cuda_info() — CUDA device information

#' Check whether the package was built with GDAL support
#'
#' Returns \code{TRUE} when xbioclim was compiled with GDAL and the native
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
    gdal_can_open(tempfile(fileext = ".tif"))
    TRUE  # gdal_can_open() returned without throwing — GDAL is present
  }, error = function(e) {
    if (grepl("built without GDAL", conditionMessage(e), fixed = TRUE)) {
      FALSE  # stop() was called by the no-GDAL stub
    } else {
      TRUE   # GDAL is present but the path is invalid — that's fine
    }
  })
}

#' Check CUDA GPU availability
#'
#' Returns \code{TRUE} when at least one CUDA-capable GPU is detected at
#' runtime.  Returns \code{FALSE} when the package was built without CUDA
#' support or when no CUDA-capable device is present.
#'
#' @return Logical scalar: \code{TRUE} if at least one CUDA GPU is available.
#' @seealso \code{\link{cuda_info}}, \code{\link{bioclim_engine}}
#' @export
#' @examples
#' has_cuda()
has_cuda <- function() {
  tryCatch(cuda_device_count() > 0L, error = function(e) FALSE)
}

#' Query CUDA GPU device information
#'
#' Returns hardware information about the first CUDA-capable GPU detected on
#' this machine.
#'
#' @return A named list with \code{name} (character GPU model name),
#'   \code{memory_gb} (numeric total memory in GB), and
#'   \code{compute_capability} (character, e.g. \code{"8.0"} for A100).
#'   Returns an empty list when no CUDA device is available.
#' @seealso \code{\link{has_cuda}}, \code{\link{bioclim_engine}}
#' @export
#' @examples
#' cuda_info()
cuda_info <- function() {
  tryCatch(cuda_device_info(), error = function(e) list())
}
