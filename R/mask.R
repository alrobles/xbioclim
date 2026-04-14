#' Create a binary mask raster from a polygon source
#'
#' Converts a polygon source to a binary raster mask aligned to a reference
#' raster, and optionally applies that mask to an input raster (replacing
#' pixels outside all polygons with `NA`).
#'
#' The function accepts three types of polygon input:
#' \describe{
#'   \item{`sf` object}{Written to a temporary GeoJSON file via
#'     \code{\link[sf:st_write]{sf::st_write()}} and then rasterized.
#'     Requires the \pkg{sf} package.}
#'   \item{Character file path}{Passed directly to the C++ rasterizer.
#'     Any OGR-readable format is supported (shapefile, GeoJSON,
#'     GeoPackage, etc.).}
#'   \item{`SpatRaster` object}{Used directly as a pre-made mask raster.
#'     Must already be binary (0/1) and aligned to `reference_raster`.
#'     Requires the \pkg{terra} package.}
#' }
#'
#' All GDAL-dependent operations skip gracefully (returning `NULL` invisibly
#' with a message) when the package was built without GDAL support.
#'
#' @param polygon Polygon source: an `sf` object, a character file path to an
#'   OGR-readable vector source, or a `terra::SpatRaster` used directly as
#'   a binary mask.
#' @param reference_raster Character string: path to a GDAL-readable raster
#'   used to set the spatial reference (extent, resolution, CRS) of the
#'   output mask.  Ignored when `polygon` is a `SpatRaster`.
#' @param output Character string: file path for the output mask GeoTIFF.
#'   Defaults to a temporary file.
#' @param apply_to Character string or `NULL`: if supplied, path to a
#'   GDAL-readable raster to which the mask will be applied.  Pixels outside
#'   all polygons are set to `NA` in the result.
#' @param apply_output Character string: file path for the masked output
#'   raster.  Used only when `apply_to` is not `NULL`.  Defaults to a
#'   temporary file.
#'
#' @return When `apply_to` is `NULL`, returns the path to the binary mask
#'   GeoTIFF invisibly.  When `apply_to` is provided, returns the path to the
#'   masked output raster invisibly.
#'
#' @examples
#' \donttest{
#' # Requires GDAL support at build time.
#' ref  <- system.file("extdata", "tiny.tif", package = "rxbioclim")
#' poly <- tempfile(fileext = ".geojson")
#' writeLines(
#'   paste0(
#'     '{"type":"FeatureCollection","features":[{"type":"Feature",',
#'     '"geometry":{"type":"Polygon",',
#'     '"coordinates":[[[0,0],[1,0],[1,1],[0,1],[0,0]]]},',
#'     '"properties":{}}]}'
#'   ),
#'   poly
#' )
#' mask_path <- create_mask(poly, ref)
#' }
#' @seealso [rasterize_mask_cpp()], [apply_mask_cpp()]
#' @export
create_mask <- function(polygon,
                        reference_raster = NULL,
                        output           = tempfile(fileext = ".tif"),
                        apply_to         = NULL,
                        apply_output     = tempfile(fileext = ".tif")) {

  if (!is.character(output) || length(output) != 1L || is.na(output)) {
    stop("'output' must be a single non-NA character string.")
  }

  # ── 1. Resolve the mask raster path ─────────────────────────────────────
  mask_path <- .resolve_mask(polygon, reference_raster, output)
  if (is.null(mask_path)) return(invisible(NULL))

  # ── 2. Optionally apply the mask to another raster ──────────────────────
  if (!is.null(apply_to)) {
    if (!is.character(apply_to) || length(apply_to) != 1L || is.na(apply_to)) {
      stop("'apply_to' must be a single non-NA character string.")
    }
    if (!is.character(apply_output) || length(apply_output) != 1L || is.na(apply_output)) {
      stop("'apply_output' must be a single non-NA character string.")
    }
    tryCatch(
      apply_mask_cpp(apply_to, mask_path, apply_output),
      error = function(e) {
        if (grepl("built without GDAL", conditionMessage(e), fixed = TRUE)) {
          message(
            "create_mask: GDAL not available; apply_to step skipped.")
          return(invisible(NULL))
        }
        stop(e)
      }
    )
    return(invisible(apply_output))
  }

  invisible(mask_path)
}

# Internal helper — resolves `polygon` to a mask raster path.
# Returns the path on success, NULL with a message when GDAL is absent.
#' @keywords internal
.resolve_mask <- function(polygon, reference_raster, output) {

  # ── Case 1: SpatRaster ───────────────────────────────────────────────────
  if (inherits(polygon, "SpatRaster")) {
    if (!requireNamespace("terra", quietly = TRUE)) {
      stop(
        "create_mask: the 'terra' package is required to use a SpatRaster ",
        "as the polygon argument.  Install it with install.packages('terra').")
    }
    terra::writeRaster(polygon, output, overwrite = TRUE)
    return(output)
  }

  # ── Case 2: sf object ────────────────────────────────────────────────────
  if (inherits(polygon, "sf") || inherits(polygon, "sfc")) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop(
        "create_mask: the 'sf' package is required to use an sf object as ",
        "the polygon argument.  Install it with install.packages('sf').")
    }
    vec_tmp <- tempfile(fileext = ".geojson")
    sf::st_write(polygon, vec_tmp, quiet = TRUE)
    polygon <- vec_tmp
    on.exit(unlink(vec_tmp), add = TRUE)
  }

  # ── Case 3: file path ────────────────────────────────────────────────────
  if (!is.character(polygon) || length(polygon) != 1L || is.na(polygon)) {
    stop("'polygon' must be a single non-NA character string (file path).")
  }
  if (!is.character(reference_raster) || length(reference_raster) != 1L || is.na(reference_raster)) {
    stop("'reference_raster' must be a single non-NA character string.")
  }

  tryCatch(
    {
      rasterize_mask_cpp(polygon, reference_raster, output)
      output
    },
    error = function(e) {
      if (grepl("built without GDAL", conditionMessage(e), fixed = TRUE)) {
        message("create_mask: GDAL not available; mask creation skipped.")
        return(NULL)
      }
      stop(e)
    }
  )
}
