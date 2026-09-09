// gdal_mask_rcpp.cpp — Rcpp-exported polygon masking functions.
//
// Exports:
//   rasterize_mask_cpp(vector_path, ref_raster_path, output_mask_path)
//   apply_mask_cpp(input_path, mask_path, output_path)
//
// Both functions stop() with an informative message when the package was
// built without GDAL support.

#include <Rcpp.h>
#include "gdal_mask.hpp"

using namespace Rcpp;

//' Rasterize a vector polygon layer to a binary mask raster
//'
//' Burns all polygon features from a vector source into a new single-band
//' GeoTIFF raster that is spatially aligned to a reference raster.  Pixels
//' that fall inside at least one polygon are set to \code{1}; all other
//' pixels are set to \code{0}.
//'
//' This function is the low-level C++ entry point.  Most users should call
//' the higher-level \code{\link{create_mask}} wrapper instead.
//'
//' Stops with an informative error if the package was built without GDAL
//' support.
//'
//' @param vector_path Character string: path to any OGR-readable vector
//'   source (shapefile, GeoJSON, GeoPackage, etc.).
//' @param ref_raster_path Character string: path to a GDAL-readable raster
//'   used as the spatial reference (extent, resolution, CRS).
//' @param output_mask_path Character string: file path where the output
//'   \code{GDT_Byte} GeoTIFF will be written (created or overwritten).
//' @return Invisibly returns \code{NULL}.  The side effect is the creation
//'   of the output mask raster at \code{output_mask_path}.
//' @seealso \code{\link{create_mask}}, \code{\link{apply_mask_cpp}}
//' @examples
//' \donttest{
//' if (has_gdal()) {
//'   # Requires GDAL support at build time.
//'   ref  <- system.file("extdata", "tiny.tif", package = "xbioclim")
//'   poly <- tempfile(fileext = ".geojson")
//'   mask <- tempfile(fileext = ".tif")
//' writeLines(
//'   '{"type":"FeatureCollection","features":[{"type":"Feature",
//'     "geometry":{"type":"Polygon","coordinates":[[[0,0],[1,0],[1,1],[0,1],[0,0]]]},
//'     "properties":{}}]}',
//'   poly)
//'   rasterize_mask_cpp(poly, ref, mask)
//' }
//' }
//' @export
// [[Rcpp::export]]
SEXP rasterize_mask_cpp(const std::string& vector_path,   // NOLINT
                        const std::string& ref_raster_path,
                        const std::string& output_mask_path) {
    try {
        xbioclim::rasterize_mask(vector_path,
                                  ref_raster_path,
                                  output_mask_path);
    } catch (const std::exception& e) {
        Rcpp::stop(e.what());
    }
    return R_NilValue;
}

//' Apply a binary mask raster to an input raster
//'
//' Reads \code{input_path} and \code{mask_path} tile-by-tile (one row at a
//' time) and writes the result to \code{output_path}.  Pixels where the mask
//' equals \code{0} are replaced with \code{NaN} in the output; all other
//' pixels retain their original values (with any GDAL scale/offset applied).
//'
//' This function is the low-level C++ entry point.  Most users should call
//' the higher-level \code{\link{create_mask}} wrapper instead.
//'
//' Stops with an informative error if the package was built without GDAL
//' support or if the mask and input dimensions differ.
//'
//' @param input_path Character string: path to a GDAL-readable raster
//'   (any number of bands).
//' @param mask_path Character string: path to a single-band \code{GDT_Byte}
//'   mask raster (e.g., produced by \code{\link{rasterize_mask_cpp}}).
//' @param output_path Character string: file path where the output
//'   \code{Float64} GeoTIFF will be written (created or overwritten).
//' @return Invisibly returns \code{NULL}.  The side effect is the creation
//'   of the masked raster at \code{output_path}.
//' @seealso \code{\link{create_mask}}, \code{\link{rasterize_mask_cpp}}
//' @examples
//' \donttest{
//' if (has_gdal()) {
//'   # Requires GDAL support at build time.
//'   ref    <- system.file("extdata", "tiny.tif", package = "xbioclim")
//'   poly   <- tempfile(fileext = ".geojson")
//'   mask   <- tempfile(fileext = ".tif")
//'   output <- tempfile(fileext = ".tif")
//' writeLines(
//'   '{"type":"FeatureCollection","features":[{"type":"Feature",
//'     "geometry":{"type":"Polygon","coordinates":[[[0,0],[1,0],[1,1],[0,1],[0,0]]]},
//'     "properties":{}}]}',
//'   poly)
//'   rasterize_mask_cpp(poly, ref, mask)
//'   apply_mask_cpp(ref, mask, output)
//' }
//' }
//' @export
// [[Rcpp::export]]
SEXP apply_mask_cpp(const std::string& input_path,   // NOLINT
                    const std::string& mask_path,
                    const std::string& output_path) {
    try {
        xbioclim::apply_mask(input_path, mask_path, output_path);
    } catch (const std::exception& e) {
        Rcpp::stop(e.what());
    }
    return R_NilValue;
}
