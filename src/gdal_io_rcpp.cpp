// gdal_io_rcpp.cpp — Rcpp-exported diagnostic functions for GDAL integration.
//
// Exports two functions:
//   gdal_can_open(path)  — logical: can GDAL open this file?
//   gdal_info(path)      — named list of raster metadata
//
// Both functions stop() with an informative message when GDAL is not compiled
// in (i.e. HAVE_GDAL is not defined).

#include <Rcpp.h>
#include "gdal_io.hpp"

using namespace Rcpp;

//' Check whether GDAL can open a raster file
//'
//' A lightweight diagnostic that tries to open the specified path via GDAL
//' and returns \code{TRUE} if successful, \code{FALSE} if GDAL cannot open
//' it.  Stops with an informative error if the package was built without
//' GDAL support.
//'
//' @param path Character string: path to the raster file.
//' @return Logical \code{TRUE} if GDAL can open the file, \code{FALSE}
//'   otherwise.
//' @examples
//' \donttest{
//' if (has_gdal()) {
//'   gdal_can_open(system.file("extdata", "tiny.tif", package = "xbioclim"))
//' }
//' }
//' @export
// [[Rcpp::export]]
bool gdal_can_open(const std::string& path) {  // NOLINT
#ifdef HAVE_GDAL
    try {
        xbioclim::GdalReader reader(path);
        (void)reader;  // success — suppress unused-variable warning
        return true;
    } catch (const std::exception&) {
        return false;
    }
#else
    (void)path;
    Rcpp::stop(
        "gdal_can_open(): xbioclim was built without GDAL support. "
        "Install GDAL >= 2.0.1 and reinstall the package."
    );
    return false;  // unreachable, but silences compiler warnings
#endif
}

//' Return metadata about a GDAL-readable raster
//'
//' Opens the raster at \code{path} and returns a named list with dimensions,
//' geotransform, coordinate reference system, and per-band scale/offset
//' values.
//'
//' Data values are always returned as \code{double} (\code{float64}) in
//' memory.  If the raster bands carry GDAL scale/offset metadata (e.g. packed
//' integers), those are reported here and applied automatically by
//' \code{GdalReader::read_window()}.
//'
//' Stops with an informative error if the package was built without GDAL
//' support.
//'
//' @param path Character string: path to the raster file.
//' @return Named list with the following elements:
//'   \describe{
//'     \item{\code{path}}{Character: the path as supplied.}
//'     \item{\code{nrows}}{Integer: number of rows (Y pixels).}
//'     \item{\code{ncols}}{Integer: number of columns (X pixels).}
//'     \item{\code{nbands}}{Integer: number of raster bands.}
//'     \item{\code{geotransform}}{Numeric vector of length 6 (GDAL
//'       convention): \code{[x_origin, pixel_width, rotation_x,
//'       y_origin, rotation_y, pixel_height]}.}
//'     \item{\code{crs}}{Character: WKT coordinate reference system string,
//'       or an empty string if not defined.}
//'     \item{\code{scale}}{Numeric vector (one per band): GDAL scale factor
//'       (\code{1.0} if not set).}
//'     \item{\code{offset}}{Numeric vector (one per band): GDAL offset
//'       (\code{0.0} if not set).}
//'   }
//' @examples
//' \donttest{
//' if (has_gdal()) {
//'   info <- gdal_info(
//'     system.file("extdata", "tiny.tif", package = "xbioclim")
//'   )
//'   str(info)
//' }
//' }
//' @export
// [[Rcpp::export]]
Rcpp::List gdal_info(const std::string& path) {
#ifdef HAVE_GDAL
    xbioclim::GdalReader reader(path);

    const int nb = reader.nbands();

    // Per-band scale and offset vectors.
    Rcpp::NumericVector scale_vec(nb), offset_vec(nb);
    for (int b = 1; b <= nb; ++b) {
        scale_vec[b - 1]  = reader.scale(b);
        offset_vec[b - 1] = reader.offset(b);
    }

    return Rcpp::List::create(
        Rcpp::Named("path")         = path,
        Rcpp::Named("nrows")        = reader.nrows(),
        Rcpp::Named("ncols")        = reader.ncols(),
        Rcpp::Named("nbands")       = nb,
        Rcpp::Named("geotransform") = Rcpp::wrap(reader.geotransform()),
        Rcpp::Named("crs")          = reader.crs(),
        Rcpp::Named("scale")        = scale_vec,
        Rcpp::Named("offset")       = offset_vec
    );
#else
    (void)path;
    Rcpp::stop(
        "gdal_info() requires GDAL but xbioclim was built without it. "
        "Install GDAL >= 2.0.1 and reinstall the package."
    );
    return Rcpp::List();  // unreachable
#endif
}
