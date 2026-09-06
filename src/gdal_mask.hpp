// gdal_mask.hpp — Polygon-to-raster masking for xbioclim.
//
// Provides two operations:
//
//   rasterize_mask  — burn a vector polygon layer into a binary byte raster
//                     (1 = inside polygon, 0 = outside) that is spatially
//                     aligned to a reference raster.
//
//   apply_mask      — copy an input raster to an output raster, replacing
//                     every pixel where the mask == 0 with NaN.
//
// All GDAL-dependent code is inside #ifdef HAVE_GDAL blocks.  When GDAL is
// absent both functions throw std::runtime_error with a clear message so the
// package still compiles and links cleanly on platforms without GDAL.
//
// Thread safety
// -------------
// These functions are standalone (stateless); they are safe to call from
// concurrent threads provided each call uses distinct file paths.

#pragma once

#include <string>

#ifdef HAVE_GDAL
#include <gdal_priv.h>
#include <ogr_api.h>
#endif

namespace xbioclim {

// Rasterize a vector polygon source into a temporary binary mask GeoTIFF.
//
// Parameters:
//   vector_path   — path to any OGR-readable vector source (shapefile,
//                   GeoJSON, GPKG, …); all layers are burned.
//   ref_raster    — path to a GDAL-readable raster used as the spatial
//                   reference (extent, resolution, CRS).
//   output_path   — path for the output GeoTIFF (created / overwritten).
//                   Pixel type is GDT_Byte; value 1 = inside, 0 = outside.
//
// Throws std::runtime_error on any error or when GDAL is absent.
void rasterize_mask(const std::string& vector_path,
                    const std::string& ref_raster,
                    const std::string& output_path);

// Apply a binary mask raster to an input raster, writing results to an
// output raster.
//
// Pixels where mask == 0 are replaced with NaN in the output.  Pixels where
// mask != 0 keep the original input value.  Scale / offset metadata from the
// input raster is respected (values are read as physical float64).
//
// Parameters:
//   input_path    — path to a GDAL-readable raster (any number of bands).
//   mask_path     — path to a single-band GDT_Byte mask raster produced by
//                   rasterize_mask() (or any compatible binary mask).
//   output_path   — path for the output GeoTIFF (Float64, same number of
//                   bands as the input).
//
// Throws std::runtime_error on any error or when GDAL is absent.
void apply_mask(const std::string& input_path,
                const std::string& mask_path,
                const std::string& output_path);

} // namespace xbioclim
