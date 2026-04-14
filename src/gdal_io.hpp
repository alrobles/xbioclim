// gdal_io.hpp — GDAL-backed tiled raster I/O for rxbioclim.
//
// Provides GdalReader and GdalWriter C++ classes that implement windowed
// (tiled) reads and writes so that only 1–2 tiles need to reside in memory
// at any time.  All public API is compiled unconditionally; the implementation
// bodies are guarded by #ifdef HAVE_GDAL so the package still compiles (and
// links) when GDAL is absent.
//
// Data model
// ----------
// * GdalReader always returns double-precision values.  If the on-disk type
//   is a packed integer and the band carries GDAL scale/offset metadata, those
//   are applied before returning the data.
// * GdalWriter accepts double-precision input and writes Float64 GTiff bands
//   by default.  Pass cog_compatible = true to add GTiff tiling + LZW
//   compression options compatible with Cloud-Optimised GeoTIFF.
//
// Thread safety
// -------------
// Each GdalReader / GdalWriter object owns its GDALDataset exclusively.
// Do not share a single object across threads; open one per thread instead.

#pragma once

#include <stdexcept>
#include <string>
#include <vector>

#ifdef HAVE_GDAL
#include <gdal_priv.h>
#include <cpl_conv.h>   // CPLMalloc / CPLFree
#include <ogr_spatialref.h>
#endif

namespace rxbioclim {

// ── GdalReader ──────────────────────────────────────────────────────────────

class GdalReader {
public:
    // Open a raster for reading.
    // Throws std::runtime_error if the file cannot be opened or if GDAL is
    // not compiled in.
    explicit GdalReader(const std::string& path);

    // The destructor closes the underlying dataset.
    ~GdalReader();

    // Non-copyable, movable.
    GdalReader(const GdalReader&)            = delete;
    GdalReader& operator=(const GdalReader&) = delete;
    GdalReader(GdalReader&&)                 = default;
    GdalReader& operator=(GdalReader&&)      = default;

    // Raster dimensions (pixels).
    int nrows() const;
    int ncols() const;
    int nbands() const;

    // 6-element GDAL geotransform vector.
    std::vector<double> geotransform() const;

    // WKT coordinate reference system string (empty if not defined).
    std::string crs() const;

    // Per-band scale / offset used when reading packed integer rasters.
    // Band index is 1-based.  Returns 1.0 / 0.0 if not set.
    double scale(int band) const;
    double offset(int band) const;

    // Read a rectangular window from one band into buf (resized to
    // xsize * ysize).  Scale and offset are applied when the band carries
    // GDAL metadata for them, so the returned values are always physical
    // (float64) quantities.
    //
    // Parameters (0-based pixel coordinates):
    //   xoff, yoff  — upper-left corner of the window
    //   xsize, ysize — width and height of the window in pixels
    //   band        — 1-based band index
    void read_window(int xoff, int yoff, int xsize, int ysize,
                     int band, std::vector<double>& buf) const;

private:
#ifdef HAVE_GDAL
    GDALDataset* ds_ = nullptr;
#endif
    std::string  path_;
};

// ── GdalWriter ──────────────────────────────────────────────────────────────

class GdalWriter {
public:
    // Create a new GTiff for writing.
    //
    // Parameters:
    //   path          — output file path (will be overwritten if it exists)
    //   nrows, ncols  — raster dimensions in pixels
    //   nbands        — number of bands
    //   geotransform  — 6-element GDAL geotransform (may be empty → no-op)
    //   crs           — WKT CRS string (may be empty → no-op)
    //   cog_compatible — when true, adds TILED=YES COMPRESS=LZW options
    //                    suitable for a Cloud-Optimised GeoTIFF workflow
    //
    // Throws std::runtime_error if the dataset cannot be created or if GDAL
    // is not compiled in.
    GdalWriter(const std::string& path,
               int nrows, int ncols, int nbands,
               const std::vector<double>& geotransform,
               const std::string& crs,
               bool cog_compatible = false);

    // The destructor flushes and closes the dataset if close() was not already
    // called.
    ~GdalWriter();

    // Non-copyable, movable.
    GdalWriter(const GdalWriter&)            = delete;
    GdalWriter& operator=(const GdalWriter&) = delete;
    GdalWriter(GdalWriter&&)                 = default;
    GdalWriter& operator=(GdalWriter&&)      = default;

    // Write a rectangular window to one band from buf.  buf must contain
    // exactly xsize * ysize elements (row-major).
    //
    // Parameters (0-based pixel coordinates):
    //   xoff, yoff   — upper-left corner of the window
    //   xsize, ysize — width and height of the window in pixels
    //   band         — 1-based band index
    void write_window(int xoff, int yoff, int xsize, int ysize,
                      int band, const std::vector<double>& buf);

    // Flush caches and close the dataset.  Safe to call more than once.
    void close();

private:
#ifdef HAVE_GDAL
    GDALDataset* ds_ = nullptr;
#endif
    bool closed_ = false;
};

} // namespace rxbioclim
