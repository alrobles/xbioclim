// gdal_io.hpp — GDAL-backed tiled raster I/O for xbioclim.
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
//   by default.  Pass dtype = GDT_Float32 for Float32 output, or
//   cog_compatible = true to add GTiff tiling + LZW compression options
//   compatible with Cloud-Optimised GeoTIFF.
//
// Thread safety
// -------------
// Each GdalReader / GdalWriter object owns its GDALDataset exclusively.
// Do not share a single object across threads; open one per thread instead.

#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef HAVE_GDAL
// Minimal stand-in for GDALDataType when GDAL is not available, so the
// public API and BioclimEngine members still compile.
enum GDALDataType {
    GDT_Unknown = -1,
    GDT_Float64 = 0,
    GDT_Float32 = 1
};
#else
#include <gdal_priv.h>
#include <cpl_conv.h>   // CPLMalloc / CPLFree
#include <ogr_spatialref.h>
#endif

namespace xbioclim {

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

    // Same as above, but writes directly into a caller-provided buffer.
    //
    // pixel_stride is the distance between consecutive pixels in `out`,
    // measured in doubles (default 1 = packed).  A stride greater than 1
    // interleaves this band with other data: the value of pixel i is
    // written to out[i * pixel_stride], so the buffer must hold at least
    // (xsize * ysize - 1) * pixel_stride + 1 doubles.  This is used to
    // write directly into pixel-major (band-interleaved) tile buffers.
    void read_window(int xoff, int yoff, int xsize, int ysize,
                     int band, double* out,
                     std::size_t pixel_stride = 1) const;

    // Read a rectangular window from multiple bands in one GDAL call.
    // Bands are 1-based.  Scale and offset are applied per band.
    //
    // When pixel_major is false (default) the output is stored band-major:
    // buf[band * n_pix + pixel].  When pixel_major is true the output is
    // pixel-major (band-interleaved): buf[pixel * bands.size() + band].
    void read_bands_window(int xoff, int yoff, int xsize, int ysize,
                           const std::vector<int>& bands,
                           std::vector<double>& buf,
                           bool pixel_major = false) const;

    // Same, but writes into a caller-provided buffer of at least
    // bands.size() * xsize * ysize doubles.
    void read_bands_window(int xoff, int yoff, int xsize, int ysize,
                           const std::vector<int>& bands,
                           double* out,
                           bool pixel_major = false) const;

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
    //   dtype         — GDAL band data type (GDT_Float64 or GDT_Float32)
    //
    // Throws std::runtime_error if the dataset cannot be created or if GDAL
    // is not compiled in.
    GdalWriter(const std::string& path,
               int nrows, int ncols, int nbands,
               const std::vector<double>& geotransform,
               const std::string& crs,
               bool cog_compatible = false,
               GDALDataType dtype = GDT_Float64);

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

    // Same, but with explicit buffer/output GDAL data type.
    void write_window(int xoff, int yoff, int xsize, int ysize,
                      int band, const std::vector<double>& buf,
                      GDALDataType dtype);

    // Write a rectangular window to multiple bands in one GDAL call.
    // buf must contain bands.size() * xsize * ysize elements.  bands are
    // 1-based output band indices.
    //
    // When pixel_major is false (default) buf is band-major:
    // buf[band * n_pix + pixel].  When pixel_major is true buf is
    // pixel-major (band-interleaved): buf[pixel * bands.size() + band].
    void write_bands_window(int xoff, int yoff, int xsize, int ysize,
                            const std::vector<int>& bands,
                            const std::vector<double>& buf,
                            bool pixel_major = false);

    // Same, but with explicit buffer/output GDAL data type.
    void write_bands_window(int xoff, int yoff, int xsize, int ysize,
                            const std::vector<int>& bands,
                            const std::vector<double>& buf,
                            GDALDataType dtype,
                            bool pixel_major = false);

    // Flush caches and close the dataset.  Safe to call more than once.
    void close();

private:
#ifdef HAVE_GDAL
    GDALDataset* ds_ = nullptr;
#endif
    GDALDataType dtype_ = GDT_Float64;
    bool closed_ = false;
};

} // namespace xbioclim
