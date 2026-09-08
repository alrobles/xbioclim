// BioclimEngine.hpp — tiled bioclimatic-variable engine for xbioclim.
//
// BioclimEngine orchestrates the full pipeline:
//   open() → set_output() → [set_mask()] → [set_threads()] → compute()
//
// All data I/O and computation are handled in C++; R sends only file paths
// and scalar parameters.  Peak memory is proportional to tile size, not
// raster size.
//
// GDAL dependency
// ---------------
// compute() requires GDAL.  When the package is built without GDAL it throws
// a clear std::runtime_error (re-raised as Rcpp::stop() by the Rcpp wrapper).
// All other methods and the constructor compile and run on every platform.

#pragma once

#include "gdal_io.hpp"

#include <string>
#include <vector>

namespace xbioclim {

class BioclimEngine {
public:
    BioclimEngine()  = default;
    ~BioclimEngine() = default;

    // Non-copyable; moving is disabled to keep the API simple.
    BioclimEngine(const BioclimEngine&)            = delete;
    BioclimEngine& operator=(const BioclimEngine&) = delete;
    BioclimEngine(BioclimEngine&&)                 = delete;
    BioclimEngine& operator=(BioclimEngine&&)      = delete;

    // ── Input ──────────────────────────────────────────────────────────────

    // Set monthly climate input files for the four variables.
    //
    // Each vector must contain either:
    //   * 1 element  — a single multi-band file with bands 1..12, or
    //   * 12 elements — one single-band file per calendar month.
    //
    // Throws std::runtime_error from compute() (not here) if the sizes are
    // wrong or the files cannot be opened.
    void open(const std::vector<std::string>& tas_files,
              const std::vector<std::string>& tasmax_files,
              const std::vector<std::string>& tasmin_files,
              const std::vector<std::string>& pr_files);

    // ── Output ─────────────────────────────────────────────────────────────

    // Set the output directory.  compute() creates (or overwrites) a single
    // 19-band GeoTIFF named bio.tif inside this directory.
    void set_output(const std::string& path);

    // ── Optional settings ──────────────────────────────────────────────────

    // Set an optional mask raster.  Pixels where the mask value is 0 or NaN
    // receive all-NaN output.  Pass an empty string to disable masking.
    void set_mask(const std::string& path);

    // Set the number of OpenMP threads used in the per-pixel inner loop.
    // Values < 1 are clamped to 1.
    void set_threads(int n);

    // Override the tile size (width and height in pixels).  Default: 256.
    // Mostly useful for tests.  Values < 1 are clamped to 1.
    void set_tile_size(int tile_size);

    // Set the compute device.  Accepted values: "auto", "cpu", "gpu".
    // "auto" (default) uses the GPU when at least one CUDA device is present,
    // otherwise falls back to the CPU.  "gpu" without CUDA compiled in or
    // without a visible device silently falls back to CPU.
    void set_device(const std::string& device);

    // Set the output GDAL data type.  Accepted values: "Float64", "Float32".
    // Default is "Float64".
    void set_dtype(const std::string& dtype);

    // Select which bioclimatic variables to return.  Each element must be in
    // [1, 19].  An empty vector (the default) means "all 19".
    // The engine always computes all 19 internally (they share intermediate
    // values), and writes a single 19-band output file; the R wrapper
    // subsets the returned SpatRaster when a selection is requested.
    void set_variables(const std::vector<int>& variables);

    // ── Execution ──────────────────────────────────────────────────────────

    // Run the tiled pipeline and return the output directory path.
    //
    // Throws std::runtime_error when:
    //   * GDAL is not compiled in
    //   * input or output paths are not configured
    //   * any file cannot be opened or read
    std::string compute();

private:
    std::vector<std::string> tas_files_;
    std::vector<std::string> tasmax_files_;
    std::vector<std::string> tasmin_files_;
    std::vector<std::string> pr_files_;

    std::string output_path_;
    std::string mask_path_;

    int n_threads_ = 1;
    int tile_size_ = 256;

    // Compute device: 0 = auto, 1 = CPU, 2 = GPU
    enum class Device { Auto, CPU, GPU };
    Device device_ = Device::Auto;

    // GDAL band data type for the output file (Float64 by default).
    GDALDataType output_dtype_ = GDT_Float64;

    // 1-based indices of variables selected by the R wrapper; empty = all 19.
    // compute() always writes all 19 bands, so this only affects the
    // post-processing done by bioclim_engine().
    std::vector<int> variables_;
};

}  // namespace xbioclim
