// BioclimEngine.cpp — implementation of the BioclimEngine tiled pipeline.
//
// Structure
// ---------
// 1. Member-function bodies that always compile (open, set_*, etc.).
// 2. An anonymous namespace (inside rxbioclim, guarded by #ifdef HAVE_GDAL)
//    with per-pixel bioclim helpers and GDAL I/O helpers.
// 3. BioclimEngine::compute() — thin always-compiled shell that delegates to
//    the GDAL path or throws a clear error.
// 4. Rcpp::export wrappers — always compiled, no GDAL dependency.

#ifdef _OPENMP
#include <omp.h>
#endif
#include <Rcpp.h>

#include "BioclimEngine.hpp"
#include "gdal_io.hpp"

#ifdef HAVE_CUDA
#include <cuda_runtime.h>
#include "bioclim_cuda.hpp"
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace rxbioclim {

// ── Always-compiled member implementations ───────────────────────────────────

void BioclimEngine::open(const std::vector<std::string>& tas_files,
                         const std::vector<std::string>& tasmax_files,
                         const std::vector<std::string>& tasmin_files,
                         const std::vector<std::string>& pr_files) {
    tas_files_    = tas_files;
    tasmax_files_ = tasmax_files;
    tasmin_files_ = tasmin_files;
    pr_files_     = pr_files;
}

void BioclimEngine::set_output(const std::string& path) {
    output_path_ = path;
}

void BioclimEngine::set_mask(const std::string& path) {
    mask_path_ = path;
}

void BioclimEngine::set_threads(int n) {
    n_threads_ = (n < 1) ? 1 : n;
}

void BioclimEngine::set_tile_size(int tile_size) {
    tile_size_ = (tile_size < 1) ? 1 : tile_size;
}

void BioclimEngine::set_device(const std::string& device) {
    if (device == "auto" || device == "Auto") {
        device_ = Device::Auto;
    } else if (device == "cpu" || device == "CPU") {
        device_ = Device::CPU;
    } else if (device == "gpu" || device == "GPU") {
        device_ = Device::GPU;
    } else {
        throw std::runtime_error(
            "BioclimEngine::set_device: unknown device '" + device +
            "'. Use \"auto\", \"cpu\", or \"gpu\".");
    }
}

// ── GDAL-dependent helpers (anonymous namespace, internal linkage) ────────────

#ifdef HAVE_GDAL

namespace {

// ── Per-pixel bioclim primitives (stack-only, OpenMP-safe) ──────────────────

// Rolling 3-month circular sums: qs[k] = x[k] + x[(k+1)%12] + x[(k+2)%12]
inline void rolling_quarter_sum(const double* x, double* qs) {
    for (int k = 0; k < 12; ++k)
        qs[k] = x[k] + x[(k + 1) % 12] + x[(k + 2) % 12];
}

// Index of maximum in a 12-element array (0-based)
inline int argmax12(const double* x) {
    int best = 0;
    for (int k = 1; k < 12; ++k) if (x[k] > x[best]) best = k;
    return best;
}

// Index of minimum in a 12-element array (0-based)
inline int argmin12(const double* x) {
    int best = 0;
    for (int k = 1; k < 12; ++k) if (x[k] < x[best]) best = k;
    return best;
}

// Population standard deviation (denominator N) for a 12-element array
inline double sd_pop(const double* x) {
    double s = 0.0, ss = 0.0;
    for (int i = 0; i < 12; ++i) { s += x[i]; ss += x[i] * x[i]; }
    const double m = s / 12.0;
    return std::sqrt(ss / 12.0 - m * m);
}

// Compute 19 bioclimatic variables for one pixel.
// Writes NaN to all 19 outputs when any input is NaN (covers R's NA_real_).
void compute_pixel(const double* t, const double* tmx,
                   const double* tmn, const double* p,
                   double* bio) {
    // NA / NaN guard: R's NA_real_ is a specific NaN, caught by std::isnan.
    for (int m = 0; m < 12; ++m) {
        if (std::isnan(t[m]) || std::isnan(tmx[m]) ||
            std::isnan(tmn[m]) || std::isnan(p[m])) {
            for (int j = 0; j < 19; ++j)
                bio[j] = std::numeric_limits<double>::quiet_NaN();
            return;
        }
    }

    // ── Temperature aggregates ──────────────────────────────────────────────
    double t_sum = 0.0, diurnal_sum = 0.0;
    double tmx_max = tmx[0], tmn_min = tmn[0];
    for (int m = 0; m < 12; ++m) {
        t_sum       += t[m];
        diurnal_sum += tmx[m] - tmn[m];
        if (tmx[m] > tmx_max) tmx_max = tmx[m];
        if (tmn[m] < tmn_min) tmn_min = tmn[m];
    }

    const double b01 = t_sum / 12.0;                                // BIO01
    const double b02 = diurnal_sum / 12.0;                          // BIO02
    const double b05 = tmx_max;                                      // BIO05
    const double b06 = tmn_min;                                      // BIO06
    const double b07 = b05 - b06;                                    // BIO07
    const double b03 = (b07 > 0.0) ? 100.0 * b02 / b07 : 0.0;     // BIO03
    const double b04 = 100.0 * sd_pop(t);                           // BIO04

    // ── Precipitation aggregates ────────────────────────────────────────────
    double p_sum = 0.0, p_max = p[0], p_min = p[0];
    for (int m = 0; m < 12; ++m) {
        p_sum += p[m];
        if (p[m] > p_max) p_max = p[m];
        if (p[m] < p_min) p_min = p[m];
    }
    const double b12 = p_sum;                                        // BIO12
    const double b13 = p_max;                                        // BIO13
    const double b14 = p_min;                                        // BIO14

    // BIO15: Precipitation Seasonality (CV = 100 * sd / (mean + 1))
    const double p_mean = p_sum / 12.0;
    double p_ssq = 0.0;
    for (int m = 0; m < 12; ++m) {
        const double d = p[m] - p_mean;
        p_ssq += d * d;
    }
    const double b15 = 100.0 * std::sqrt(p_ssq / 12.0) / (p_mean + 1.0); // BIO15

    // ── Rolling quarter sums ────────────────────────────────────────────────
    double pr_qs[12], t_qs[12];
    rolling_quarter_sum(p, pr_qs);
    rolling_quarter_sum(t, t_qs);

    const int wet_q  = argmax12(pr_qs);
    const int dry_q  = argmin12(pr_qs);
    const int warm_q = argmax12(t_qs);
    const int cold_q = argmin12(t_qs);

    bio[ 0] = b01;
    bio[ 1] = b02;
    bio[ 2] = b03;
    bio[ 3] = b04;
    bio[ 4] = b05;
    bio[ 5] = b06;
    bio[ 6] = b07;
    bio[ 7] = t_qs[wet_q]  / 3.0;   // BIO08
    bio[ 8] = t_qs[dry_q]  / 3.0;   // BIO09
    bio[ 9] = t_qs[warm_q] / 3.0;   // BIO10
    bio[10] = t_qs[cold_q] / 3.0;   // BIO11
    bio[11] = b12;
    bio[12] = b13;
    bio[13] = b14;
    bio[14] = b15;
    bio[15] = pr_qs[wet_q];          // BIO16
    bio[16] = pr_qs[dry_q];          // BIO17
    bio[17] = pr_qs[warm_q];         // BIO18
    bio[18] = pr_qs[cold_q];         // BIO19
}

// ── GDAL I/O helpers ─────────────────────────────────────────────────────────

// Validate that a file vector is acceptable (1 multi-band or 12 single-band).
void validate_file_vector(const std::vector<std::string>& files,
                          const char* varname) {
    if (files.empty()) {
        std::ostringstream oss;
        oss << "BioclimEngine: " << varname << " file list is empty.";
        throw std::runtime_error(oss.str());
    }
    if (files.size() != 1 && files.size() != 12) {
        std::ostringstream oss;
        oss << "BioclimEngine: " << varname
            << " must be 1 multi-band file or 12 single-band files, got "
            << files.size() << ".";
        throw std::runtime_error(oss.str());
    }
}

// Open readers for a variable's file list.
// Returns a vector of unique_ptr<GdalReader>:
//   size 1  when files.size() == 1 (multi-band)
//   size 12 when files.size() == 12 (one per month)
std::vector<std::unique_ptr<GdalReader>>
open_readers(const std::vector<std::string>& files) {
    std::vector<std::unique_ptr<GdalReader>> readers;
    if (files.size() == 1) {
        readers.push_back(std::make_unique<GdalReader>(files[0]));
    } else {
        readers.reserve(files.size());
        for (const auto& f : files)
            readers.push_back(std::make_unique<GdalReader>(f));
    }
    return readers;
}

// Read one month's tile window.
// multi_band == true  → readers[0], band = month + 1
// multi_band == false → readers[month], band = 1
void read_month_window(
        const std::vector<std::unique_ptr<GdalReader>>& readers,
        bool multi_band, int month,
        int xoff, int yoff, int xsize, int ysize,
        std::vector<double>& buf) {
    if (multi_band) {
        readers[0]->read_window(xoff, yoff, xsize, ysize, month + 1, buf);
    } else {
        readers[static_cast<std::size_t>(month)]->read_window(
            xoff, yoff, xsize, ysize, 1, buf);
    }
}

}  // anonymous namespace

#endif  // HAVE_GDAL

// ── BioclimEngine::compute() ─────────────────────────────────────────────────

std::string BioclimEngine::compute() {
#ifdef HAVE_GDAL
    // ── Validate configuration ──────────────────────────────────────────────
    validate_file_vector(tas_files_,    "tas");
    validate_file_vector(tasmax_files_, "tasmax");
    validate_file_vector(tasmin_files_, "tasmin");
    validate_file_vector(pr_files_,     "pr");

    if (output_path_.empty())
        throw std::runtime_error("BioclimEngine::compute: output path not set.");

    // ── Open readers ────────────────────────────────────────────────────────
    const bool tas_multi    = (tas_files_.size()    == 1);
    const bool tasmax_multi = (tasmax_files_.size() == 1);
    const bool tasmin_multi = (tasmin_files_.size() == 1);
    const bool pr_multi     = (pr_files_.size()     == 1);

    auto tas_readers    = open_readers(tas_files_);
    auto tasmax_readers = open_readers(tasmax_files_);
    auto tasmin_readers = open_readers(tasmin_files_);
    auto pr_readers     = open_readers(pr_files_);

    // Reference dimensions come from the first tas reader.
    const int nrows = tas_readers[0]->nrows();
    const int ncols = tas_readers[0]->ncols();
    const auto gt   = tas_readers[0]->geotransform();
    const auto crs  = tas_readers[0]->crs();

    // ── Open mask reader (optional) ─────────────────────────────────────────
    std::unique_ptr<GdalReader> mask_reader;
    if (!mask_path_.empty())
        mask_reader = std::make_unique<GdalReader>(mask_path_);

    // ── Create output dataset ───────────────────────────────────────────────
    GdalWriter writer(output_path_, nrows, ncols, 19, gt, crs);

    // ── Determine compute device ────────────────────────────────────────────
    bool use_gpu = false;
#ifdef HAVE_CUDA
    if (device_ != Device::CPU) {
        int gpu_count = 0;
        cudaError_t cuda_err = cudaGetDeviceCount(&gpu_count);
        if (cuda_err == cudaSuccess && gpu_count > 0) {
            use_gpu = true;
        }
    }
#endif

    // Scratch buffer reused across tiles and bands.
    std::vector<double> band_buf;

    // ── Tiled processing loop ───────────────────────────────────────────────
    const int ts = tile_size_;

    for (int yoff = 0; yoff < nrows; yoff += ts) {
        const int ysize = std::min(ts, nrows - yoff);
        for (int xoff = 0; xoff < ncols; xoff += ts) {
            const int xsize = std::min(ts, ncols - xoff);
            const int n_pix = xsize * ysize;

            // ── Read 4 × 12 monthly bands into tile buffers ─────────────────
            // Layout: tile[month * n_pix + pixel_index]
            std::vector<double> tas_tile(   static_cast<std::size_t>(12 * n_pix));
            std::vector<double> tasmax_tile(static_cast<std::size_t>(12 * n_pix));
            std::vector<double> tasmin_tile(static_cast<std::size_t>(12 * n_pix));
            std::vector<double> pr_tile(    static_cast<std::size_t>(12 * n_pix));

            for (int m = 0; m < 12; ++m) {
                read_month_window(tas_readers, tas_multi, m,
                                  xoff, yoff, xsize, ysize, band_buf);
                std::copy(band_buf.begin(), band_buf.end(),
                          tas_tile.begin() + m * n_pix);

                read_month_window(tasmax_readers, tasmax_multi, m,
                                  xoff, yoff, xsize, ysize, band_buf);
                std::copy(band_buf.begin(), band_buf.end(),
                          tasmax_tile.begin() + m * n_pix);

                read_month_window(tasmin_readers, tasmin_multi, m,
                                  xoff, yoff, xsize, ysize, band_buf);
                std::copy(band_buf.begin(), band_buf.end(),
                          tasmin_tile.begin() + m * n_pix);

                read_month_window(pr_readers, pr_multi, m,
                                  xoff, yoff, xsize, ysize, band_buf);
                std::copy(band_buf.begin(), band_buf.end(),
                          pr_tile.begin() + m * n_pix);
            }

            // ── Read mask tile (optional) ────────────────────────────────────
            std::vector<double> mask_buf;
            if (mask_reader)
                mask_reader->read_window(xoff, yoff, xsize, ysize, 1, mask_buf);

            // ── Compute 19 bio variables per pixel ───────────────────────────
            // Output layout: bio_tile[bio_index * n_pix + pixel_index]
            std::vector<double> bio_tile(static_cast<std::size_t>(19 * n_pix));

#ifdef HAVE_CUDA
            if (use_gpu) {
                launch_bioclim_cuda(
                    tas_tile.data(), tasmax_tile.data(),
                    tasmin_tile.data(), pr_tile.data(),
                    mask_buf.empty() ? nullptr : mask_buf.data(),
                    bio_tile.data(),
                    n_pix
                );
            } else
#endif
            {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads_)
#endif
            for (int i = 0; i < n_pix; ++i) {
                // Gather 12-month vectors from the tile buffers.
                double t[12], tmx[12], tmn[12], p[12], bio[19];
                for (int m = 0; m < 12; ++m) {
                    t[m]   = tas_tile[   static_cast<std::size_t>(m * n_pix + i)];
                    tmx[m] = tasmax_tile[static_cast<std::size_t>(m * n_pix + i)];
                    tmn[m] = tasmin_tile[static_cast<std::size_t>(m * n_pix + i)];
                    p[m]   = pr_tile[    static_cast<std::size_t>(m * n_pix + i)];
                }

                // Apply mask: 0 or NaN mask → all-NaN output.
                bool masked = false;
                if (!mask_buf.empty()) {
                    const double mv = mask_buf[static_cast<std::size_t>(i)];
                    masked = (std::isnan(mv) || mv == 0.0);
                }

                if (masked) {
                    for (int j = 0; j < 19; ++j)
                        bio_tile[static_cast<std::size_t>(j * n_pix + i)] =
                            std::numeric_limits<double>::quiet_NaN();
                } else {
                    compute_pixel(t, tmx, tmn, p, bio);
                    for (int j = 0; j < 19; ++j)
                        bio_tile[static_cast<std::size_t>(j * n_pix + i)] = bio[j];
                }
            }
            }  // end CPU branch

            // ── Write 19 output bands ────────────────────────────────────────
            std::vector<double> out_band(static_cast<std::size_t>(n_pix));
            for (int j = 0; j < 19; ++j) {
                for (int i = 0; i < n_pix; ++i)
                    out_band[static_cast<std::size_t>(i)] =
                        bio_tile[static_cast<std::size_t>(j * n_pix + i)];
                writer.write_window(xoff, yoff, xsize, ysize, j + 1, out_band);
            }
        }
    }

    writer.close();
    return output_path_;

#else
    throw std::runtime_error(
        "GDAL is required for BioclimEngine. "
        "Rebuild the package with GDAL support."
    );
    return "";  // unreachable — silences compiler warning
#endif
}

}  // namespace rxbioclim

// ── Rcpp XPtr wrappers ────────────────────────────────────────────────────────
//
// These functions are always compiled.  They create / manipulate
// BioclimEngine objects via opaque external pointers (Rcpp::XPtr).

//' Create a new BioclimEngine instance
//'
//' Allocates a new \code{BioclimEngine} C++ object and returns an opaque
//' external pointer to it.  Use the companion \code{engine_*()} functions to
//' configure and run the engine.
//'
//' @return An \code{externalptr} to a new \code{BioclimEngine} object.
//' @seealso \code{\link{engine_open}}, \code{\link{engine_compute}}
// [[Rcpp::export]]
SEXP engine_create() {
    Rcpp::XPtr<rxbioclim::BioclimEngine> ptr(
        new rxbioclim::BioclimEngine(), true);
    return ptr;
}

//' Configure monthly climate input files
//'
//' Associates four sets of raster file paths with the engine.  Each vector
//' must contain either one multi-band file (12 bands) or twelve single-band
//' files (one per calendar month).
//'
//' @param xptr   External pointer returned by \code{\link{engine_create}}.
//' @param tas_files    Character vector (length 1 or 12): mean temperature.
//' @param tasmax_files Character vector (length 1 or 12): maximum temperature.
//' @param tasmin_files Character vector (length 1 or 12): minimum temperature.
//' @param pr_files     Character vector (length 1 or 12): precipitation.
//' @return \code{NULL} invisibly.
//' @seealso \code{\link{engine_create}}, \code{\link{engine_compute}}
// [[Rcpp::export]]
void engine_open(SEXP xptr,
                 Rcpp::CharacterVector tas_files,
                 Rcpp::CharacterVector tasmax_files,
                 Rcpp::CharacterVector tasmin_files,
                 Rcpp::CharacterVector pr_files) {
    Rcpp::XPtr<rxbioclim::BioclimEngine> eng(xptr);
    eng->open(Rcpp::as<std::vector<std::string>>(tas_files),
              Rcpp::as<std::vector<std::string>>(tasmax_files),
              Rcpp::as<std::vector<std::string>>(tasmin_files),
              Rcpp::as<std::vector<std::string>>(pr_files));
}

//' Set the output raster path
//'
//' The engine will create (or overwrite) a Float64 GeoTIFF with 19 bands at
//' this path when \code{\link{engine_compute}} is called.
//'
//' @param xptr External pointer returned by \code{\link{engine_create}}.
//' @param path Character scalar: output file path.
//' @return \code{NULL} invisibly.
//' @seealso \code{\link{engine_create}}, \code{\link{engine_compute}}
// [[Rcpp::export]]
void engine_set_output(SEXP xptr, std::string path) {
    Rcpp::XPtr<rxbioclim::BioclimEngine>(xptr)->set_output(path);
}

//' Set an optional mask raster
//'
//' Pixels where the mask band equals 0 or \code{NaN} receive \code{NaN}
//' (no-data) in every output band.  Pass an empty string to disable masking.
//'
//' @param xptr      External pointer returned by \code{\link{engine_create}}.
//' @param mask_path Character scalar: mask raster path, or \code{""} for none.
//' @return \code{NULL} invisibly.
//' @seealso \code{\link{engine_create}}, \code{\link{engine_compute}}
// [[Rcpp::export]]
void engine_set_mask(SEXP xptr, std::string mask_path) {
    Rcpp::XPtr<rxbioclim::BioclimEngine>(xptr)->set_mask(mask_path);
}

//' Set the number of OpenMP threads
//'
//' Controls the number of threads used in the per-pixel inner loop inside
//' each tile.  Values less than 1 are clamped to 1.
//'
//' @param xptr External pointer returned by \code{\link{engine_create}}.
//' @param n    Integer scalar: number of threads.
//' @return \code{NULL} invisibly.
//' @seealso \code{\link{engine_create}}, \code{\link{engine_compute}}
// [[Rcpp::export]]
void engine_set_threads(SEXP xptr, int n) {
    Rcpp::XPtr<rxbioclim::BioclimEngine>(xptr)->set_threads(n);
}

//' Set the tile size used during tiled processing
//'
//' Width and height of each processing tile in pixels.  Default is 256.
//' Mostly useful for testing with small rasters.  Values less than 1 are
//' clamped to 1.
//'
//' @param xptr      External pointer returned by \code{\link{engine_create}}.
//' @param tile_size Integer scalar: tile width and height in pixels.
//' @return \code{NULL} invisibly.
//' @seealso \code{\link{engine_create}}, \code{\link{engine_compute}}
//' @keywords internal
// [[Rcpp::export]]
void engine_set_tile_size(SEXP xptr, int tile_size) {
    Rcpp::XPtr<rxbioclim::BioclimEngine>(xptr)->set_tile_size(tile_size);
}

//' Run the bioclimatic-variable computation pipeline
//'
//' Reads all monthly climate input rasters tile by tile, computes the 19
//' bioclimatic variables (BIO01–BIO19) for every pixel, and writes the results
//' to the output raster.  Peak memory is proportional to the tile size, not
//' the full raster size.
//'
//' Requires GDAL support.  Stops with an informative error when the package
//' was built without GDAL.
//'
//' @param xptr External pointer returned by \code{\link{engine_create}}.
//' @return Character scalar: the output file path (same as the value passed to
//'   \code{\link{engine_set_output}}).
//' @seealso \code{\link{engine_create}}, \code{\link{engine_set_output}},
//'   \code{\link{has_gdal}}
// [[Rcpp::export]]
std::string engine_compute(SEXP xptr) {
    Rcpp::XPtr<rxbioclim::BioclimEngine> eng(xptr);
    return eng->compute();
}
