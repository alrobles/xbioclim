// BioclimEngine.cpp — implementation of the BioclimEngine tiled pipeline.
//
// Structure
// ---------
// 1. Member-function bodies that always compile (open, set_*, etc.).
// 2. An anonymous namespace (inside xbioclim, guarded by #ifdef HAVE_GDAL)
//    with per-pixel bioclim helpers and GDAL I/O helpers.
// 3. BioclimEngine::compute() — thin always-compiled shell that delegates to
//    the GDAL path or throws a clear error.
// 4. Rcpp::export wrappers — always compiled, no GDAL dependency.

#ifdef _OPENMP
#include <omp.h>
#endif
#include "xbioclim_omp.h"
#include <Rcpp.h>

#include "BioclimEngine.hpp"
#include "gdal_io.hpp"
#include "xbioclim_omp.h"

#ifdef HAVE_CUDA
#include <cuda_runtime.h>
#include "bioclim_cuda.hpp"
#endif

#include <algorithm>
#include <atomic>
#include <cmath>
#include <condition_variable>
#include <deque>
#include <exception>
#include <iomanip>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <vector>

namespace xbioclim {

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
    n_threads_ = xbioclim_safe_threads(n);
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

void BioclimEngine::set_dtype(const std::string& dtype) {
    if (dtype == "Float64" || dtype == "float64") {
        output_dtype_ = GDT_Float64;
    } else if (dtype == "Float32" || dtype == "float32") {
        output_dtype_ = GDT_Float32;
    } else {
        throw std::runtime_error(
            "BioclimEngine::set_dtype: unsupported dtype '" + dtype +
            "'. Use \"Float64\" or \"Float32\".");
    }
}

void BioclimEngine::set_variables(const std::vector<int>& variables) {
    for (int v : variables) {
        if (v < 1 || v > 19) {
            std::ostringstream oss;
            oss << "BioclimEngine::set_variables: variable index must be "
                   "1-19, got " << v << ".";
            throw std::runtime_error(oss.str());
        }
    }
    variables_ = variables;
}

void BioclimEngine::set_pipeline(bool use_pipeline) {
    pipeline_ = use_pipeline;
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
    double s = 0.0;
    for (int i = 0; i < 12; ++i) s += x[i];
    const double m = s / 12.0;
    double ss = 0.0;
    for (int i = 0; i < 12; ++i) {
        const double d = x[i] - m;
        ss += d * d;
    }
    return std::sqrt(ss / 12.0);
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

    // BIO15: Precipitation Seasonality (CV = 100 * sd / mean; NaN when mean == 0)
    const double p_mean = p_sum / 12.0;
    double p_ssq = 0.0;
    for (int m = 0; m < 12; ++m) {
        const double d = p[m] - p_mean;
        p_ssq += d * d;
    }
    const double b15 = (p_mean == 0.0)
        ? R_NaN
        : 100.0 * std::sqrt(p_ssq / 12.0) / p_mean;                 // BIO15

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

// Compute a single pixel from the month-major tile buffers and write the
// result into the band-major output tile.  Used as the scalar tail for tiles
// whose pixel count is not a multiple of the SIMD width.
inline void compute_pixel_at(const double* tas, const double* tasmax,
                             const double* tasmin, const double* pr,
                             const double* mask, double* bio,
                             int n_pix, int i) {
    double t[12], tmx[12], tmn[12], p[12], out[19];
    for (int m = 0; m < 12; ++m) {
        const std::size_t off = static_cast<std::size_t>(m) * n_pix + i;
        t[m]   = tas[off];
        tmx[m] = tasmax[off];
        tmn[m] = tasmin[off];
        p[m]   = pr[off];
    }

    bool masked = false;
    if (mask != nullptr) {
        const double mv = mask[static_cast<std::size_t>(i)];
        masked = (std::isnan(mv) || mv == 0.0);
    }

    if (masked) {
        for (int j = 0; j < 19; ++j)
            bio[static_cast<std::size_t>(j) * n_pix + i] =
                std::numeric_limits<double>::quiet_NaN();
    } else {
        compute_pixel(t, tmx, tmn, p, out);
        for (int j = 0; j < 19; ++j)
            bio[static_cast<std::size_t>(j) * n_pix + i] = out[j];
    }
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
// Always returns 12 unique_ptr<GdalReader>, one per calendar month:
//   * files.size() == 1  — 12 readers on the same multi-band file
//                          (month m reads band m+1)
//   * files.size() == 12 — one reader per single-band file
//                          (every month reads band 1)
// Per-month handles let the monthly window reads run in parallel —
// GDAL datasets are not thread-safe, so each concurrent read needs its
// own GdalReader.
std::vector<std::unique_ptr<GdalReader>>
open_readers(const std::vector<std::string>& files) {
    std::vector<std::unique_ptr<GdalReader>> readers;
    readers.reserve(12);
    if (files.size() == 1) {
        for (int m = 0; m < 12; ++m)
            readers.push_back(std::make_unique<GdalReader>(files[0]));
    } else {
        for (const auto& f : files)
            readers.push_back(std::make_unique<GdalReader>(f));
    }
    return readers;
}

// ── Tile-slot / pipeline helpers ────────────────────────────────────────────

// One reusable tile buffer set.  The vectors are resized per tile, but each
// slot is allocated once to the maximum tile size to avoid repeated heap calls.
struct TileSlot {
    int xoff = 0;
    int yoff = 0;
    int xsize = 0;
    int ysize = 0;
    int n_pix = 0;

    std::vector<double> tas;
    std::vector<double> tasmax;
    std::vector<double> tasmin;
    std::vector<double> pr;
    std::vector<double> mask;
    std::vector<double> bio;

    void reserve(int max_n_pix) {
        tas.reserve(static_cast<std::size_t>(12) * max_n_pix);
        tasmax.reserve(static_cast<std::size_t>(12) * max_n_pix);
        tasmin.reserve(static_cast<std::size_t>(12) * max_n_pix);
        pr.reserve(static_cast<std::size_t>(12) * max_n_pix);
        mask.reserve(static_cast<std::size_t>(max_n_pix));
        bio.reserve(static_cast<std::size_t>(19) * max_n_pix);
    }

    void resize(int np) {
        if (n_pix == np) return;
        n_pix = np;
        tas.resize(static_cast<std::size_t>(12) * n_pix);
        tasmax.resize(static_cast<std::size_t>(12) * n_pix);
        tasmin.resize(static_cast<std::size_t>(12) * n_pix);
        pr.resize(static_cast<std::size_t>(12) * n_pix);
        // mask is resized by the reader only when a mask is configured.
        bio.resize(static_cast<std::size_t>(19) * n_pix);
    }
};

// Simple mutex/CV queue used to hand off TileSlot* pointers between the
// reader, compute, and writer threads.  A nullptr is used as an end-of-stream
// sentinel.  Once cancel() has been called, pop() drains any remaining items
// and then returns nullptr, allowing workers to exit cleanly after an error.
class TileQueue {
public:
    void push(TileSlot* slot) {
        {
            std::lock_guard<std::mutex> lk(mtx_);
            q_.push_back(slot);
        }
        cv_.notify_one();
    }

    TileSlot* pop() {
        std::unique_lock<std::mutex> lk(mtx_);
        cv_.wait(lk, [this] { return !q_.empty() || cancelled_; });
        if (q_.empty()) return nullptr;
        TileSlot* slot = q_.front();
        q_.pop_front();
        return slot;
    }

    void cancel() {
        {
            std::lock_guard<std::mutex> lk(mtx_);
            cancelled_ = true;
        }
        cv_.notify_all();
    }

private:
    std::deque<TileSlot*> q_;
    std::mutex mtx_;
    std::condition_variable cv_;
    bool cancelled_ = false;
};

// Shared error/cancellation state for the three pipeline workers.
struct PipelineContext {
    std::atomic<bool> has_error{false};
    std::string error_message;
    std::mutex error_mtx;

    TileQueue free_queue;
    TileQueue ready_compute;
    TileQueue ready_write;

    void report_error(const std::string& msg) {
        std::lock_guard<std::mutex> lk(error_mtx);
        if (!has_error.load()) {
            has_error.store(true);
            error_message = msg;
        }
        free_queue.cancel();
        ready_compute.cancel();
        ready_write.cancel();
    }
};

// Read one tile into a slot.  This function is always called from the reader
// thread so no other thread touches the GdalReader objects used here.
void read_tile_into_slot(
    TileSlot* slot,
    const std::vector<std::unique_ptr<GdalReader>>& tas_readers,
    const std::vector<std::unique_ptr<GdalReader>>& tasmax_readers,
    const std::vector<std::unique_ptr<GdalReader>>& tasmin_readers,
    const std::vector<std::unique_ptr<GdalReader>>& pr_readers,
    bool tas_multi, bool tasmax_multi, bool tasmin_multi, bool pr_multi,
    const GdalReader* mask_reader,
    const std::vector<int>& input_bands,
    int xoff, int yoff, int xsize, int ysize) {

    const int n_pix = xsize * ysize;
    slot->xoff = xoff;
    slot->yoff = yoff;
    slot->xsize = xsize;
    slot->ysize = ysize;
    slot->resize(n_pix);

    if (tas_multi) {
        tas_readers[0]->read_bands_window(
            xoff, yoff, xsize, ysize, input_bands, slot->tas);
    } else {
        for (int m = 0; m < 12; ++m) {
            tas_readers[static_cast<std::size_t>(m)]->read_window(
                xoff, yoff, xsize, ysize, 1,
                slot->tas.data() + static_cast<std::size_t>(m) * n_pix);
        }
    }

    if (tasmax_multi) {
        tasmax_readers[0]->read_bands_window(
            xoff, yoff, xsize, ysize, input_bands, slot->tasmax);
    } else {
        for (int m = 0; m < 12; ++m) {
            tasmax_readers[static_cast<std::size_t>(m)]->read_window(
                xoff, yoff, xsize, ysize, 1,
                slot->tasmax.data() + static_cast<std::size_t>(m) * n_pix);
        }
    }

    if (tasmin_multi) {
        tasmin_readers[0]->read_bands_window(
            xoff, yoff, xsize, ysize, input_bands, slot->tasmin);
    } else {
        for (int m = 0; m < 12; ++m) {
            tasmin_readers[static_cast<std::size_t>(m)]->read_window(
                xoff, yoff, xsize, ysize, 1,
                slot->tasmin.data() + static_cast<std::size_t>(m) * n_pix);
        }
    }

    if (pr_multi) {
        pr_readers[0]->read_bands_window(
            xoff, yoff, xsize, ysize, input_bands, slot->pr);
    } else {
        for (int m = 0; m < 12; ++m) {
            pr_readers[static_cast<std::size_t>(m)]->read_window(
                xoff, yoff, xsize, ysize, 1,
                slot->pr.data() + static_cast<std::size_t>(m) * n_pix);
        }
    }

    if (mask_reader) {
        slot->mask.resize(static_cast<std::size_t>(n_pix));
        mask_reader->read_window(xoff, yoff, xsize, ysize, 1, slot->mask);
    } else {
        slot->mask.clear();
    }
}


// Scalar tile driver (used when AVX2 is unavailable).
void compute_tile_buffers_scalar(const double* tas, const double* tasmax,
                                 const double* tasmin, const double* pr,
                                 const double* mask, double* bio,
                                 int n_pix) {
    for (int i = 0; i < n_pix; ++i) {
        compute_pixel_at(tas, tasmax, tasmin, pr, mask, bio, n_pix, i);
    }
}

#ifdef __AVX2__
#include <immintrin.h>

// Compute four contiguous pixels in one AVX2 pass.
// Data layout: buf[month * n_pix + pixel], so pixels i0..i0+3 for a fixed
// month are contiguous and can be loaded into one __m256d.
inline void compute_pixel4(const double* tas, const double* tasmax,
                           const double* tasmin, const double* pr,
                           const double* mask, double* bio,
                           int n_pix, int i0) {
    const __m256d zero = _mm256_setzero_pd();
    const __m256d inv12 = _mm256_set1_pd(1.0 / 12.0);
    const __m256d inv3  = _mm256_set1_pd(1.0 / 3.0);
    const __m256d hundred = _mm256_set1_pd(100.0);
    const __m256d nan = _mm256_set1_pd(
        std::numeric_limits<double>::quiet_NaN());
    // -1.0 has the sign bit set, which is all blendv_pd needs for "true".
    __m256d valid = _mm256_set1_pd(-1.0);

    // Month 0: initialise accumulators and the 12-month caches.
    __m256d t0   = _mm256_loadu_pd(tas    + static_cast<std::size_t>(i0));
    __m256d tmx0 = _mm256_loadu_pd(tasmax + static_cast<std::size_t>(i0));
    __m256d tmn0 = _mm256_loadu_pd(tasmin + static_cast<std::size_t>(i0));
    __m256d p0   = _mm256_loadu_pd(pr     + static_cast<std::size_t>(i0));

    __m256d month_t[12];
    __m256d month_p[12];
    month_t[0] = t0;
    month_p[0] = p0;

    // NaN guard for month 0.
    __m256d nan_t   = _mm256_cmp_pd(t0,   t0,   _CMP_UNORD_Q);
    __m256d nan_tmx = _mm256_cmp_pd(tmx0, tmx0, _CMP_UNORD_Q);
    __m256d nan_tmn = _mm256_cmp_pd(tmn0, tmn0, _CMP_UNORD_Q);
    __m256d nan_p   = _mm256_cmp_pd(p0,   p0,   _CMP_UNORD_Q);
    __m256d any_nan = _mm256_or_pd(_mm256_or_pd(nan_t, nan_tmx),
                                   _mm256_or_pd(nan_tmn, nan_p));
    valid = _mm256_andnot_pd(any_nan, valid);

    __m256d t_sum   = t0;
    __m256d diurnal = _mm256_sub_pd(tmx0, tmn0);
    __m256d tmx_max = tmx0;
    __m256d tmn_min = tmn0;
    __m256d p_sum   = p0;
    __m256d p_max   = p0;
    __m256d p_min   = p0;

    // Months 1..11.
    for (int m = 1; m < 12; ++m) {
        const std::size_t off = static_cast<std::size_t>(m) * n_pix + i0;
        __m256d t   = _mm256_loadu_pd(tas    + off);
        __m256d tmx = _mm256_loadu_pd(tasmax + off);
        __m256d tmn = _mm256_loadu_pd(tasmin + off);
        __m256d p   = _mm256_loadu_pd(pr     + off);

        month_t[m] = t;
        month_p[m] = p;

        nan_t   = _mm256_cmp_pd(t,   t,   _CMP_UNORD_Q);
        nan_tmx = _mm256_cmp_pd(tmx, tmx, _CMP_UNORD_Q);
        nan_tmn = _mm256_cmp_pd(tmn, tmn, _CMP_UNORD_Q);
        nan_p   = _mm256_cmp_pd(p,   p,   _CMP_UNORD_Q);
        any_nan = _mm256_or_pd(_mm256_or_pd(nan_t, nan_tmx),
                               _mm256_or_pd(nan_tmn, nan_p));
        valid = _mm256_andnot_pd(any_nan, valid);

        t_sum   = _mm256_add_pd(t_sum, t);
        diurnal = _mm256_add_pd(diurnal, _mm256_sub_pd(tmx, tmn));
        tmx_max = _mm256_max_pd(tmx_max, tmx);
        tmn_min = _mm256_min_pd(tmn_min, tmn);
        p_sum   = _mm256_add_pd(p_sum, p);
        p_max   = _mm256_max_pd(p_max, p);
        p_min   = _mm256_min_pd(p_min, p);
    }

    if (mask != nullptr) {
        __m256d m = _mm256_loadu_pd(mask + i0);
        __m256d nan_m  = _mm256_cmp_pd(m, m, _CMP_UNORD_Q);
        __m256d zero_m = _mm256_cmp_pd(m, zero, _CMP_EQ_OQ);
        __m256d bad = _mm256_or_pd(nan_m, zero_m);
        valid = _mm256_andnot_pd(bad, valid);
    }

    // BIO01--BIO07.
    __m256d b01 = _mm256_mul_pd(t_sum, inv12);
    __m256d b02 = _mm256_mul_pd(diurnal, inv12);
    __m256d b05 = tmx_max;
    __m256d b06 = tmn_min;
    __m256d b07 = _mm256_sub_pd(b05, b06);

    __m256d b07_gt0 = _mm256_cmp_pd(b07, zero, _CMP_GT_OQ);
    __m256d b03 = _mm256_mul_pd(hundred, _mm256_div_pd(b02, b07));
    b03 = _mm256_blendv_pd(zero, b03, b07_gt0);

    // BIO04: use a two-pass variance (matches the R bioclim() reference).
    __m256d t_mean = b01;
    __m256d t_ssq  = zero;
    for (int m = 0; m < 12; ++m) {
        __m256d d = _mm256_sub_pd(month_t[m], t_mean);
        t_ssq = _mm256_add_pd(t_ssq, _mm256_mul_pd(d, d));
    }
    __m256d t_var = _mm256_mul_pd(t_ssq, inv12);
    t_var = _mm256_max_pd(t_var, zero);
    __m256d b04 = _mm256_mul_pd(hundred, _mm256_sqrt_pd(t_var));

    // BIO12--BIO15.
    __m256d b12 = p_sum;
    __m256d b13 = p_max;
    __m256d b14 = p_min;

    __m256d p_mean = _mm256_mul_pd(p_sum, inv12);
    __m256d p_ssq = zero;
    for (int m = 0; m < 12; ++m) {
        __m256d d = _mm256_sub_pd(month_p[m], p_mean);
        p_ssq = _mm256_add_pd(p_ssq, _mm256_mul_pd(d, d));
    }
    __m256d p_var = _mm256_mul_pd(p_ssq, inv12);
    p_var = _mm256_max_pd(p_var, zero);
    __m256d b15 = _mm256_mul_pd(hundred,
                                _mm256_div_pd(_mm256_sqrt_pd(p_var), p_mean));
    __m256d p_mean_eq0 = _mm256_cmp_pd(p_mean, zero, _CMP_EQ_OQ);
    b15 = _mm256_blendv_pd(b15, nan, p_mean_eq0);

    // Rolling quarter sums, re-using the cached month vectors.
    __m256d t_qs[12];
    __m256d p_qs[12];
    for (int k = 0; k < 12; ++k) {
        int k1 = k + 1;
        if (k1 == 12) k1 = 0;
        int k2 = k + 2;
        if (k2 == 12) k2 = 0;
        if (k2 == 13) k2 = 1;

        t_qs[k] = _mm256_add_pd(month_t[k],
                     _mm256_add_pd(month_t[k1], month_t[k2]));
        p_qs[k] = _mm256_add_pd(month_p[k],
                     _mm256_add_pd(month_p[k1], month_p[k2]));
    }

    // Argmax / argmin over the 12 quarter-sum vectors, keeping both the
    // temperature and precipitation values at the chosen index.
    __m256d wet_pr  = p_qs[0], dry_pr  = p_qs[0];
    __m256d warm_t  = t_qs[0], cold_t  = t_qs[0];
    __m256d wet_t   = t_qs[0], dry_t   = t_qs[0];
    __m256d warm_pr = p_qs[0], cold_pr = p_qs[0];

    for (int k = 1; k < 12; ++k) {
        __m256d pk = p_qs[k];
        __m256d tk = t_qs[k];

        __m256d gt_pr = _mm256_cmp_pd(pk, wet_pr, _CMP_GT_OQ);
        wet_pr = _mm256_blendv_pd(wet_pr, pk, gt_pr);
        wet_t  = _mm256_blendv_pd(wet_t,  tk, gt_pr);

        __m256d lt_pr = _mm256_cmp_pd(pk, dry_pr, _CMP_LT_OQ);
        dry_pr = _mm256_blendv_pd(dry_pr, pk, lt_pr);
        dry_t  = _mm256_blendv_pd(dry_t,  tk, lt_pr);

        __m256d gt_t = _mm256_cmp_pd(tk, warm_t, _CMP_GT_OQ);
        warm_t  = _mm256_blendv_pd(warm_t,  tk, gt_t);
        warm_pr = _mm256_blendv_pd(warm_pr, pk, gt_t);

        __m256d lt_t = _mm256_cmp_pd(tk, cold_t, _CMP_LT_OQ);
        cold_t  = _mm256_blendv_pd(cold_t,  tk, lt_t);
        cold_pr = _mm256_blendv_pd(cold_pr, pk, lt_t);
    }

    __m256d b08 = _mm256_mul_pd(wet_t,  inv3);
    __m256d b09 = _mm256_mul_pd(dry_t,  inv3);
    __m256d b10 = _mm256_mul_pd(warm_t, inv3);
    __m256d b11 = _mm256_mul_pd(cold_t, inv3);

    __m256d b16 = wet_pr;
    __m256d b17 = dry_pr;
    __m256d b18 = warm_pr;
    __m256d b19 = cold_pr;

    // Apply the validity mask and store 19 outputs for 4 pixels.
    __m256d* out[19] = {
        &b01, &b02, &b03, &b04, &b05, &b06, &b07, &b08, &b09,
        &b10, &b11, &b12, &b13, &b14, &b15, &b16, &b17, &b18, &b19
    };
    for (int j = 0; j < 19; ++j) {
        __m256d v = _mm256_blendv_pd(nan, *out[j], valid);
        _mm256_storeu_pd(bio + static_cast<std::size_t>(j) * n_pix + i0, v);
    }
}
#endif  // __AVX2__

// Compute one tile on the CPU.
void compute_tile_buffers(const double* tas, const double* tasmax,
                          const double* tasmin, const double* pr,
                          const double* mask, double* bio,
                          int n_pix, int n_threads) {
#ifdef __AVX2__
    const int full_end = n_pix - 3;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
    for (int i = 0; i < full_end; i += 4) {
        compute_pixel4(tas, tasmax, tasmin, pr, mask, bio, n_pix, i);
    }

    // Scalar tail for the last 0-3 pixels.
    const int tail = n_pix - (n_pix % 4);
    for (int k = tail; k < n_pix; ++k) {
        compute_pixel_at(tas, tasmax, tasmin, pr, mask, bio, n_pix, k);
    }
#else
    compute_tile_buffers_scalar(tas, tasmax, tasmin, pr, mask, bio, n_pix);
    (void)n_threads;  // unused when scalar
#endif
}

// TileSlot wrapper around the raw-buffer compute kernel.
void compute_tile_cpu(TileSlot* slot, int n_threads) {
    const double* mask = slot->mask.empty() ? nullptr : slot->mask.data();
    compute_tile_buffers(slot->tas.data(), slot->tasmax.data(),
                         slot->tasmin.data(), slot->pr.data(),
                         mask, slot->bio.data(),
                         slot->n_pix, n_threads);
}

}  // anonymous namespace

#endif  // HAVE_GDAL

// ── BioclimEngine::compute() ─────────────────────────────────────────────────

std::string BioclimEngine::compute() {
#ifdef HAVE_GDAL
    // ── Validate configuration ──────────────────────────────────────────────
    if (output_path_.empty())
        throw std::runtime_error("BioclimEngine::compute: output path not set.");

    validate_file_vector(tas_files_,    "tas");
    validate_file_vector(tasmax_files_, "tasmax");
    validate_file_vector(tasmin_files_, "tasmin");
    validate_file_vector(pr_files_,     "pr");

    // Use the overlapped pipeline if requested.
    if (pipeline_) {
        return compute_pipelined();
    }

    // ── Determine effective variable set ─────────────────────────────────────
    std::vector<int> eff_vars = variables_;
    if (eff_vars.empty()) {
        eff_vars.resize(19);
        for (int i = 0; i < 19; ++i) eff_vars[i] = i + 1;
    }

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

    // ── Create one 19-band output GeoTIFF ───────────────────────────────────
    // output_path_ is a directory; all variables are written to a single
    // multi-band GeoTIFF named bio.tif.
    std::ostringstream out_fname;
    out_fname << output_path_ << "/bio.tif";
    GdalWriter writer(out_fname.str(), nrows, ncols, 19, gt, crs,
                      false, output_dtype_);

    // Band map for writing all 19 bands in one call (1-based).
    std::vector<int> write_bands(19);
    std::iota(write_bands.begin(), write_bands.end(), 1);

    // ── Determine compute device ────────────────────────────────────────────
#ifdef HAVE_CUDA
    bool use_gpu = false;
    if (device_ != Device::CPU) {
        int gpu_count = 0;
        cudaError_t cuda_err = cudaGetDeviceCount(&gpu_count);
        if (cuda_err == cudaSuccess && gpu_count > 0) {
            use_gpu = true;
        }
    }
#endif

    // ── Tiled processing loop ───────────────────────────────────────────────
    const int ts = tile_size_;
    const int max_n_pix = ts * ts;

    // Reusable tile buffers; capacity is the largest possible tile.
    std::vector<double> tas_tile;
    std::vector<double> tasmax_tile;
    std::vector<double> tasmin_tile;
    std::vector<double> pr_tile;
    std::vector<double> bio_tile;
    std::vector<double> mask_buf;

    tas_tile.reserve(static_cast<std::size_t>(12) * max_n_pix);
    tasmax_tile.reserve(static_cast<std::size_t>(12) * max_n_pix);
    tasmin_tile.reserve(static_cast<std::size_t>(12) * max_n_pix);
    pr_tile.reserve(static_cast<std::size_t>(12) * max_n_pix);
    bio_tile.reserve(static_cast<std::size_t>(19) * max_n_pix);
    mask_buf.reserve(static_cast<std::size_t>(max_n_pix));

    for (int yoff = 0; yoff < nrows; yoff += ts) {
        const int ysize = std::min(ts, nrows - yoff);
        for (int xoff = 0; xoff < ncols; xoff += ts) {
            const int xsize = std::min(ts, ncols - xoff);
            const int n_pix = xsize * ysize;

            // ── Resize / reuse the tile buffers for this tile ───────────────
            const std::size_t n_in  = static_cast<std::size_t>(12) * n_pix;
            const std::size_t n_out = static_cast<std::size_t>(19) * n_pix;
            if (tas_tile.size() != n_in)     tas_tile.resize(n_in);
            if (tasmax_tile.size() != n_in)  tasmax_tile.resize(n_in);
            if (tasmin_tile.size() != n_in)  tasmin_tile.resize(n_in);
            if (pr_tile.size() != n_in)      pr_tile.resize(n_in);
            if (bio_tile.size() != n_out)    bio_tile.resize(n_out);

            // ── Read 4 × 12 monthly bands into band-major tile buffers ─────
            // Layout: tile[var][month * n_pix + pixel] (band-major).
            // The 48 band reads run in parallel: each month has its own
            // GdalReader (GDAL datasets are not thread-safe) and writes
            // directly into the contiguous month slice.
            std::vector<std::unique_ptr<GdalReader>>* var_readers[4] = {
                &tas_readers, &tasmax_readers, &tasmin_readers, &pr_readers
            };
            const bool var_multi[4] = {
                tas_multi, tasmax_multi, tasmin_multi, pr_multi
            };
            double* var_tiles[4] = {
                tas_tile.data(), tasmax_tile.data(),
                tasmin_tile.data(), pr_tile.data()
            };

            std::exception_ptr read_error;
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static) num_threads(n_threads_)
#endif
            for (int v = 0; v < 4; ++v) {
                for (int m = 0; m < 12; ++m) {
                    try {
                        (*var_readers[v])[static_cast<std::size_t>(m)]
                            ->read_window(
                                xoff, yoff, xsize, ysize,
                                var_multi[v] ? m + 1 : 1,
                                var_tiles[v] +
                                    static_cast<std::size_t>(m) * n_pix);
                    } catch (...) {
                        // Exceptions must not escape an OpenMP region;
                        // record the first one and rethrow below.
#ifdef _OPENMP
#pragma omp critical
#endif
                        {
                            if (!read_error)
                                read_error = std::current_exception();
                        }
                    }
                }
            }
            if (read_error) std::rethrow_exception(read_error);

            // ── Read mask tile (optional) ────────────────────────────────────
            if (mask_reader) {
                mask_buf.resize(static_cast<std::size_t>(n_pix));
                mask_reader->read_window(xoff, yoff, xsize, ysize, 1, mask_buf);
            } else {
                mask_buf.clear();
            }

            // ── Compute 19 bio variables per pixel ───────────────────────────
            // Output layout: bio_tile[bio_index * n_pix + pixel_index]
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
                compute_tile_buffers(
                    tas_tile.data(), tasmax_tile.data(),
                    tasmin_tile.data(), pr_tile.data(),
                    mask_buf.empty() ? nullptr : mask_buf.data(),
                    bio_tile.data(),
                    n_pix, n_threads_);
            }  // end CPU branch

            // ── Write all 19 output bands in one multi-band call ─────────────
            writer.write_bands_window(xoff, yoff, xsize, ysize,
                                      write_bands, bio_tile, output_dtype_);
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

// ── BioclimEngine::compute_pipelined() ───────────────────────────────────────

std::string BioclimEngine::compute_pipelined() {
#ifdef HAVE_GDAL
    // ── Validate configuration ──────────────────────────────────────────────
    if (output_path_.empty())
        throw std::runtime_error("BioclimEngine::compute_pipelined: output path not set.");

    validate_file_vector(tas_files_,    "tas");
    validate_file_vector(tasmax_files_, "tasmax");
    validate_file_vector(tasmin_files_, "tasmin");
    validate_file_vector(pr_files_,     "pr");

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

    // ── Create one 19-band output GeoTIFF ───────────────────────────────────
    std::ostringstream out_fname;
    out_fname << output_path_ << "/bio.tif";
    GdalWriter writer(out_fname.str(), nrows, ncols, 19, gt, crs,
                      false, output_dtype_);

    // Band map for writing all 19 bands in one call (1-based).
    std::vector<int> write_bands(19);
    std::iota(write_bands.begin(), write_bands.end(), 1);

    // Input band map for multi-band reads (1..12).
    std::vector<int> input_bands(12);
    std::iota(input_bands.begin(), input_bands.end(), 1);

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

    // ── Set up the overlapped pipeline ──────────────────────────────────────
    const int ts = tile_size_;
    const int max_n_pix = ts * ts;

    PipelineContext ctx;
    std::vector<std::unique_ptr<TileSlot>> slots;
    slots.reserve(3);
    for (int i = 0; i < 3; ++i) {
        auto s = std::make_unique<TileSlot>();
        s->reserve(max_n_pix);
        slots.push_back(std::move(s));
        ctx.free_queue.push(slots.back().get());
    }

    // Reader thread: fills free slots and pushes them to the compute queue.
    std::thread reader_thread([&]() {
        try {
            for (int yoff = 0; yoff < nrows; yoff += ts) {
                const int ysize = std::min(ts, nrows - yoff);
                for (int xoff = 0; xoff < ncols; xoff += ts) {
                    const int xsize = std::min(ts, ncols - xoff);

                    TileSlot* slot = ctx.free_queue.pop();
                    if (slot == nullptr) break;

                    read_tile_into_slot(
                        slot,
                        tas_readers, tasmax_readers, tasmin_readers, pr_readers,
                        tas_multi, tasmax_multi, tasmin_multi, pr_multi,
                        mask_reader.get(), input_bands,
                        xoff, yoff, xsize, ysize);

                    ctx.ready_compute.push(slot);
                }
            }
            ctx.ready_compute.push(nullptr);  // end-of-stream sentinel
        } catch (const std::exception& e) {
            ctx.report_error(e.what());
        }
    });

    // Compute thread: consumes filled slots, runs the per-pixel BIOCLIM loop,
    // and pushes results to the writer queue.
    std::thread compute_thread([&]() {
        try {
            while (true) {
                TileSlot* slot = ctx.ready_compute.pop();
                if (slot == nullptr) {
                    ctx.ready_write.push(nullptr);
                    break;
                }

#ifdef HAVE_CUDA
                if (use_gpu) {
                    launch_bioclim_cuda(
                        slot->tas.data(), slot->tasmax.data(),
                        slot->tasmin.data(), slot->pr.data(),
                        slot->mask.empty() ? nullptr : slot->mask.data(),
                        slot->bio.data(),
                        slot->n_pix
                    );
                } else
#endif
                {
                    compute_tile_cpu(slot, n_threads_);
                }

                ctx.ready_write.push(slot);
            }
        } catch (const std::exception& e) {
            ctx.report_error(e.what());
        }
    });

    // Writer thread: drains computed tiles in order and returns slots to the
    // free queue for reuse.
    std::thread writer_thread([&]() {
        try {
            while (true) {
                TileSlot* slot = ctx.ready_write.pop();
                if (slot == nullptr) break;

                writer.write_bands_window(
                    slot->xoff, slot->yoff, slot->xsize, slot->ysize,
                    write_bands, slot->bio, output_dtype_);

                ctx.free_queue.push(slot);
            }
        } catch (const std::exception& e) {
            ctx.report_error(e.what());
        }
    });

    reader_thread.join();
    compute_thread.join();
    writer_thread.join();

    if (ctx.has_error.load()) {
        throw std::runtime_error(ctx.error_message);
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

}  // namespace xbioclim

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
    Rcpp::XPtr<xbioclim::BioclimEngine> ptr(
        new xbioclim::BioclimEngine(), true);
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
    Rcpp::XPtr<xbioclim::BioclimEngine> eng(xptr);
    eng->open(Rcpp::as<std::vector<std::string>>(tas_files),
              Rcpp::as<std::vector<std::string>>(tasmax_files),
              Rcpp::as<std::vector<std::string>>(tasmin_files),
              Rcpp::as<std::vector<std::string>>(pr_files));
}

//' Set the output raster path
//'
//' The engine will create (or overwrite) a multi-band GeoTIFF named
//' \code{bio.tif} inside this directory when \code{\link{engine_compute}} is
//' called.
//'
//' @param xptr External pointer returned by \code{\link{engine_create}}.
//' @param path Character scalar: output directory path.
//' @return \code{NULL} invisibly.
//' @seealso \code{\link{engine_create}}, \code{\link{engine_compute}}
// [[Rcpp::export]]
void engine_set_output(SEXP xptr, std::string path) {
    Rcpp::XPtr<xbioclim::BioclimEngine>(xptr)->set_output(path);
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
    Rcpp::XPtr<xbioclim::BioclimEngine>(xptr)->set_mask(mask_path);
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
    Rcpp::XPtr<xbioclim::BioclimEngine>(xptr)->set_threads(n);
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
    Rcpp::XPtr<xbioclim::BioclimEngine>(xptr)->set_tile_size(tile_size);
}

//' Set the output data type
//'
//' Controls the on-disk data type of the output \code{bio.tif} file.
//'\describe{
//'   \item{"Float64"}{IEEE 754 double precision (default).}
//'   \item{"Float32"}{IEEE 754 single precision — half the file size with
//'     negligible loss for most climate data.}
//' }
//'
//' @param xptr  External pointer returned by \code{\link{engine_create}}.
//' @param dtype Character scalar: one of \code{"Float64"} or \code{"Float32"}.
//' @return \code{NULL} invisibly.
//' @seealso \code{\link{engine_create}}, \code{\link{engine_compute}}
//' @keywords internal
// [[Rcpp::export]]
void engine_set_dtype(SEXP xptr, std::string dtype) {
    Rcpp::XPtr<xbioclim::BioclimEngine>(xptr)->set_dtype(dtype);
}

//' Select which bioclimatic variables to write
//'
//' Restricts the output to a subset of the 19 standard bioclimatic variables.
//' The engine always computes all 19 internally (they share intermediate
//' values).  With the multi-band output file, all 19 bands are written and
//' the \code{\link{bioclim_engine}} R wrapper subsets the returned
//' \code{SpatRaster}.
//'
//' @param xptr      External pointer returned by \code{\link{engine_create}}.
//' @param variables Integer vector with elements in 1..19.
//' @return \code{NULL} invisibly.
//' @seealso \code{\link{engine_create}}, \code{\link{engine_compute}}
//' @export
// [[Rcpp::export]]
void engine_set_variables(SEXP xptr, Rcpp::IntegerVector variables) {
    Rcpp::XPtr<xbioclim::BioclimEngine> eng(xptr);
    eng->set_variables(Rcpp::as<std::vector<int>>(variables));
}

//' Enable or disable the overlapped read/compute/write pipeline
//'
//' This is an internal, opt-in flag.  When \code{TRUE}, the next call to
//' \code{\link{engine_compute}} uses three background threads to overlap
//' the GDAL read, BIOCLIM computation, and GDAL write stages for each tile.
//' When \code{FALSE} (the default) the engine uses the original serial loop.
//'
//' @param xptr External pointer returned by \code{\link{engine_create}}.
//' @param use_pipeline Logical scalar: \code{TRUE} to enable the pipeline.
//' @return \code{NULL} invisibly.
//' @seealso \code{\link{engine_create}}, \code{\link{engine_compute}}
//' @keywords internal
// [[Rcpp::export]]
void engine_set_pipeline(SEXP xptr, bool use_pipeline) {
    Rcpp::XPtr<xbioclim::BioclimEngine> eng(xptr);
    eng->set_pipeline(use_pipeline);
}

//' Run the bioclimatic-variable computation pipeline
//'
//' Reads all monthly climate input rasters tile by tile, computes the
//' bioclimatic variables for every pixel, and writes all 19 variables to a
//' single multi-band GeoTIFF named \code{bio.tif} inside the output
//' directory.  Peak memory is proportional to the tile size, not the full
//' raster size.
//'
//' Requires GDAL support.  Stops with an informative error when the package
//' was built without GDAL.
//'
//' @param xptr External pointer returned by \code{\link{engine_create}}.
//' @return Character scalar: the output directory path (same as the value
//'   passed to \code{\link{engine_set_output}}).
//' @seealso \code{\link{engine_create}}, \code{\link{engine_set_output}},
//'   \code{\link{has_gdal}}
// [[Rcpp::export]]
std::string engine_compute(SEXP xptr) {
    Rcpp::XPtr<xbioclim::BioclimEngine> eng(xptr);
    return eng->compute();
}
