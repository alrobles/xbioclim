// bioclim_cuda.cpp — GPU kernel for bioclimatic variable computation.
//
// This file is compiled with nvcc -x cu only when HAVE_CUDA is defined
// (detected by configure.ac).  All CUDA-specific code is guarded by
// #ifdef HAVE_CUDA so that the file compiles to nothing when CUDA is absent.
//
// Data layout (pixel-major, matches BioclimEngine.cpp tile buffers):
//   Input:  var[i * 12 + month]    — 12 months × n_pix pixels
//   Mask:   mask[i]                — n_pix pixels (nullptr = no mask)
//   Output: bio[i * 19 + bio_idx]  — 19 variables × n_pix pixels

#ifdef HAVE_CUDA

#include <cuda_runtime.h>
#include <cmath>
#include <stdexcept>
#include <string>

// Quiet NaN bit pattern (IEEE 754 double precision)
static constexpr unsigned long long kQuietNanBits = 0x7FF8000000000000ULL;

// Macro for CUDA error checking: throws std::runtime_error on failure.
#define CUDA_CHECK(call)                                                       \
    do {                                                                       \
        cudaError_t _e = (call);                                               \
        if (_e != cudaSuccess) {                                               \
            throw std::runtime_error(                                          \
                std::string("CUDA error in " #call ": ") +                    \
                cudaGetErrorString(_e));                                        \
        }                                                                      \
    } while (0)

// ── Device helper functions ───────────────────────────────────────────────────

// Rolling 3-month circular sums: qs[k] = x[k] + x[(k+1)%12] + x[(k+2)%12]
__device__ static void rolling_quarter_sum_dev(const double* x, double* qs) {
    for (int k = 0; k < 12; ++k)
        qs[k] = x[k] + x[(k + 1) % 12] + x[(k + 2) % 12];
}

// Index of maximum in a 12-element array (0-based)
__device__ static int argmax12_dev(const double* x) {
    int best = 0;
    for (int k = 1; k < 12; ++k) if (x[k] > x[best]) best = k;
    return best;
}

// Index of minimum in a 12-element array (0-based)
__device__ static int argmin12_dev(const double* x) {
    int best = 0;
    for (int k = 1; k < 12; ++k) if (x[k] < x[best]) best = k;
    return best;
}

// Population standard deviation (denominator N) for a 12-element array
__device__ static double sd_pop_dev(const double* x) {
    double s = 0.0, ss = 0.0;
    for (int i = 0; i < 12; ++i) { s += x[i]; ss += x[i] * x[i]; }
    const double m = s / 12.0;
    return sqrt(ss / 12.0 - m * m);
}

// ── Main kernel ───────────────────────────────────────────────────────────────

// Each CUDA thread processes one pixel (index i = blockIdx.x * blockDim.x + threadIdx.x).
// The kernel mirrors the logic of compute_pixel() in BioclimEngine.cpp exactly.
__global__ void bioclim_kernel(
    const double* __restrict__ tas,     // [12 * n_pix]
    const double* __restrict__ tasmax,  // [12 * n_pix]
    const double* __restrict__ tasmin,  // [12 * n_pix]
    const double* __restrict__ pr,      // [12 * n_pix]
    const double* __restrict__ mask,    // [n_pix] or nullptr
    double* __restrict__       bio,     // [19 * n_pix]
    int n_pix
) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_pix) return;

    // -- Apply mask (0 or NaN → all-NaN output) --------------------------------
    if (mask != nullptr) {
        const double mv = mask[i];
        if (isnan(mv) || mv == 0.0) {
            for (int j = 0; j < 19; ++j)
                bio[i * 19 + j] = __longlong_as_double(kQuietNanBits);
            return;
        }
    }

    // -- Load 12-month vectors for this pixel (pixel-major layout) -------------
    double t[12], tmx[12], tmn[12], p[12];
    for (int m = 0; m < 12; ++m) {
        t[m]   = tas[   i * 12 + m];
        tmx[m] = tasmax[i * 12 + m];
        tmn[m] = tasmin[i * 12 + m];
        p[m]   = pr[    i * 12 + m];
    }

    // -- NA / NaN guard -------------------------------------------------------
    for (int m = 0; m < 12; ++m) {
        if (isnan(t[m]) || isnan(tmx[m]) || isnan(tmn[m]) || isnan(p[m])) {
            for (int j = 0; j < 19; ++j)
                bio[i * 19 + j] = __longlong_as_double(kQuietNanBits);
            return;
        }
    }

    // -- Temperature aggregates -----------------------------------------------
    double t_sum = 0.0, diurnal_sum = 0.0;
    double tmx_max = tmx[0], tmn_min = tmn[0];
    for (int m = 0; m < 12; ++m) {
        t_sum       += t[m];
        diurnal_sum += tmx[m] - tmn[m];
        if (tmx[m] > tmx_max) tmx_max = tmx[m];
        if (tmn[m] < tmn_min) tmn_min = tmn[m];
    }

    const double b01 = t_sum / 12.0;
    const double b02 = diurnal_sum / 12.0;
    const double b05 = tmx_max;
    const double b06 = tmn_min;
    const double b07 = b05 - b06;
    const double b03 = (b07 > 0.0) ? 100.0 * b02 / b07 : 0.0;
    const double b04 = 100.0 * sd_pop_dev(t);

    // -- Precipitation aggregates ---------------------------------------------
    double p_sum = 0.0, p_max = p[0], p_min = p[0];
    for (int m = 0; m < 12; ++m) {
        p_sum += p[m];
        if (p[m] > p_max) p_max = p[m];
        if (p[m] < p_min) p_min = p[m];
    }
    const double b12 = p_sum;
    const double b13 = p_max;
    const double b14 = p_min;

    // BIO15: Precipitation Seasonality (CV = 100 * sd / mean; NaN when mean == 0)
    const double p_mean = p_sum / 12.0;
    double p_ssq = 0.0;
    for (int m = 0; m < 12; ++m) {
        const double d = p[m] - p_mean;
        p_ssq += d * d;
    }
    const double b15 = (p_mean == 0.0)
        ? __longlong_as_double(kQuietNanBits)
        : 100.0 * sqrt(p_ssq / 12.0) / p_mean;

    // -- Rolling quarter sums -------------------------------------------------
    double pr_qs[12], t_qs[12];
    rolling_quarter_sum_dev(p, pr_qs);
    rolling_quarter_sum_dev(t, t_qs);

    const int wet_q  = argmax12_dev(pr_qs);
    const int dry_q  = argmin12_dev(pr_qs);
    const int warm_q = argmax12_dev(t_qs);
    const int cold_q = argmin12_dev(t_qs);

    // -- Write output (pixel-major: bio[i * 19 + bio_idx]) ---------------------
    bio[i * 19 +  0] = b01;
    bio[i * 19 +  1] = b02;
    bio[i * 19 +  2] = b03;
    bio[i * 19 +  3] = b04;
    bio[i * 19 +  4] = b05;
    bio[i * 19 +  5] = b06;
    bio[i * 19 +  6] = b07;
    bio[i * 19 +  7] = t_qs[wet_q]  / 3.0;   // BIO08
    bio[i * 19 +  8] = t_qs[dry_q]  / 3.0;   // BIO09
    bio[i * 19 +  9] = t_qs[warm_q] / 3.0;   // BIO10
    bio[i * 19 + 10] = t_qs[cold_q] / 3.0;   // BIO11
    bio[i * 19 + 11] = b12;
    bio[i * 19 + 12] = b13;
    bio[i * 19 + 13] = b14;
    bio[i * 19 + 14] = b15;
    bio[i * 19 + 15] = pr_qs[wet_q];           // BIO16
    bio[i * 19 + 16] = pr_qs[dry_q];           // BIO17
    bio[i * 19 + 17] = pr_qs[warm_q];          // BIO18
    bio[i * 19 + 18] = pr_qs[cold_q];          // BIO19
}

// ── Host launcher ─────────────────────────────────────────────────────────────

void launch_bioclim_cuda(
    const double* h_tas,
    const double* h_tasmax,
    const double* h_tasmin,
    const double* h_pr,
    const double* h_mask,
    double*       h_bio,
    int           n_pix
) {
    const std::size_t var_bytes  = static_cast<std::size_t>(12 * n_pix) * sizeof(double);
    const std::size_t mask_bytes = static_cast<std::size_t>(n_pix)      * sizeof(double);
    const std::size_t bio_bytes  = static_cast<std::size_t>(19 * n_pix) * sizeof(double);

    double *d_tas = nullptr, *d_tasmax = nullptr;
    double *d_tasmin = nullptr, *d_pr = nullptr;
    double *d_bio = nullptr, *d_mask = nullptr;

    // Allocate and transfer input data — clean up everything on any error.
    try {
        CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_tas),    var_bytes));
        CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_tasmax), var_bytes));
        CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_tasmin), var_bytes));
        CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_pr),     var_bytes));
        CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_bio),    bio_bytes));

        CUDA_CHECK(cudaMemcpy(d_tas,    h_tas,    var_bytes, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_tasmax, h_tasmax, var_bytes, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_tasmin, h_tasmin, var_bytes, cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_pr,     h_pr,     var_bytes, cudaMemcpyHostToDevice));

        if (h_mask != nullptr) {
            CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_mask), mask_bytes));
            CUDA_CHECK(cudaMemcpy(d_mask, h_mask, mask_bytes, cudaMemcpyHostToDevice));
        }

        const int block_size = 256;
        const int grid_size  = (n_pix + block_size - 1) / block_size;
        bioclim_kernel<<<grid_size, block_size>>>(
            d_tas, d_tasmax, d_tasmin, d_pr, d_mask, d_bio, n_pix);

        CUDA_CHECK(cudaDeviceSynchronize());
        CUDA_CHECK(cudaMemcpy(h_bio, d_bio, bio_bytes, cudaMemcpyDeviceToHost));
    } catch (...) {
        // Free any successfully allocated device memory before re-throwing.
        if (d_tas)    cudaFree(d_tas);
        if (d_tasmax) cudaFree(d_tasmax);
        if (d_tasmin) cudaFree(d_tasmin);
        if (d_pr)     cudaFree(d_pr);
        if (d_bio)    cudaFree(d_bio);
        if (d_mask)   cudaFree(d_mask);
        throw;
    }

    cudaFree(d_tas);
    cudaFree(d_tasmax);
    cudaFree(d_tasmin);
    cudaFree(d_pr);
    cudaFree(d_bio);
    if (d_mask != nullptr) cudaFree(d_mask);
}

#endif  // HAVE_CUDA
