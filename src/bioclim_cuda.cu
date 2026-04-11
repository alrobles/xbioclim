// bioclim_cuda.cu — Per-pixel CUDA kernel for BIO1–BIO19 computation.
//
// Each CUDA thread processes one pixel. Input layout:
//   tas[p*12 + m], tasmax[p*12 + m], tasmin[p*12 + m], pr[p*12 + m]
// Output layout:
//   bio[p*19 + b]  (b = 0..18 for BIO01..BIO19)

#include "xbioclim/bioclim_cuda.hpp"
#include "xbioclim/primitives.hpp"

#include <cuda_runtime.h>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace xbioclim {

// ---------------------------------------------------------------------------
// Device helpers (inlined into the kernel)
// ---------------------------------------------------------------------------

/// Mean of 12 values starting at base pointer.
__device__ static float dev_mean12(const float* v) {
    float s = 0.0f;
    #pragma unroll
    for (int i = 0; i < 12; ++i) s += v[i];
    return s / 12.0f;
}

/// Sum of 12 values.
__device__ static float dev_sum12(const float* v) {
    float s = 0.0f;
    #pragma unroll
    for (int i = 0; i < 12; ++i) s += v[i];
    return s;
}

/// Maximum of 12 values.
__device__ static float dev_max12(const float* v) {
    float m = v[0];
    #pragma unroll
    for (int i = 1; i < 12; ++i) if (v[i] > m) m = v[i];
    return m;
}

/// Minimum of 12 values.
__device__ static float dev_min12(const float* v) {
    float m = v[0];
    #pragma unroll
    for (int i = 1; i < 12; ++i) if (v[i] < m) m = v[i];
    return m;
}

/// Population standard deviation of 12 values (ddof=0).
__device__ static float dev_std12(const float* v) {
    float mu = dev_mean12(v);
    float var = 0.0f;
    #pragma unroll
    for (int i = 0; i < 12; ++i) {
        float d = v[i] - mu;
        var += d * d;
    }
    return sqrtf(var / 12.0f);
}

/// 0-based start index of the 3-month rolling window with the maximum
/// circular sum.  Returns 12 (sentinel) if any value is NaN.
__device__ static int dev_rolling_argmax(const float* v) {
    float best = -3.4028235e+38f; // FLT_MIN lowest
    int   bi   = 0;
    #pragma unroll
    for (int i = 0; i < 12; ++i) {
        float a = v[i];
        float b = v[(i + 1) % 12];
        float c = v[(i + 2) % 12];
        if (isnan(a) || isnan(b) || isnan(c)) return 12; // sentinel
        float s = a + b + c;
        if (s > best) { best = s; bi = i; }
    }
    return bi;
}

/// 0-based start index of the 3-month rolling window with the minimum
/// circular sum.  Returns 12 (sentinel) if any value is NaN.
__device__ static int dev_rolling_argmin(const float* v) {
    float best = 3.4028235e+38f; // FLT_MAX
    int   bi   = 0;
    #pragma unroll
    for (int i = 0; i < 12; ++i) {
        float a = v[i];
        float b = v[(i + 1) % 12];
        float c = v[(i + 2) % 12];
        if (isnan(a) || isnan(b) || isnan(c)) return 12; // sentinel
        float s = a + b + c;
        if (s < best) { best = s; bi = i; }
    }
    return bi;
}

/// Mean of the circular 3-month window starting at index `start`.
/// Returns NaN when start == 12 (sentinel).
__device__ static float dev_quarter_mean(const float* v, int start) {
    if (start == 12) return __int_as_float(0x7FC00000); // NaN
    return (v[start % 12] + v[(start + 1) % 12] + v[(start + 2) % 12]) / 3.0f;
}

// ---------------------------------------------------------------------------
// Main CUDA kernel — one thread per pixel
// ---------------------------------------------------------------------------

__global__ void bioclim_kernel(
    const float* __restrict__ tas,    // [N, 12]
    const float* __restrict__ tasmax, // [N, 12]
    const float* __restrict__ tasmin, // [N, 12]
    const float* __restrict__ pr,     // [N, 12]
    float* __restrict__       bio,    // [N, 19]
    int N)
{
    const int p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= N) return;

    const float* t    = tas    + p * 12;
    const float* tmax = tasmax + p * 12;
    const float* tmin = tasmin + p * 12;
    const float* rain = pr     + p * 12;
    float*       out  = bio    + p * 19;

    // Diurnal range per month (local registers)
    float diurnal[12];
    #pragma unroll
    for (int m = 0; m < 12; ++m) diurnal[m] = tmax[m] - tmin[m];

    // Simple statistics
    const float bio01 = dev_mean12(t);
    const float bio02 = dev_mean12(diurnal);
    const float bio05 = dev_max12(tmax);
    const float bio06 = dev_min12(tmin);
    const float bio07 = bio05 - bio06;

    const float nan_val = __int_as_float(0x7FC00000);

    // BIO03: guard division by zero
    const float bio03 = (bio07 != 0.0f) ? (100.0f * bio02 / bio07) : nan_val;
    const float bio04 = 100.0f * dev_std12(t);

    const float bio12 = dev_sum12(rain);
    const float bio13 = dev_max12(rain);
    const float bio14 = dev_min12(rain);

    // BIO15: precipitation CV
    const float pr_mean = dev_mean12(rain);
    const float pr_std  = dev_std12(rain);
    const float bio15 = (pr_mean > 0.0f) ? (100.0f * pr_std / pr_mean) : nan_val;

    // Rolling quarter indices
    const int wet_q  = dev_rolling_argmax(rain);
    const int dry_q  = dev_rolling_argmin(rain);
    const int warm_q = dev_rolling_argmax(t);
    const int cold_q = dev_rolling_argmin(t);

    // Quarter means
    const float bio08 = dev_quarter_mean(t,    wet_q);
    const float bio09 = dev_quarter_mean(t,    dry_q);
    const float bio10 = dev_quarter_mean(t,    warm_q);
    const float bio11 = dev_quarter_mean(t,    cold_q);
    const float bio16 = dev_quarter_mean(rain, wet_q);
    const float bio17 = dev_quarter_mean(rain, dry_q);
    const float bio18 = dev_quarter_mean(rain, warm_q);
    const float bio19 = dev_quarter_mean(rain, cold_q);

    // Write results (index 0..18 = BIO01..BIO19)
    out[ 0] = bio01;
    out[ 1] = bio02;
    out[ 2] = bio03;
    out[ 3] = bio04;
    out[ 4] = bio05;
    out[ 5] = bio06;
    out[ 6] = bio07;
    out[ 7] = bio08;
    out[ 8] = bio09;
    out[ 9] = bio10;
    out[10] = bio11;
    out[11] = bio12;
    out[12] = bio13;
    out[13] = bio14;
    out[14] = bio15;
    out[15] = bio16;
    out[16] = bio17;
    out[17] = bio18;
    out[18] = bio19;
}

// ---------------------------------------------------------------------------
// Host wrapper
// ---------------------------------------------------------------------------

BioBlock compute_bioclim_cuda(const ClimateBlock& data) {
    // --- Validate ---
    const std::size_t N = data.tas.shape(0);
    if (data.tasmax.shape(0) != N ||
        data.tasmin.shape(0) != N ||
        data.pr.shape(0)     != N) {
        throw std::invalid_argument(
            "ClimateBlock: all fields must have the same n_pixels (shape(0))");
    }
    if (data.tas.shape(1)    != 12 ||
        data.tasmax.shape(1) != 12 ||
        data.tasmin.shape(1) != 12 ||
        data.pr.shape(1)     != 12) {
        throw std::invalid_argument(
            "ClimateBlock: all fields must have exactly 12 months (shape(1))");
    }

    const int n = static_cast<int>(N);
    const std::size_t input_bytes  = N * 12 * sizeof(float);
    const std::size_t output_bytes = N * 19 * sizeof(float);

    // --- Allocate device memory ---
    float *d_tas = nullptr, *d_tasmax = nullptr,
          *d_tasmin = nullptr, *d_pr = nullptr, *d_bio = nullptr;

    auto cuda_check = [](cudaError_t e, const char* msg) {
        if (e != cudaSuccess) {
            throw std::runtime_error(
                std::string(msg) + ": " + cudaGetErrorString(e));
        }
    };

    cuda_check(cudaMalloc(&d_tas,    input_bytes),  "cudaMalloc tas");
    cuda_check(cudaMalloc(&d_tasmax, input_bytes),  "cudaMalloc tasmax");
    cuda_check(cudaMalloc(&d_tasmin, input_bytes),  "cudaMalloc tasmin");
    cuda_check(cudaMalloc(&d_pr,     input_bytes),  "cudaMalloc pr");
    cuda_check(cudaMalloc(&d_bio,    output_bytes), "cudaMalloc bio");

    // --- Copy inputs to device (xtensor data is row-major contiguous) ---
    cuda_check(cudaMemcpy(d_tas,    data.tas.data(),    input_bytes, cudaMemcpyHostToDevice), "cudaMemcpy tas");
    cuda_check(cudaMemcpy(d_tasmax, data.tasmax.data(), input_bytes, cudaMemcpyHostToDevice), "cudaMemcpy tasmax");
    cuda_check(cudaMemcpy(d_tasmin, data.tasmin.data(), input_bytes, cudaMemcpyHostToDevice), "cudaMemcpy tasmin");
    cuda_check(cudaMemcpy(d_pr,     data.pr.data(),     input_bytes, cudaMemcpyHostToDevice), "cudaMemcpy pr");

    // --- Launch kernel ---
    const int threads_per_block = 256;
    const int blocks = (n + threads_per_block - 1) / threads_per_block;
    bioclim_kernel<<<blocks, threads_per_block>>>(
        d_tas, d_tasmax, d_tasmin, d_pr, d_bio, n);
    cuda_check(cudaGetLastError(), "bioclim_kernel launch");
    cuda_check(cudaDeviceSynchronize(), "cudaDeviceSynchronize");

    // --- Copy results to host ---
    std::vector<float> h_bio(N * 19);
    cuda_check(cudaMemcpy(h_bio.data(), d_bio, output_bytes, cudaMemcpyDeviceToHost), "cudaMemcpy bio back");

    // --- Free device memory ---
    cudaFree(d_tas);
    cudaFree(d_tasmax);
    cudaFree(d_tasmin);
    cudaFree(d_pr);
    cudaFree(d_bio);

    // --- Pack results into BioBlock ---
    BioBlock bio;
    bio.bio01 = Array1D::from_shape({N});
    bio.bio02 = Array1D::from_shape({N});
    bio.bio03 = Array1D::from_shape({N});
    bio.bio04 = Array1D::from_shape({N});
    bio.bio05 = Array1D::from_shape({N});
    bio.bio06 = Array1D::from_shape({N});
    bio.bio07 = Array1D::from_shape({N});
    bio.bio08 = Array1D::from_shape({N});
    bio.bio09 = Array1D::from_shape({N});
    bio.bio10 = Array1D::from_shape({N});
    bio.bio11 = Array1D::from_shape({N});
    bio.bio12 = Array1D::from_shape({N});
    bio.bio13 = Array1D::from_shape({N});
    bio.bio14 = Array1D::from_shape({N});
    bio.bio15 = Array1D::from_shape({N});
    bio.bio16 = Array1D::from_shape({N});
    bio.bio17 = Array1D::from_shape({N});
    bio.bio18 = Array1D::from_shape({N});
    bio.bio19 = Array1D::from_shape({N});

    for (std::size_t p = 0; p < N; ++p) {
        const float* out = h_bio.data() + p * 19;
        bio.bio01(p) = out[ 0];
        bio.bio02(p) = out[ 1];
        bio.bio03(p) = out[ 2];
        bio.bio04(p) = out[ 3];
        bio.bio05(p) = out[ 4];
        bio.bio06(p) = out[ 5];
        bio.bio07(p) = out[ 6];
        bio.bio08(p) = out[ 7];
        bio.bio09(p) = out[ 8];
        bio.bio10(p) = out[ 9];
        bio.bio11(p) = out[10];
        bio.bio12(p) = out[11];
        bio.bio13(p) = out[12];
        bio.bio14(p) = out[13];
        bio.bio15(p) = out[14];
        bio.bio16(p) = out[15];
        bio.bio17(p) = out[16];
        bio.bio18(p) = out[17];
        bio.bio19(p) = out[18];
    }

    return bio;
}

} // namespace xbioclim
