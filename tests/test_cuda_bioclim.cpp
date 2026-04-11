// test_cuda_bioclim.cpp — Tests comparing CPU and CUDA backend output.
//
// These tests are compiled and run only when XBIOCLIM_USE_CUDA is defined.
// When no CUDA-capable GPU is available at runtime the tests are skipped
// gracefully using Catch2's SKIP() macro.

#ifdef XBIOCLIM_USE_CUDA

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "xbioclim/bioclim.hpp"
#include "xbioclim/bioclim_cuda.hpp"

#include <cmath>
#include <cuda_runtime.h>
#include <limits>
#include <stdexcept>

using Catch::Matchers::WithinAbs;
using namespace xbioclim;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Returns true when at least one CUDA device is available.
static bool cuda_device_available() {
    int count = 0;
    return cudaGetDeviceCount(&count) == cudaSuccess && count > 0;
}

/// Synthetic test block (same as test_bioclim.cpp).
static ClimateBlock make_test_block(std::size_t N = 1) {
    ClimateBlock data;
    data.tas    = Array2D::from_shape({N, 12});
    data.tasmax = Array2D::from_shape({N, 12});
    data.tasmin = Array2D::from_shape({N, 12});
    data.pr     = Array2D::from_shape({N, 12});

    for (std::size_t p = 0; p < N; ++p) {
        for (std::size_t m = 0; m < 12; ++m) {
            data.tas   (p, m) = static_cast<float>(m + 1);
            data.tasmax(p, m) = static_cast<float>(m + 2);
            data.tasmin(p, m) = static_cast<float>(m);
            data.pr    (p, m) = static_cast<float>(m + 1);
        }
    }
    return data;
}

/// Tolerance for CPU vs CUDA comparison (single-precision rounding).
static constexpr float kTol = 1e-3f;

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

TEST_CASE("CUDA: output matches CPU for single pixel", "[cuda][bioclim]") {
    if (!cuda_device_available()) {
        SKIP("No CUDA device available");
    }

    auto data    = make_test_block(1);
    auto cpu_bio = compute_bioclim(data);
    auto gpu_bio = compute_bioclim_cuda(data);

    // BIO01–BIO19 must match the CPU reference to within kTol
    CHECK_THAT(gpu_bio.bio01(0), WithinAbs(cpu_bio.bio01(0), kTol));
    CHECK_THAT(gpu_bio.bio02(0), WithinAbs(cpu_bio.bio02(0), kTol));
    CHECK_THAT(gpu_bio.bio03(0), WithinAbs(cpu_bio.bio03(0), kTol));
    CHECK_THAT(gpu_bio.bio04(0), WithinAbs(cpu_bio.bio04(0), kTol));
    CHECK_THAT(gpu_bio.bio05(0), WithinAbs(cpu_bio.bio05(0), kTol));
    CHECK_THAT(gpu_bio.bio06(0), WithinAbs(cpu_bio.bio06(0), kTol));
    CHECK_THAT(gpu_bio.bio07(0), WithinAbs(cpu_bio.bio07(0), kTol));
    CHECK_THAT(gpu_bio.bio08(0), WithinAbs(cpu_bio.bio08(0), kTol));
    CHECK_THAT(gpu_bio.bio09(0), WithinAbs(cpu_bio.bio09(0), kTol));
    CHECK_THAT(gpu_bio.bio10(0), WithinAbs(cpu_bio.bio10(0), kTol));
    CHECK_THAT(gpu_bio.bio11(0), WithinAbs(cpu_bio.bio11(0), kTol));
    CHECK_THAT(gpu_bio.bio12(0), WithinAbs(cpu_bio.bio12(0), kTol));
    CHECK_THAT(gpu_bio.bio13(0), WithinAbs(cpu_bio.bio13(0), kTol));
    CHECK_THAT(gpu_bio.bio14(0), WithinAbs(cpu_bio.bio14(0), kTol));
    CHECK_THAT(gpu_bio.bio15(0), WithinAbs(cpu_bio.bio15(0), kTol));
    CHECK_THAT(gpu_bio.bio16(0), WithinAbs(cpu_bio.bio16(0), kTol));
    CHECK_THAT(gpu_bio.bio17(0), WithinAbs(cpu_bio.bio17(0), kTol));
    CHECK_THAT(gpu_bio.bio18(0), WithinAbs(cpu_bio.bio18(0), kTol));
    CHECK_THAT(gpu_bio.bio19(0), WithinAbs(cpu_bio.bio19(0), kTol));
}

TEST_CASE("CUDA: output matches CPU for multiple pixels", "[cuda][bioclim]") {
    if (!cuda_device_available()) {
        SKIP("No CUDA device available");
    }

    const std::size_t N = 1024;
    auto data    = make_test_block(N);
    auto cpu_bio = compute_bioclim(data);
    auto gpu_bio = compute_bioclim_cuda(data);

    for (std::size_t p = 0; p < N; ++p) {
        CHECK_THAT(gpu_bio.bio01(p), WithinAbs(cpu_bio.bio01(p), kTol));
        CHECK_THAT(gpu_bio.bio04(p), WithinAbs(cpu_bio.bio04(p), kTol));
        CHECK_THAT(gpu_bio.bio12(p), WithinAbs(cpu_bio.bio12(p), kTol));
        CHECK_THAT(gpu_bio.bio15(p), WithinAbs(cpu_bio.bio15(p), kTol));
    }
}

TEST_CASE("CUDA: NaN inputs propagate correctly", "[cuda][bioclim]") {
    if (!cuda_device_available()) {
        SKIP("No CUDA device available");
    }

    const float nan = std::numeric_limits<float>::quiet_NaN();
    ClimateBlock data;
    data.tas    = Array2D::from_shape({std::size_t(1), std::size_t(12)});
    data.tasmax = Array2D::from_shape({std::size_t(1), std::size_t(12)});
    data.tasmin = Array2D::from_shape({std::size_t(1), std::size_t(12)});
    data.pr     = Array2D::from_shape({std::size_t(1), std::size_t(12)});
    data.tas.fill(nan);
    data.tasmax.fill(nan);
    data.tasmin.fill(nan);
    data.pr.fill(nan);

    auto gpu_bio = compute_bioclim_cuda(data);

    CHECK(std::isnan(gpu_bio.bio01(0)));
    CHECK(std::isnan(gpu_bio.bio02(0)));
    CHECK(std::isnan(gpu_bio.bio04(0)));
    CHECK(std::isnan(gpu_bio.bio05(0)));
    CHECK(std::isnan(gpu_bio.bio06(0)));
    CHECK(std::isnan(gpu_bio.bio07(0)));
    CHECK(std::isnan(gpu_bio.bio08(0)));
    CHECK(std::isnan(gpu_bio.bio09(0)));
    CHECK(std::isnan(gpu_bio.bio10(0)));
    CHECK(std::isnan(gpu_bio.bio11(0)));
    CHECK(std::isnan(gpu_bio.bio12(0)));
    CHECK(std::isnan(gpu_bio.bio13(0)));
    CHECK(std::isnan(gpu_bio.bio14(0)));
    CHECK(std::isnan(gpu_bio.bio16(0)));
    CHECK(std::isnan(gpu_bio.bio17(0)));
    CHECK(std::isnan(gpu_bio.bio18(0)));
    CHECK(std::isnan(gpu_bio.bio19(0)));
}

TEST_CASE("CUDA: BIO03 is NoData when BIO07 is zero", "[cuda][bioclim]") {
    if (!cuda_device_available()) {
        SKIP("No CUDA device available");
    }

    ClimateBlock data;
    data.tas    = xt::ones<float>({std::size_t(1), std::size_t(12)});
    data.tasmax = xt::ones<float>({std::size_t(1), std::size_t(12)});
    data.tasmin = xt::ones<float>({std::size_t(1), std::size_t(12)});
    data.pr     = xt::ones<float>({std::size_t(1), std::size_t(12)});

    auto gpu_bio = compute_bioclim_cuda(data);
    CHECK(std::isnan(gpu_bio.bio03(0)));
}

TEST_CASE("CUDA: BIO15 is NoData when mean precipitation is zero", "[cuda][bioclim]") {
    if (!cuda_device_available()) {
        SKIP("No CUDA device available");
    }

    ClimateBlock data;
    data.tas    = xt::ones<float>({std::size_t(1), std::size_t(12)});
    data.tasmax = xt::ones<float>({std::size_t(1), std::size_t(12)});
    data.tasmin = xt::zeros<float>({std::size_t(1), std::size_t(12)});
    data.pr     = xt::zeros<float>({std::size_t(1), std::size_t(12)});

    auto gpu_bio = compute_bioclim_cuda(data);
    CHECK(std::isnan(gpu_bio.bio15(0)));
}

TEST_CASE("CUDA: rejects mismatched n_pixels", "[cuda][bioclim]") {
    if (!cuda_device_available()) {
        SKIP("No CUDA device available");
    }

    ClimateBlock data;
    data.tas    = Array2D::from_shape({std::size_t(2), std::size_t(12)});
    data.tasmax = Array2D::from_shape({std::size_t(3), std::size_t(12)});
    data.tasmin = Array2D::from_shape({std::size_t(2), std::size_t(12)});
    data.pr     = Array2D::from_shape({std::size_t(2), std::size_t(12)});
    CHECK_THROWS_AS(compute_bioclim_cuda(data), std::invalid_argument);
}

#endif // XBIOCLIM_USE_CUDA
