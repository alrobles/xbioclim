// xbioclim.h — Self-contained C++ implementation of the xbioclim core.
//
// This header mirrors the API of the xbioclim C++ library
// (github.com/alrobles/xbioclim) but uses only the C++ standard library,
// making it embeddable in an R package without extra dependencies.
//
// Key types
// ---------
//   ClimateBlock  — input bundle: four arrays of shape [n_pixels × 12]
//   BioBlock      — output bundle: 19 arrays of shape [n_pixels]
//
// Key functions
// -------------
//   compute_pixel(tas, tasmax, tasmin, pr)  — single-pixel path
//   compute_bioclim(ClimateBlock)           — multi-pixel batch path

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

namespace xbioclim {

// ---------------------------------------------------------------------------
// Primitive helpers (operate on raw double arrays of length 12)
// ---------------------------------------------------------------------------

/// Population standard deviation (denominator N, not N-1).
inline double sd_pop(const double* x, std::size_t n) {
    if (n == 0) return std::numeric_limits<double>::quiet_NaN();
    double m = 0.0;
    for (std::size_t i = 0; i < n; ++i) m += x[i];
    m /= static_cast<double>(n);
    double sq = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double d = x[i] - m;
        sq += d * d;
    }
    return std::sqrt(sq / static_cast<double>(n));
}

/// 0-based index of the rolling 3-month window with maximum circular sum.
inline std::size_t rolling_quarter_argmax(const double* x, std::size_t n = 12) {
    std::size_t best_i = 0;
    double best = std::numeric_limits<double>::lowest();
    for (std::size_t i = 0; i < n; ++i) {
        double s = x[i] + x[(i + 1) % n] + x[(i + 2) % n];
        if (s > best) {
            best   = s;
            best_i = i;
        }
    }
    return best_i;
}

/// 0-based index of the rolling 3-month window with minimum circular sum.
inline std::size_t rolling_quarter_argmin(const double* x, std::size_t n = 12) {
    std::size_t best_i = 0;
    double best = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < n; ++i) {
        double s = x[i] + x[(i + 1) % n] + x[(i + 2) % n];
        if (s < best) {
            best   = s;
            best_i = i;
        }
    }
    return best_i;
}

/// Mean of 3 consecutive months (circular wrap) starting at 0-based index.
inline double quarter_mean(const double* x, std::size_t start,
                           std::size_t n = 12) {
    return (x[start % n] + x[(start + 1) % n] + x[(start + 2) % n]) / 3.0;
}

/// Sum of 3 consecutive months (circular wrap) starting at 0-based index.
inline double quarter_sum(const double* x, std::size_t start,
                          std::size_t n = 12) {
    return x[start % n] + x[(start + 1) % n] + x[(start + 2) % n];
}

// ---------------------------------------------------------------------------
// Single-pixel computation
// ---------------------------------------------------------------------------

/// Compute all 19 bioclimatic variables for a single pixel.
///
/// @param tas    Pointer to 12 monthly mean temperatures.
/// @param tasmax Pointer to 12 monthly maximum temperatures.
/// @param tasmin Pointer to 12 monthly minimum temperatures.
/// @param pr     Pointer to 12 monthly precipitation values.
/// @return std::array<double, 19> ordered bio01 … bio19.
inline std::array<double, 19> compute_pixel(const double* tas,
                                            const double* tasmax,
                                            const double* tasmin,
                                            const double* pr) {
    // BIO01 — Mean Annual Temperature
    double b01 = 0.0;
    for (int m = 0; m < 12; ++m) b01 += tas[m];
    b01 /= 12.0;

    // BIO02 — Mean Diurnal Range
    double b02 = 0.0;
    for (int m = 0; m < 12; ++m) b02 += (tasmax[m] - tasmin[m]);
    b02 /= 12.0;

    // BIO05 / BIO06 / BIO07 — Max/Min temp and Annual Range
    double b05 = *std::max_element(tasmax, tasmax + 12);
    double b06 = *std::min_element(tasmin, tasmin + 12);
    double b07 = b05 - b06;

    // BIO03 — Isothermality (guard against zero range)
    double b03 = (b07 == 0.0)
        ? std::numeric_limits<double>::quiet_NaN()
        : 100.0 * b02 / b07;

    // BIO04 — Temperature Seasonality
    double b04 = 100.0 * sd_pop(tas, 12);

    // Rolling quarter indices (0-based)
    std::size_t wet_q  = rolling_quarter_argmax(pr);
    std::size_t dry_q  = rolling_quarter_argmin(pr);
    std::size_t warm_q = rolling_quarter_argmax(tas);
    std::size_t cold_q = rolling_quarter_argmin(tas);

    // BIO08 / BIO09 — Mean temp of wettest/driest quarter
    double b08 = quarter_mean(tas, wet_q);
    double b09 = quarter_mean(tas, dry_q);

    // BIO10 / BIO11 — Mean temp of warmest/coldest quarter
    double b10 = quarter_mean(tas, warm_q);
    double b11 = quarter_mean(tas, cold_q);

    // BIO12 — Annual Precipitation
    double b12 = 0.0;
    for (int m = 0; m < 12; ++m) b12 += pr[m];

    // BIO13 / BIO14 — Wettest/Driest month precipitation
    double b13 = *std::max_element(pr, pr + 12);
    double b14 = *std::min_element(pr, pr + 12);

    // BIO15 — Precipitation Seasonality (CV; guard against zero mean)
    double pr_mean = b12 / 12.0;
    double b15 = (pr_mean == 0.0)
        ? std::numeric_limits<double>::quiet_NaN()
        : 100.0 * sd_pop(pr, 12) / pr_mean;

    // BIO16 / BIO17 — Precipitation of Wettest/Driest Quarter (total sum)
    double b16 = quarter_sum(pr, wet_q);
    double b17 = quarter_sum(pr, dry_q);

    // BIO18 / BIO19 — Precipitation of Warmest/Coldest Quarter (total sum)
    double b18 = quarter_sum(pr, warm_q);
    double b19 = quarter_sum(pr, cold_q);

    return {b01, b02, b03, b04, b05, b06, b07,
            b08, b09, b10, b11,
            b12, b13, b14, b15,
            b16, b17, b18, b19};
}

// ---------------------------------------------------------------------------
// Multi-pixel data structures (mirrors xbioclim ClimateBlock / BioBlock)
// ---------------------------------------------------------------------------

/// Input bundle: four climate variable arrays stored row-major
/// (pixel index × 12 months).  data[pixel * 12 + month].
struct ClimateBlock {
    std::vector<double> tas;    ///< Monthly mean temperature  [n_pixels × 12]
    std::vector<double> tasmax; ///< Monthly max temperature   [n_pixels × 12]
    std::vector<double> tasmin; ///< Monthly min temperature   [n_pixels × 12]
    std::vector<double> pr;     ///< Monthly precipitation     [n_pixels × 12]
    std::size_t n_pixels;       ///< Number of pixels

    ClimateBlock() : n_pixels(0) {}

    ClimateBlock(std::vector<double> tas_in,
                 std::vector<double> tasmax_in,
                 std::vector<double> tasmin_in,
                 std::vector<double> pr_in,
                 std::size_t npix)
        : tas(std::move(tas_in))
        , tasmax(std::move(tasmax_in))
        , tasmin(std::move(tasmin_in))
        , pr(std::move(pr_in))
        , n_pixels(npix) {}
};

/// Output bundle: one value per pixel for each of the 19 bio variables.
struct BioBlock {
    std::vector<double> bio01, bio02, bio03, bio04, bio05,
                        bio06, bio07, bio08, bio09, bio10,
                        bio11, bio12, bio13, bio14, bio15,
                        bio16, bio17, bio18, bio19;
};

// ---------------------------------------------------------------------------
// Multi-pixel batch computation
// ---------------------------------------------------------------------------

/// Compute all 19 bioclimatic variables for every pixel in a ClimateBlock.
inline BioBlock compute_bioclim(const ClimateBlock& data) {
    std::size_t N = data.n_pixels;

    if (data.tas.size()    != N * 12 ||
        data.tasmax.size() != N * 12 ||
        data.tasmin.size() != N * 12 ||
        data.pr.size()     != N * 12) {
        throw std::invalid_argument(
            "ClimateBlock: each vector must have n_pixels * 12 elements");
    }

    BioBlock bio;
    bio.bio01.resize(N); bio.bio02.resize(N); bio.bio03.resize(N);
    bio.bio04.resize(N); bio.bio05.resize(N); bio.bio06.resize(N);
    bio.bio07.resize(N); bio.bio08.resize(N); bio.bio09.resize(N);
    bio.bio10.resize(N); bio.bio11.resize(N); bio.bio12.resize(N);
    bio.bio13.resize(N); bio.bio14.resize(N); bio.bio15.resize(N);
    bio.bio16.resize(N); bio.bio17.resize(N); bio.bio18.resize(N);
    bio.bio19.resize(N);

    for (std::size_t p = 0; p < N; ++p) {
        const double* tas_p    = data.tas.data()    + p * 12;
        const double* tasmax_p = data.tasmax.data() + p * 12;
        const double* tasmin_p = data.tasmin.data() + p * 12;
        const double* pr_p     = data.pr.data()     + p * 12;

        auto res = compute_pixel(tas_p, tasmax_p, tasmin_p, pr_p);

        bio.bio01[p] = res[ 0]; bio.bio02[p] = res[ 1];
        bio.bio03[p] = res[ 2]; bio.bio04[p] = res[ 3];
        bio.bio05[p] = res[ 4]; bio.bio06[p] = res[ 5];
        bio.bio07[p] = res[ 6]; bio.bio08[p] = res[ 7];
        bio.bio09[p] = res[ 8]; bio.bio10[p] = res[ 9];
        bio.bio11[p] = res[10]; bio.bio12[p] = res[11];
        bio.bio13[p] = res[12]; bio.bio14[p] = res[13];
        bio.bio15[p] = res[14]; bio.bio16[p] = res[15];
        bio.bio17[p] = res[16]; bio.bio18[p] = res[17];
        bio.bio19[p] = res[18];
    }
    return bio;
}

} // namespace xbioclim
