// era5_monthly.h — ERA5-Land hourly-to-monthly aggregation for xbioclim.
//
// Converts ERA5-Land hourly reanalysis data to CHELSA-compatible monthly
// climate variables suitable for bioclimatic variable computation.
//
// ERA5-Land input variables
// -------------------------
//   t2m — 2-metre temperature (K), hourly instantaneous
//   tp  — total precipitation (m), hourly accumulation
//
// CHELSA-compatible output variables
// ----------------------------------
//   tas    — monthly mean of daily mean temperature (K or °C)
//   tasmax — monthly mean of daily maximum temperature (K or °C)
//   tasmin — monthly mean of daily minimum temperature (K or °C)
//   pr     — monthly total precipitation (kg m-2 month-1 = mm)
//
// Aggregation method
// ------------------
//   1. Group hourly t2m by day (24 hours per day).
//   2. Per day: mean → daily mean; max → daily max; min → daily min.
//   3. Per month: mean of daily means → tas;
//                 mean of daily maxima → tasmax;
//                 mean of daily minima → tasmin.
//   4. Sum hourly tp over the month, multiply by 1000 (m → mm) → pr.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace xbioclim {
namespace era5 {

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Kelvin-to-Celsius offset.
static constexpr double K_TO_C = -273.15;

/// Metres of water to mm (= kg m-2).
static constexpr double M_TO_MM = 1000.0;

// ---------------------------------------------------------------------------
// Single-pixel aggregation (operates on raw double arrays)
// ---------------------------------------------------------------------------

/// Aggregate hourly temperature for one pixel to monthly statistics.
///
/// @param hourly  Pointer to n_hours hourly t2m values (K).
/// @param n_hours Number of hourly values (must equal 24 × n_days).
/// @param n_days  Number of days in the month.
/// @param[out] tas    Monthly mean of daily means.
/// @param[out] tasmax Monthly mean of daily maxima.
/// @param[out] tasmin Monthly mean of daily minima.
/// @param to_celsius  If true, subtract 273.15 to convert K → °C.
inline void t2m_to_monthly(const double* hourly,
                           int n_hours,
                           int n_days,
                           double& tas,
                           double& tasmax,
                           double& tasmin,
                           bool to_celsius = false) {
    const double offset = to_celsius ? K_TO_C : 0.0;
    double sum_dmean = 0.0;
    double sum_dmax  = 0.0;
    double sum_dmin  = 0.0;

    for (int d = 0; d < n_days; ++d) {
        const int base = d * 24;
        double dsum = 0.0;
        double dmax = hourly[base];
        double dmin = hourly[base];
        int count = 0;

        for (int h = 0; h < 24 && (base + h) < n_hours; ++h) {
            double v = hourly[base + h];
            if (std::isnan(v)) continue;
            dsum += v;
            if (v > dmax) dmax = v;
            if (v < dmin) dmin = v;
            ++count;
        }

        if (count > 0) {
            sum_dmean += dsum / static_cast<double>(count);
            sum_dmax  += dmax;
            sum_dmin  += dmin;
        } else {
            sum_dmean += std::numeric_limits<double>::quiet_NaN();
            sum_dmax  += std::numeric_limits<double>::quiet_NaN();
            sum_dmin  += std::numeric_limits<double>::quiet_NaN();
        }
    }

    tas    = (sum_dmean / static_cast<double>(n_days)) + offset;
    tasmax = (sum_dmax  / static_cast<double>(n_days)) + offset;
    tasmin = (sum_dmin  / static_cast<double>(n_days)) + offset;
}

/// Aggregate hourly precipitation for one pixel to monthly total.
///
/// @param hourly  Pointer to n_hours hourly tp values (m).
/// @param n_hours Number of hourly values.
/// @param[out] pr Monthly total precipitation (kg m-2 month-1 = mm).
inline void tp_to_monthly(const double* hourly,
                          int n_hours,
                          double& pr) {
    double total = 0.0;
    for (int h = 0; h < n_hours; ++h) {
        double v = hourly[h];
        if (!std::isnan(v) && v > 0.0) {
            total += v;
        }
    }
    pr = total * M_TO_MM;
}

// ---------------------------------------------------------------------------
// Multi-pixel batch aggregation
// ---------------------------------------------------------------------------

/// Aggregate hourly t2m for n_pixels to monthly statistics.
///
/// @param hourly  Column-major matrix [n_pixels × n_hours].
///                Element [i, h] is at hourly[i + h * n_pixels].
/// @param n_pixels Number of grid cells.
/// @param n_hours  Total hourly steps (must equal 24 × n_days).
/// @param n_days   Number of days in the month.
/// @param[out] tas    Array of length n_pixels.
/// @param[out] tasmax Array of length n_pixels.
/// @param[out] tasmin Array of length n_pixels.
/// @param to_celsius If true, convert K → °C.
inline void t2m_to_monthly_batch(const double* hourly,
                                 int n_pixels,
                                 int n_hours,
                                 int n_days,
                                 double* tas,
                                 double* tasmax,
                                 double* tasmin,
                                 bool to_celsius = false) {
    const double offset = to_celsius ? K_TO_C : 0.0;

    for (int i = 0; i < n_pixels; ++i) {
        double sum_dmean = 0.0;
        double sum_dmax  = 0.0;
        double sum_dmin  = 0.0;

        for (int d = 0; d < n_days; ++d) {
            const int h_start = d * 24;
            double dsum = 0.0;
            double dmax = hourly[i + h_start * n_pixels];
            double dmin = dmax;
            int count = 0;

            for (int h = 0; h < 24 && (h_start + h) < n_hours; ++h) {
                double v = hourly[i + (h_start + h) * n_pixels];
                if (std::isnan(v)) continue;
                dsum += v;
                if (v > dmax) dmax = v;
                if (v < dmin) dmin = v;
                ++count;
            }

            if (count > 0) {
                sum_dmean += dsum / static_cast<double>(count);
                sum_dmax  += dmax;
                sum_dmin  += dmin;
            }
        }

        tas[i]    = (sum_dmean / static_cast<double>(n_days)) + offset;
        tasmax[i] = (sum_dmax  / static_cast<double>(n_days)) + offset;
        tasmin[i] = (sum_dmin  / static_cast<double>(n_days)) + offset;
    }
}

/// Aggregate hourly tp for n_pixels to monthly totals.
///
/// @param hourly  Column-major matrix [n_pixels × n_hours].
/// @param n_pixels Number of grid cells.
/// @param n_hours  Total hourly steps.
/// @param[out] pr Array of length n_pixels (mm / month).
inline void tp_to_monthly_batch(const double* hourly,
                                int n_pixels,
                                int n_hours,
                                double* pr) {
    for (int i = 0; i < n_pixels; ++i) {
        double total = 0.0;
        for (int h = 0; h < n_hours; ++h) {
            double v = hourly[i + h * n_pixels];
            if (!std::isnan(v) && v > 0.0) {
                total += v;
            }
        }
        pr[i] = total * M_TO_MM;
    }
}

}  // namespace era5
}  // namespace xbioclim
