// bioclim_xt.cpp — vectorized whole-array bioclim computation with a
// zero-copy R↔C++ bridge.
//
// This file implements `bioclim_xt()`, an alternative to `bioclim_cpp()` that:
//   1. Maps R NumericMatrix memory directly onto raw double* arrays (zero-copy
//      on input via REAL() pointer aliasing).
//   2. Pre-allocates the output R NumericMatrix and writes into its REAL()
//      memory directly — zero-copy on output.
//   3. Computes all 19 bioclimatic variables using whole-array, column-oriented
//      loops (batch processing) with OpenMP parallelism.
//
// The computation core is identical to bioclim_cpp() but restructured so that
// per-pixel work is tightly packed: temporary arrays are stack-allocated
// (no heap), and all inner loops operate on contiguous memory.
//
// Part A (Issue #20): The key primitives are:
//   - rolling_quarter_sum_xt   — rolling 3-month sums for all 12 months
//   - rolling_quarter_argmax   — which quarter has max precipitation
//   - rolling_quarter_argmin   — which quarter has min precipitation
//   - column mean/max/min/sum  over the 12-month axis
//
// Part B (Issue #21): Zero-copy bridge — REAL(x) is used directly; no
// intermediate std::vector or matrix copies are made.

#ifdef _OPENMP
#include <omp.h>
#endif
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include "rxbioclim_omp.h"
using namespace Rcpp;

// ── Rolling-quarter primitives (stack-based, OpenMP-safe) ────────────────────

// Isothermality (BIO03) when annual temperature range is zero
static constexpr double BIO03_ZERO_RANGE_FALLBACK = 0.0;

// Compute rolling 3-month sums for all 12 starting months (circular).
// qs[k] = x[k] + x[(k+1)%12] + x[(k+2)%12]
static inline void rolling_quarter_sum_xt(const double* x, double* qs) {
  for (int k = 0; k < 12; k++) {
    qs[k] = x[k] + x[(k + 1) % 12] + x[(k + 2) % 12];
  }
}

// Index of maximum in a 12-element array
static inline int argmax12(const double* x) {
  int best = 0;
  for (int k = 1; k < 12; k++) {
    if (x[k] > x[best]) best = k;
  }
  return best;
}

// Index of minimum in a 12-element array
static inline int argmin12(const double* x) {
  int best = 0;
  for (int k = 1; k < 12; k++) {
    if (x[k] < x[best]) best = k;
  }
  return best;
}

// Sum of a 3-month quarter starting at 0-based index (circular)
static inline double quarter_sum_xt(const double* x, int start) {
  return x[start] + x[(start + 1) % 12] + x[(start + 2) % 12];
}

// Mean of a 3-month quarter starting at 0-based index (circular)
static inline double quarter_mean_xt(const double* x, int start) {
  return quarter_sum_xt(x, start) / 3.0;
}

// Population standard deviation (denominator N) for a 12-element array
static inline double sd_pop_xt(const double* x) {
  double s = 0.0, ss = 0.0;
  for (int i = 0; i < 12; i++) {
    s  += x[i];
    ss += x[i] * x[i];
  }
  double m = s / 12.0;
  return std::sqrt(ss / 12.0 - m * m);
}

// ── bioclim_xt: vectorized whole-array computation with zero-copy bridge ──────

//' Compute all 19 bioclimatic variables (vectorized, zero-copy bridge)
//'
//' A faster alternative to \code{bioclim_cpp()} that uses whole-array
//' vectorized operations and maps R matrix memory directly onto C++ pointers
//' (zero-copy on both input and output).
//'
//' @param tas    Numeric matrix (n_pixels x 12): monthly mean temperature.
//' @param tasmax Numeric matrix (n_pixels x 12): monthly max temperature.
//' @param tasmin Numeric matrix (n_pixels x 12): monthly min temperature.
//' @param pr     Numeric matrix (n_pixels x 12): monthly precipitation.
//' @param ncores Integer: number of OpenMP threads (default 1).
//' @return Numeric matrix (n_pixels x 19) with one column per variable
//'   (bio01..bio19), named accordingly.  Rows with any NA input are returned
//'   as all-NA.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix bioclim_xt(NumericMatrix tas,
                         NumericMatrix tasmax,
                         NumericMatrix tasmin,
                         NumericMatrix pr,
                         int ncores = 1) {
  const int n = tas.nrow();

  // Zero-copy output: pre-allocate R matrix and take its REAL() pointer.
  NumericMatrix result(n, 19);
  CharacterVector cnames = CharacterVector::create(
    "bio01", "bio02", "bio03", "bio04", "bio05",
    "bio06", "bio07", "bio08", "bio09", "bio10",
    "bio11", "bio12", "bio13", "bio14", "bio15",
    "bio16", "bio17", "bio18", "bio19"
  );
  colnames(result) = cnames;

  // Zero-copy input: alias R matrix memory via REAL().
  // Column-major layout: element [i, m] is at ptr[i + m * n].
  const double* tas_p    = REAL(tas);
  const double* tasmax_p = REAL(tasmax);
  const double* tasmin_p = REAL(tasmin);
  const double* pr_p     = REAL(pr);
  double*       res_p    = REAL(result);

#ifdef _OPENMP
  int prev_threads = omp_get_max_threads();
  omp_set_num_threads(safe_omp_threads(ncores));
#endif

  // ── Whole-array parallel loop over pixels ──────────────────────────────────
  //
  // Each iteration processes one pixel.  All intermediate storage is on the
  // stack (private per thread under OpenMP), so no heap allocation occurs
  // inside the loop.  The output is written directly into the R matrix memory.
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < n; i++) {
    // ── Part B (zero-copy): gather column i from each input matrix ───────────
    // Stack buffers — private per thread, no heap allocation.
    double t[12], tmx[12], tmn[12], p[12];
    for (int m = 0; m < 12; m++) {
      t[m]   = tas_p[i   + m * n];
      tmx[m] = tasmax_p[i + m * n];
      tmn[m] = tasmin_p[i + m * n];
      p[m]   = pr_p[i    + m * n];
    }

    // NA guard — write all-NA row and skip if any input is NA.
    bool has_na = false;
    for (int m = 0; m < 12; m++) {
      if (ISNA(t[m]) || ISNA(tmx[m]) || ISNA(tmn[m]) || ISNA(p[m])) {
        has_na = true;
        break;
      }
    }
    if (has_na) {
      for (int j = 0; j < 19; j++) res_p[i + j * n] = NA_REAL;
      continue;
    }

    // ── Part A (vectorized primitives) ────────────────────────────────────────

    // Temperature scalars
    double t_sum = 0.0, t_sumsq = 0.0;
    double tmx_max = tmx[0], tmn_min = tmn[0];
    double diurnal_sum = 0.0;
    for (int m = 0; m < 12; m++) {
      t_sum      += t[m];
      t_sumsq    += t[m] * t[m];
      if (tmx[m] > tmx_max) tmx_max = tmx[m];
      if (tmn[m] < tmn_min) tmn_min = tmn[m];
      diurnal_sum += tmx[m] - tmn[m];
    }
    double b01 = t_sum / 12.0;                              // BIO01
    double b02 = diurnal_sum / 12.0;                        // BIO02
    double b05 = tmx_max;                                   // BIO05
    double b06 = tmn_min;                                   // BIO06
    double b07 = b05 - b06;                                 // BIO07
    double b03 = (b07 > 0.0) ? 100.0 * b02 / b07 : BIO03_ZERO_RANGE_FALLBACK; // BIO03
    double b04 = 100.0 * sd_pop_xt(t);                      // BIO04

    // Precipitation scalars
    double p_sum = 0.0, p_max = p[0], p_min = p[0];
    for (int m = 0; m < 12; m++) {
      p_sum += p[m];
      if (p[m] > p_max) p_max = p[m];
      if (p[m] < p_min) p_min = p[m];
    }
    double b12 = p_sum;  // BIO12
    double b13 = p_max;  // BIO13
    double b14 = p_min;  // BIO14

    // BIO15: Precipitation Seasonality (CV = 100 * sd / mean; NaN when mean == 0)
    double p_mean = p_sum / 12.0;
    double p_ssq  = 0.0;
    for (int m = 0; m < 12; m++) {
      double d = p[m] - p_mean;
      p_ssq += d * d;
    }
    double b15 = (p_mean == 0.0)
                   ? R_NaN
                   : 100.0 * std::sqrt(p_ssq / 12.0) / p_mean;  // BIO15

    // Rolling quarter sums for temperature and precipitation
    double pr_qs[12], t_qs[12];
    rolling_quarter_sum_xt(p, pr_qs);   // precipitation quarter sums
    rolling_quarter_sum_xt(t, t_qs);    // temperature quarter sums

    // Quarter extremes for precipitation
    int wet_q  = argmax12(pr_qs);  // wettest quarter
    int dry_q  = argmin12(pr_qs);  // driest quarter
    double b16 = pr_qs[wet_q];     // BIO16
    double b17 = pr_qs[dry_q];     // BIO17

    // Quarter extremes for temperature
    int warm_q = argmax12(t_qs);                       // warmest quarter
    int cold_q = argmin12(t_qs);                       // coldest quarter
    double b10 = t_qs[warm_q] / 3.0;                  // BIO10
    double b11 = t_qs[cold_q] / 3.0;                  // BIO11

    // BIO08/BIO09: Mean Temp of Wettest/Driest Quarter
    double b08 = t_qs[wet_q] / 3.0;                   // BIO08
    double b09 = t_qs[dry_q] / 3.0;                   // BIO09

    // BIO18/BIO19: Precip of Warmest/Coldest Quarter
    double b18 = pr_qs[warm_q];                        // BIO18
    double b19 = pr_qs[cold_q];                        // BIO19

    // ── Zero-copy output: write directly to R matrix memory ──────────────────
    res_p[i +  0 * n] = b01;
    res_p[i +  1 * n] = b02;
    res_p[i +  2 * n] = b03;
    res_p[i +  3 * n] = b04;
    res_p[i +  4 * n] = b05;
    res_p[i +  5 * n] = b06;
    res_p[i +  6 * n] = b07;
    res_p[i +  7 * n] = b08;
    res_p[i +  8 * n] = b09;
    res_p[i +  9 * n] = b10;
    res_p[i + 10 * n] = b11;
    res_p[i + 11 * n] = b12;
    res_p[i + 12 * n] = b13;
    res_p[i + 13 * n] = b14;
    res_p[i + 14 * n] = b15;
    res_p[i + 15 * n] = b16;
    res_p[i + 16 * n] = b17;
    res_p[i + 17 * n] = b18;
    res_p[i + 18 * n] = b19;
  }

#ifdef _OPENMP
  omp_set_num_threads(prev_threads);
#endif

  return result;
}
