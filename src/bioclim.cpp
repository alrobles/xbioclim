#ifdef _OPENMP
#include <omp.h>
#endif
#include "xbioclim_omp.h"
#include <Rcpp.h>
#include <array>
#include <cstddef>
#include <limits>
#include <cmath>
#include "xbioclim.h"
using namespace Rcpp;

// ── Primitive helpers (NumericVector) ────────────────────────────────────────

// Population standard deviation (denominator N, matching xbioclim convention)
static inline double sd_pop_c(const NumericVector& x) {
  int n = x.size();
  double s = 0.0, ss = 0.0;
  for (int i = 0; i < n; i++) {
    s  += x[i];
    ss += x[i] * x[i];
  }
  double m = s / n;
  return std::sqrt(ss / n - m * m);
}

// Rolling quarter sums for 12 monthly values (circular wrapping)
static inline NumericVector rolling_quarter_sum_c(const NumericVector& x) {
  NumericVector result(12);
  for (int i = 0; i < 12; i++) {
    result[i] = x[i] + x[(i + 1) % 12] + x[(i + 2) % 12];
  }
  return result;
}

// 0-based index of quarter with maximum sum
static inline int quarter_argmax_c(const NumericVector& x) {
  NumericVector qs = rolling_quarter_sum_c(x);
  return which_max(qs);
}

// 0-based index of quarter with minimum sum
static inline int quarter_argmin_c(const NumericVector& x) {
  NumericVector qs = rolling_quarter_sum_c(x);
  return which_min(qs);
}

// Sum of 3-month quarter starting at 0-based index
static inline double quarter_sum_c(const NumericVector& x, int start) {
  return x[start] + x[(start + 1) % 12] + x[(start + 2) % 12];
}

// Mean of 3-month quarter starting at 0-based index
static inline double quarter_mean_c(const NumericVector& x, int start) {
  return quarter_sum_c(x, start) / 3.0;
}

// ── Raw double* helpers (OpenMP-safe, no heap allocation) ────────────────────

// Population standard deviation (denominator N) for a 12-element array
static inline double sd_pop_ptr(const double* x) {
  double s = 0.0, ss = 0.0;
  for (int i = 0; i < 12; i++) {
    s  += x[i];
    ss += x[i] * x[i];
  }
  double m = s / 12.0;
  return std::sqrt(ss / 12.0 - m * m);
}

// Sum of 3-month quarter starting at 0-based index (circular)
static inline double quarter_sum_ptr(const double* x, int start) {
  return x[start] + x[(start + 1) % 12] + x[(start + 2) % 12];
}

// Mean of 3-month quarter starting at 0-based index
static inline double quarter_mean_ptr(const double* x, int start) {
  return quarter_sum_ptr(x, start) / 3.0;
}

// 0-based index of quarter with maximum rolling sum
static inline int quarter_argmax_ptr(const double* x) {
  int best = 0;
  double best_sum = x[0] + x[1] + x[2];
  for (int k = 1; k < 12; k++) {
    double s = x[k] + x[(k + 1) % 12] + x[(k + 2) % 12];
    if (s > best_sum) { best_sum = s; best = k; }
  }
  return best;
}

// 0-based index of quarter with minimum rolling sum
static inline int quarter_argmin_ptr(const double* x) {
  int best = 0;
  double best_sum = x[0] + x[1] + x[2];
  for (int k = 1; k < 12; k++) {
    double s = x[k] + x[(k + 1) % 12] + x[(k + 2) % 12];
    if (s < best_sum) { best_sum = s; best = k; }
  }
  return best;
}

// ── BIO01: Mean Annual Temperature ───────────────────────────────────────────

//' Compute BIO01 (Mean Annual Temperature) for a raster block
//'
//' @param tas Numeric matrix with 12 columns (one per month); rows are pixels.
//' @return Numeric vector with one value per pixel.
//' @keywords internal
// [[Rcpp::export]]
NumericVector bio01_cpp(NumericMatrix tas) {
  int n = tas.nrow();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = mean(NumericVector(tas.row(i)));
  }
  return result;
}

// ── BIO02: Mean Diurnal Range ─────────────────────────────────────────────────

//' Compute BIO02 (Mean Diurnal Range) for a raster block
//'
//' @param tasmax Numeric matrix (pixels x 12): monthly max temperature.
//' @param tasmin Numeric matrix (pixels x 12): monthly min temperature.
//' @return Numeric vector with one value per pixel.
//' @keywords internal
// [[Rcpp::export]]
NumericVector bio02_cpp(NumericMatrix tasmax, NumericMatrix tasmin) {
  int n = tasmax.nrow();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = mean(NumericVector(tasmax.row(i)) - NumericVector(tasmin.row(i)));
  }
  return result;
}

// ── BIO03: Isothermality (100 * BIO02 / BIO07) ───────────────────────────────

//' Compute BIO03 (Isothermality) for a raster block
//'
//' @param tasmax Numeric matrix (pixels x 12): monthly max temperature.
//' @param tasmin Numeric matrix (pixels x 12): monthly min temperature.
//' @return Numeric vector with one value per pixel (NaN where BIO07 == 0).
//' @keywords internal
// [[Rcpp::export]]
NumericVector bio03_cpp(NumericMatrix tasmax, NumericMatrix tasmin) {
  int n = tasmax.nrow();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    NumericVector mx(tasmax.row(i));
    NumericVector mn(tasmin.row(i));
    double b02 = mean(mx - mn);
    double b07 = max(mx) - min(mn);
    result[i] = (b07 == 0.0) ? R_NaN : 100.0 * b02 / b07;
  }
  return result;
}

// ── BIO04: Temperature Seasonality ───────────────────────────────────────────

//' Compute BIO04 (Temperature Seasonality) for a raster block
//'
//' @param tas Numeric matrix (pixels x 12): monthly mean temperature.
//' @return Numeric vector with one value per pixel.
//' @keywords internal
// [[Rcpp::export]]
NumericVector bio04_cpp(NumericMatrix tas) {
  int n = tas.nrow();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = 100.0 * sd_pop_c(NumericVector(tas.row(i)));
  }
  return result;
}

// ── BIO05: Max Temperature of Warmest Month ──────────────────────────────────

//' Compute BIO05 (Max Temperature of Warmest Month) for a raster block
//'
//' @param tasmax Numeric matrix (pixels x 12): monthly max temperature.
//' @return Numeric vector with one value per pixel.
//' @keywords internal
// [[Rcpp::export]]
NumericVector bio05_cpp(NumericMatrix tasmax) {
  int n = tasmax.nrow();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = max(NumericVector(tasmax.row(i)));
  }
  return result;
}

// ── BIO06: Min Temperature of Coldest Month ──────────────────────────────────

//' Compute BIO06 (Min Temperature of Coldest Month) for a raster block
//'
//' @param tasmin Numeric matrix (pixels x 12): monthly min temperature.
//' @return Numeric vector with one value per pixel.
//' @keywords internal
// [[Rcpp::export]]
NumericVector bio06_cpp(NumericMatrix tasmin) {
  int n = tasmin.nrow();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = min(NumericVector(tasmin.row(i)));
  }
  return result;
}

// ── BIO07: Temperature Annual Range (BIO05 - BIO06) ──────────────────────────

//' Compute BIO07 (Temperature Annual Range) for a raster block
//'
//' @param tasmax Numeric matrix (pixels x 12): monthly max temperature.
//' @param tasmin Numeric matrix (pixels x 12): monthly min temperature.
//' @return Numeric vector with one value per pixel.
//' @keywords internal
// [[Rcpp::export]]
NumericVector bio07_cpp(NumericMatrix tasmax, NumericMatrix tasmin) {
  int n = tasmax.nrow();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = max(NumericVector(tasmax.row(i))) - min(NumericVector(tasmin.row(i)));
  }
  return result;
}

// ── BIO08: Mean Temperature of Wettest Quarter ───────────────────────────────

//' Compute BIO08 (Mean Temperature of Wettest Quarter) for a raster block
//'
//' @param tas Numeric matrix (pixels x 12): monthly mean temperature.
//' @param pr  Numeric matrix (pixels x 12): monthly precipitation.
//' @return Numeric vector with one value per pixel.
//' @keywords internal
// [[Rcpp::export]]
NumericVector bio08_cpp(NumericMatrix tas, NumericMatrix pr) {
  int n = tas.nrow();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    int ws = quarter_argmax_c(NumericVector(pr.row(i)));
    result[i] = quarter_mean_c(NumericVector(tas.row(i)), ws);
  }
  return result;
}

// ── BIO09: Mean Temperature of Driest Quarter ────────────────────────────────

//' Compute BIO09 (Mean Temperature of Driest Quarter) for a raster block
//'
//' @param tas Numeric matrix (pixels x 12): monthly mean temperature.
//' @param pr  Numeric matrix (pixels x 12): monthly precipitation.
//' @return Numeric vector with one value per pixel.
//' @keywords internal
// [[Rcpp::export]]
NumericVector bio09_cpp(NumericMatrix tas, NumericMatrix pr) {
  int n = tas.nrow();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    int ds = quarter_argmin_c(NumericVector(pr.row(i)));
    result[i] = quarter_mean_c(NumericVector(tas.row(i)), ds);
  }
  return result;
}

// ── BIO10: Mean Temperature of Warmest Quarter ───────────────────────────────

//' Compute BIO10 (Mean Temperature of Warmest Quarter) for a raster block
//'
//' @param tas Numeric matrix (pixels x 12): monthly mean temperature.
//' @return Numeric vector with one value per pixel.
//' @keywords internal
// [[Rcpp::export]]
NumericVector bio10_cpp(NumericMatrix tas) {
  int n = tas.nrow();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    NumericVector row(tas.row(i));
    int ws = quarter_argmax_c(row);
    result[i] = quarter_mean_c(row, ws);
  }
  return result;
}

// ── BIO11: Mean Temperature of Coldest Quarter ───────────────────────────────

//' Compute BIO11 (Mean Temperature of Coldest Quarter) for a raster block
//'
//' @param tas Numeric matrix (pixels x 12): monthly mean temperature.
//' @return Numeric vector with one value per pixel.
//' @keywords internal
// [[Rcpp::export]]
NumericVector bio11_cpp(NumericMatrix tas) {
  int n = tas.nrow();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    NumericVector row(tas.row(i));
    int cs = quarter_argmin_c(row);
    result[i] = quarter_mean_c(row, cs);
  }
  return result;
}

// ── BIO12: Annual Precipitation ──────────────────────────────────────────────

//' Compute BIO12 (Annual Precipitation) for a raster block
//'
//' @param pr Numeric matrix (pixels x 12): monthly precipitation.
//' @return Numeric vector with one value per pixel.
//' @keywords internal
// [[Rcpp::export]]
NumericVector bio12_cpp(NumericMatrix pr) {
  int n = pr.nrow();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = sum(NumericVector(pr.row(i)));
  }
  return result;
}

// ── BIO13: Precipitation of Wettest Month ────────────────────────────────────

//' Compute BIO13 (Precipitation of Wettest Month) for a raster block
//'
//' @param pr Numeric matrix (pixels x 12): monthly precipitation.
//' @return Numeric vector with one value per pixel.
//' @keywords internal
// [[Rcpp::export]]
NumericVector bio13_cpp(NumericMatrix pr) {
  int n = pr.nrow();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = max(NumericVector(pr.row(i)));
  }
  return result;
}

// ── BIO14: Precipitation of Driest Month ─────────────────────────────────────

//' Compute BIO14 (Precipitation of Driest Month) for a raster block
//'
//' @param pr Numeric matrix (pixels x 12): monthly precipitation.
//' @return Numeric vector with one value per pixel.
//' @keywords internal
// [[Rcpp::export]]
NumericVector bio14_cpp(NumericMatrix pr) {
  int n = pr.nrow();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = min(NumericVector(pr.row(i)));
  }
  return result;
}

// ── BIO15: Precipitation Seasonality (CV) ────────────────────────────────────

//' Compute BIO15 (Precipitation Seasonality) for a raster block
//'
//' @param pr Numeric matrix (pixels x 12): monthly precipitation.
//' @return Numeric vector with one value per pixel (NaN where mean precip == 0).
//' @keywords internal
// [[Rcpp::export]]
NumericVector bio15_cpp(NumericMatrix pr) {
  int n = pr.nrow();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    NumericVector row(pr.row(i));
    double pr_mean = mean(row);
    result[i] = (pr_mean == 0.0) ? R_NaN : 100.0 * sd_pop_c(row) / pr_mean;
  }
  return result;
}

// ── BIO16: Precipitation of Wettest Quarter ──────────────────────────────────

//' Compute BIO16 (Precipitation of Wettest Quarter) for a raster block
//'
//' @param pr Numeric matrix (pixels x 12): monthly precipitation.
//' @return Numeric vector with one value per pixel.
//' @keywords internal
// [[Rcpp::export]]
NumericVector bio16_cpp(NumericMatrix pr) {
  int n = pr.nrow();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    NumericVector row(pr.row(i));
    int ws = quarter_argmax_c(row);
    result[i] = quarter_sum_c(row, ws);
  }
  return result;
}

// ── BIO17: Precipitation of Driest Quarter ───────────────────────────────────

//' Compute BIO17 (Precipitation of Driest Quarter) for a raster block
//'
//' @param pr Numeric matrix (pixels x 12): monthly precipitation.
//' @return Numeric vector with one value per pixel.
//' @keywords internal
// [[Rcpp::export]]
NumericVector bio17_cpp(NumericMatrix pr) {
  int n = pr.nrow();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    NumericVector row(pr.row(i));
    int ds = quarter_argmin_c(row);
    result[i] = quarter_sum_c(row, ds);
  }
  return result;
}

// ── BIO18: Precipitation of Warmest Quarter ──────────────────────────────────

//' Compute BIO18 (Precipitation of Warmest Quarter) for a raster block
//'
//' @param tas Numeric matrix (pixels x 12): monthly mean temperature.
//' @param pr  Numeric matrix (pixels x 12): monthly precipitation.
//' @return Numeric vector with one value per pixel.
//' @keywords internal
// [[Rcpp::export]]
NumericVector bio18_cpp(NumericMatrix tas, NumericMatrix pr) {
  int n = tas.nrow();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    int ws = quarter_argmax_c(NumericVector(tas.row(i)));
    result[i] = quarter_sum_c(NumericVector(pr.row(i)), ws);
  }
  return result;
}

// ── BIO19: Precipitation of Coldest Quarter ──────────────────────────────────

//' Compute BIO19 (Precipitation of Coldest Quarter) for a raster block
//'
//' @param tas Numeric matrix (pixels x 12): monthly mean temperature.
//' @param pr  Numeric matrix (pixels x 12): monthly precipitation.
//' @return Numeric vector with one value per pixel.
//' @keywords internal
// [[Rcpp::export]]
NumericVector bio19_cpp(NumericMatrix tas, NumericMatrix pr) {
  int n = tas.nrow();
  NumericVector result(n);
  for (int i = 0; i < n; i++) {
    int cs = quarter_argmin_c(NumericVector(tas.row(i)));
    result[i] = quarter_sum_c(NumericVector(pr.row(i)), cs);
  }
  return result;
}

namespace {

// Helpers for na.rm = TRUE: operate on raw 12-element arrays.
// R/BioclimData.R is the reference for semantics.

static const double XB_NAN = std::numeric_limits<double>::quiet_NaN();

inline bool all_na_r(const double* x, int n) {
  for (int i = 0; i < n; ++i) {
    if (!ISNA(x[i])) return false;
  }
  return true;
}

inline int count_valid_r(const double* x, int n) {
  int c = 0;
  for (int i = 0; i < n; ++i) if (!ISNA(x[i])) ++c;
  return c;
}

inline double sum_valid_or_nan(const double* x, int n) {
  double s = 0.0;
  int c = 0;
  for (int i = 0; i < n; ++i) {
    if (!ISNA(x[i])) { s += x[i]; ++c; }
  }
  return c == 0 ? XB_NAN : s;
}

inline double mean_valid_or_nan(const double* x, int n) {
  double s = 0.0;
  int c = 0;
  for (int i = 0; i < n; ++i) {
    if (!ISNA(x[i])) { s += x[i]; ++c; }
  }
  return c == 0 ? XB_NAN : s / c;
}

inline double max_valid_or_nan(const double* x, int n) {
  bool first = true;
  double best = 0.0;
  for (int i = 0; i < n; ++i) {
    if (!ISNA(x[i])) {
      if (first || x[i] > best) { best = x[i]; first = false; }
    }
  }
  return first ? XB_NAN : best;
}

inline double min_valid_or_nan(const double* x, int n) {
  bool first = true;
  double best = 0.0;
  for (int i = 0; i < n; ++i) {
    if (!ISNA(x[i])) {
      if (first || x[i] < best) { best = x[i]; first = false; }
    }
  }
  return first ? XB_NAN : best;
}

// Population standard deviation of non-NA values; NaN if none.
inline double sd_pop_valid(const double* x, int n) {
  int c = 0;
  for (int i = 0; i < n; ++i) if (!ISNA(x[i])) ++c;
  if (c == 0) return XB_NAN;
  double s = 0.0;
  for (int i = 0; i < n; ++i) if (!ISNA(x[i])) s += x[i];
  double m = s / c;
  double sq = 0.0;
  for (int i = 0; i < n; ++i) if (!ISNA(x[i])) { double d = x[i] - m; sq += d * d; }
  return std::sqrt(sq / c);
}

// Sum of a circular 3-month quarter starting at 0-based month, skipping NA.
// Returns 0 if no month is valid (matches R sum(..., na.rm = TRUE)).
inline double quarter_sum_valid(const double* x, int start) {
  double s = 0.0;
  for (int k = 0; k < 3; ++k) {
    double v = x[(start + k) % 12];
    if (!ISNA(v)) s += v;
  }
  return s;
}

// Mean of a circular 3-month quarter starting at 0-based month, skipping NA.
// Returns NaN if no month is valid (matches R mean(..., na.rm = TRUE)).
inline double quarter_mean_valid(const double* x, int start) {
  double s = 0.0;
  int c = 0;
  for (int k = 0; k < 3; ++k) {
    double v = x[(start + k) % 12];
    if (!ISNA(v)) { s += v; ++c; }
  }
  return c == 0 ? XB_NAN : s / c;
}

// Quarter with maximum sum of valid months; quarters with 0 valid months get -Inf.
// If all quarters are empty, returns 0 (matches which.max on all -Inf).
inline int quarter_argmax_valid(const double* x) {
  int best = 0;
  double best_sum = std::numeric_limits<double>::lowest();
  for (int start = 0; start < 12; ++start) {
    double s = 0.0;
    int c = 0;
    for (int k = 0; k < 3; ++k) {
      double v = x[(start + k) % 12];
      if (!ISNA(v)) { s += v; ++c; }
    }
    if (c == 0) s = std::numeric_limits<double>::lowest();
    if (s > best_sum) { best_sum = s; best = start; }
  }
  return best;
}

// Quarter with minimum sum of valid months; quarters with 0 valid months get +Inf.
inline int quarter_argmin_valid(const double* x) {
  int best = 0;
  double best_sum = std::numeric_limits<double>::max();
  for (int start = 0; start < 12; ++start) {
    double s = 0.0;
    int c = 0;
    for (int k = 0; k < 3; ++k) {
      double v = x[(start + k) % 12];
      if (!ISNA(v)) { s += v; ++c; }
    }
    if (c == 0) s = std::numeric_limits<double>::max();
    if (s < best_sum) { best_sum = s; best = start; }
  }
  return best;
}

// Single-pixel computation for na.rm = TRUE, following R/BioclimData.R.
inline std::array<double, 19> compute_pixel_na_rm(const double* tas,
                                                  const double* tasmax,
                                                  const double* tasmin,
                                                  const double* pr) {
  // b01: mean annual temperature of valid months
  double b01 = all_na_r(tas, 12) ? XB_NAN : mean_valid_or_nan(tas, 12);

  // b02: mean diurnal range (tasmax - tasmin) over valid months
  double diff_sum = 0.0;
  int diff_count = 0;
  for (int m = 0; m < 12; ++m) {
    if (!ISNA(tasmax[m]) && !ISNA(tasmin[m])) {
      diff_sum += tasmax[m] - tasmin[m];
      ++diff_count;
    }
  }
  double b02 = diff_count == 0 ? XB_NAN : diff_sum / diff_count;

  // b05-b07
  double b05 = max_valid_or_nan(tasmax, 12);
  double b06 = min_valid_or_nan(tasmin, 12);
  double b07 = b05 - b06;
  double b03 = (std::isnan(b07) || b07 == 0.0) ? XB_NAN : 100.0 * b02 / b07;

  // b04: temperature seasonality
  double b04 = 100.0 * sd_pop_valid(tas, 12);

  // b12-b15
  bool pr_all_na = all_na_r(pr, 12);
  double b12 = pr_all_na ? XB_NAN : sum_valid_or_nan(pr, 12);
  double b13 = pr_all_na ? XB_NAN : max_valid_or_nan(pr, 12);
  double b14 = pr_all_na ? XB_NAN : min_valid_or_nan(pr, 12);
  double pr_mean = pr_all_na ? XB_NAN : mean_valid_or_nan(pr, 12);
  double b15 = (std::isnan(pr_mean) || pr_mean <= 0.0)
                 ? XB_NAN
                 : 100.0 * sd_pop_valid(pr, 12) / pr_mean;

  // Quarter indices (0-based)
  int wet_q  = quarter_argmax_valid(pr);
  int dry_q  = quarter_argmin_valid(pr);
  int warm_q = quarter_argmax_valid(tas);
  int cold_q = quarter_argmin_valid(tas);

  // b08-b11: mean temperature of the selected quarter
  double b08 = quarter_mean_valid(tas, wet_q);
  double b09 = quarter_mean_valid(tas, dry_q);
  double b10 = quarter_mean_valid(tas, warm_q);
  double b11 = quarter_mean_valid(tas, cold_q);

  // b16-b19: precipitation of selected quarter
  double b16 = pr_all_na ? XB_NAN : quarter_sum_valid(pr, wet_q);
  double b17 = pr_all_na ? XB_NAN : quarter_sum_valid(pr, dry_q);
  double b18 = pr_all_na ? XB_NAN : quarter_sum_valid(pr, warm_q);
  double b19 = pr_all_na ? XB_NAN : quarter_sum_valid(pr, cold_q);

  return {b01, b02, b03, b04, b05, b06, b07,
          b08, b09, b10, b11,
          b12, b13, b14, b15,
          b16, b17, b18, b19};
}

} // namespace

// ── Batch: All 19 variables ───────────────────────────────────────────────────

//' Compute all 19 bioclimatic variables for a raster block
//'
//' @param tas    Numeric matrix (pixels x 12): monthly mean temperature.
//' @param tasmax Numeric matrix (pixels x 12): monthly max temperature.
//' @param tasmin Numeric matrix (pixels x 12): monthly min temperature.
//' @param pr     Numeric matrix (pixels x 12): monthly precipitation.
//' @param ncores Integer: number of OpenMP threads (default 1).
//' @param na_rm  Logical: if TRUE, treat NA as missing and compute each BIO
//'   from the available months (quarters need >=1 valid month). If FALSE,
//'   a single NA in any input for a pixel gives an all-NA row (default).
//' @return Numeric matrix (pixels x 19) with one column per variable
//'   (bio01..bio19), named accordingly.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix bioclim_cpp(NumericMatrix tas,
                          NumericMatrix tasmax,
                          NumericMatrix tasmin,
                          NumericMatrix pr,
                          int ncores = 1,
                          bool na_rm = false) {
  int n = tas.nrow();
  NumericMatrix result(n, 19);

  CharacterVector cnames = CharacterVector::create(
    "bio01", "bio02", "bio03", "bio04", "bio05",
    "bio06", "bio07", "bio08", "bio09", "bio10",
    "bio11", "bio12", "bio13", "bio14", "bio15",
    "bio16", "bio17", "bio18", "bio19"
  );
  colnames(result) = cnames;

#ifdef _OPENMP
  int prev_threads = omp_get_max_threads();
  omp_set_num_threads(xbioclim_safe_threads(ncores));
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < n; i++) {
    // Gather one pixel into contiguous stack arrays.
    double t[12], tmx[12], tmn[12], p[12];
    bool has_na = false;
    for (int m = 0; m < 12; m++) {
      t[m]   = tas(i, m);
      tmx[m] = tasmax(i, m);
      tmn[m] = tasmin(i, m);
      p[m]   = pr(i, m);
      if (ISNA(t[m]) || ISNA(tmx[m]) || ISNA(tmn[m]) || ISNA(p[m])) {
        has_na = true;
      }
    }

    std::array<double, 19> bio;
    if (!na_rm && has_na) {
      // Conservative semantics: a single NA makes the whole pixel NA.
      for (int j = 0; j < 19; ++j) bio[j] = NA_REAL;
    } else if (na_rm) {
      bio = compute_pixel_na_rm(t, tmx, tmn, p);
    } else {
      bio = xbioclim::compute_pixel(t, tmx, tmn, p);
    }

    for (int j = 0; j < 19; j++) {
      double v = bio[j];
      result(i, j) = std::isnan(v) ? NA_REAL : v;
    }
  }

#ifdef _OPENMP
  omp_set_num_threads(prev_threads);
#endif

  return result;
}

// ── Quarterly / seasonal variables for arbitrary months ──────────────────────

static inline bool has_na_sel(const NumericVector& x, const IntegerVector& months) {
  for (int j = 0; j < months.size(); j++) {
    if (NumericVector::is_na(x[months[j] - 1])) return true;
  }
  return false;
}

static inline int n_valid_sel(const NumericVector& x, const IntegerVector& months) {
  int n = 0;
  for (int j = 0; j < months.size(); j++) {
    if (!NumericVector::is_na(x[months[j] - 1])) n++;
  }
  return n;
}

static inline double mean_sel(const NumericVector& x, const IntegerVector& months, bool na_rm) {
  double s = 0.0;
  int n = 0;
  for (int j = 0; j < months.size(); j++) {
    double v = x[months[j] - 1];
    if (na_rm && NumericVector::is_na(v)) continue;
    s += v;
    n++;
  }
  if (n == 0) return R_NaN;
  return s / n;
}

static inline double sum_sel(const NumericVector& x, const IntegerVector& months, bool na_rm) {
  double s = 0.0;
  int n = 0;
  for (int j = 0; j < months.size(); j++) {
    double v = x[months[j] - 1];
    if (na_rm && NumericVector::is_na(v)) continue;
    s += v;
    n++;
  }
  if (n == 0) return R_NaN;
  return s;
}

static inline double max_sel(const NumericVector& x, const IntegerVector& months, bool na_rm) {
  bool first = true;
  double m = 0.0;
  for (int j = 0; j < months.size(); j++) {
    double v = x[months[j] - 1];
    if (na_rm && NumericVector::is_na(v)) continue;
    if (first || v > m) { m = v; first = false; }
  }
  if (first) return R_NaN;
  return m;
}

static inline double min_sel(const NumericVector& x, const IntegerVector& months, bool na_rm) {
  bool first = true;
  double m = 0.0;
  for (int j = 0; j < months.size(); j++) {
    double v = x[months[j] - 1];
    if (na_rm && NumericVector::is_na(v)) continue;
    if (first || v < m) { m = v; first = false; }
  }
  if (first) return R_NaN;
  return m;
}

static inline double sd_pop_sel(const NumericVector& x, const IntegerVector& months, bool na_rm) {
  double s = 0.0, ss = 0.0;
  int n = 0;
  for (int j = 0; j < months.size(); j++) {
    double v = x[months[j] - 1];
    if (na_rm && NumericVector::is_na(v)) continue;
    s += v;
    ss += v * v;
    n++;
  }
  if (n == 0) return R_NaN;
  double m = s / n;
  return std::sqrt(ss / n - m * m);
}

//' Compute quarterly/seasonal climate variables for a raster block
//'
//' @param tas    Numeric matrix (pixels x 12): monthly mean temperature.
//' @param tasmax Numeric matrix (pixels x 12): monthly maximum temperature.
//' @param tasmin Numeric matrix (pixels x 12): monthly minimum temperature.
//' @param pr     Numeric matrix (pixels x 12): monthly precipitation.
//' @param months Integer vector of 1-based month indices to include.
//' @param na_rm  Logical: if `TRUE`, skip `NA` months.
//' @return Numeric matrix (pixels x 6) with columns
//'   `tmean_s`, `tmax_max`, `tmin_min`, `trange`, `pr_tot`, `pr_cv`.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix quarterly_variables_cpp(NumericMatrix tas,
                                      NumericMatrix tasmax,
                                      NumericMatrix tasmin,
                                      NumericMatrix pr,
                                      IntegerVector months,
                                      bool na_rm = false) {
  if (tas.nrow() != tasmax.nrow() || tas.nrow() != tasmin.nrow() ||
      tas.nrow() != pr.nrow()) {
    stop("All input matrices must have the same number of rows");
  }
  if (tas.ncol() != 12 || tasmax.ncol() != 12 || tasmin.ncol() != 12 ||
      pr.ncol() != 12) {
    stop("All input matrices must have 12 columns (one per month)");
  }
  for (int j = 0; j < months.size(); j++) {
    if (months[j] < 1 || months[j] > 12) {
      stop("month indices must be between 1 and 12");
    }
  }

  int n = tas.nrow();
  NumericMatrix result(n, 6);
  colnames(result) = CharacterVector::create(
    "tmean_s", "tmax_max", "tmin_min", "trange", "pr_tot", "pr_cv");

  for (int i = 0; i < n; i++) {
    NumericVector r_tas    = tas.row(i);
    NumericVector r_tasmax = tasmax.row(i);
    NumericVector r_tasmin = tasmin.row(i);
    NumericVector r_pr     = pr.row(i);

    if (!na_rm && (has_na_sel(r_tas, months) || has_na_sel(r_tasmax, months) ||
                   has_na_sel(r_tasmin, months) || has_na_sel(r_pr, months))) {
      for (int j = 0; j < 6; j++) result(i, j) = NA_REAL;
      continue;
    }

    double tmean_s = mean_sel(r_tas, months, na_rm);
    double tmax_m  = max_sel(r_tasmax, months, na_rm);
    double tmin_m  = min_sel(r_tasmin, months, na_rm);
    double pr_sum  = sum_sel(r_pr, months, na_rm);
    double pr_mean = mean_sel(r_pr, months, na_rm);

    result(i, 0) = tmean_s;
    result(i, 1) = tmax_m;
    result(i, 2) = tmin_m;
    result(i, 3) = std::isnan(tmax_m) || std::isnan(tmin_m)
                     ? R_NaN : tmax_m - tmin_m;
    result(i, 4) = pr_sum;
    result(i, 5) = (std::isnan(pr_mean) || pr_mean <= 0.0)
                     ? R_NaN : 100.0 * sd_pop_sel(r_pr, months, na_rm) / pr_mean;
  }

  return result;
}

// ── Window / rolling bioclimatic variables ───────────────────────────────────

// Build a 0-based mask for selected months (true = selected).
static inline void build_month_mask(const IntegerVector& months, bool mask[12]) {
  for (int i = 0; i < 12; i++) mask[i] = false;
  for (int j = 0; j < months.size(); j++) {
    if (months[j] >= 1 && months[j] <= 12) mask[months[j] - 1] = true;
  }
}

// Sum over a contiguous window starting at 1-based `start` (length `win`).
static inline double window_sum_sel(const NumericVector& x,
                                    int start, int win, bool na_rm) {
  double s = 0.0;
  int n = 0;
  for (int i = 0; i < win; i++) {
    int idx = (start - 1 + i) % 12;
    double v = x[idx];
    if (na_rm && NumericVector::is_na(v)) continue;
    s += v;
    n++;
  }
  if (n == 0) return R_NaN;
  return s;
}

// Mean over a contiguous window starting at 1-based `start`.
static inline double window_mean_sel(const NumericVector& x,
                                     int start, int win, bool na_rm) {
  double s = 0.0;
  int n = 0;
  for (int i = 0; i < win; i++) {
    int idx = (start - 1 + i) % 12;
    double v = x[idx];
    if (na_rm && NumericVector::is_na(v)) continue;
    s += v;
    n++;
  }
  if (n == 0) return R_NaN;
  return s / n;
}

// Find the starting month of the valid contiguous window (length `win`)
// within `months` with the maximum sum of `x`.
static inline int argmax_window_sum(const NumericVector& x,
                                    const IntegerVector& months,
                                    int win,
                                    const bool mask[12],
                                    bool na_rm) {
  int best_start = -1;
  double best_sum = 0.0;
  bool first = true;
  for (int j = 0; j < months.size(); j++) {
    int start = months[j];
    bool valid = true;
    for (int i = 0; i < win; i++) {
      int idx = (start - 1 + i) % 12;
      if (!mask[idx]) { valid = false; break; }
    }
    if (!valid) continue;

    double s = window_sum_sel(x, start, win, na_rm);
    if (NumericVector::is_na(s)) continue;
    if (first || s > best_sum) {
      best_sum = s;
      best_start = start;
      first = false;
    }
  }
  return best_start;
}

// Minimum sum.
static inline int argmin_window_sum(const NumericVector& x,
                                    const IntegerVector& months,
                                    int win,
                                    const bool mask[12],
                                    bool na_rm) {
  int best_start = -1;
  double best_sum = 0.0;
  bool first = true;
  for (int j = 0; j < months.size(); j++) {
    int start = months[j];
    bool valid = true;
    for (int i = 0; i < win; i++) {
      int idx = (start - 1 + i) % 12;
      if (!mask[idx]) { valid = false; break; }
    }
    if (!valid) continue;

    double s = window_sum_sel(x, start, win, na_rm);
    if (NumericVector::is_na(s)) continue;
    if (first || s < best_sum) {
      best_sum = s;
      best_start = start;
      first = false;
    }
  }
  return best_start;
}

// Maximum mean.
static inline int argmax_window_mean(const NumericVector& x,
                                     const IntegerVector& months,
                                     int win,
                                     const bool mask[12],
                                     bool na_rm) {
  int best_start = -1;
  double best_mean = 0.0;
  bool first = true;
  for (int j = 0; j < months.size(); j++) {
    int start = months[j];
    bool valid = true;
    for (int i = 0; i < win; i++) {
      int idx = (start - 1 + i) % 12;
      if (!mask[idx]) { valid = false; break; }
    }
    if (!valid) continue;

    double m = window_mean_sel(x, start, win, na_rm);
    if (NumericVector::is_na(m)) continue;
    if (first || m > best_mean) {
      best_mean = m;
      best_start = start;
      first = false;
    }
  }
  return best_start;
}

// Minimum mean.
static inline int argmin_window_mean(const NumericVector& x,
                                     const IntegerVector& months,
                                     int win,
                                     const bool mask[12],
                                     bool na_rm) {
  int best_start = -1;
  double best_mean = 0.0;
  bool first = true;
  for (int j = 0; j < months.size(); j++) {
    int start = months[j];
    bool valid = true;
    for (int i = 0; i < win; i++) {
      int idx = (start - 1 + i) % 12;
      if (!mask[idx]) { valid = false; break; }
    }
    if (!valid) continue;

    double m = window_mean_sel(x, start, win, na_rm);
    if (NumericVector::is_na(m)) continue;
    if (first || m < best_mean) {
      best_mean = m;
      best_start = start;
      first = false;
    }
  }
  return best_start;
}

// Compute the 19 bioclimatic variables over a selected set of months.
static inline std::array<double, 19> compute_pixel_window(
    const NumericVector& tas,
    const NumericVector& tasmax,
    const NumericVector& tasmin,
    const NumericVector& pr,
    const IntegerVector& months,
    int window,
    bool na_rm) {

  std::array<double, 19> bio;
  bio.fill(NA_REAL);

  if (!na_rm && (has_na_sel(tas, months) || has_na_sel(tasmax, months) ||
                 has_na_sel(tasmin, months) || has_na_sel(pr, months))) {
    return bio;
  }

  // Base statistics over the selected months.
  double b01 = mean_sel(tas, months, na_rm);
  double b05 = max_sel(tasmax, months, na_rm);
  double b06 = min_sel(tasmin, months, na_rm);
  double b12 = sum_sel(pr, months, na_rm);
  double b13 = max_sel(pr, months, na_rm);
  double b14 = min_sel(pr, months, na_rm);
  double pr_mean = mean_sel(pr, months, na_rm);

  double b02 = mean_sel(tasmax, months, na_rm) - mean_sel(tasmin, months, na_rm);
  double b07 = (std::isnan(b05) || std::isnan(b06))
                 ? R_NaN : b05 - b06;
  double b03 = (std::isnan(b07) || b07 == 0.0) ? R_NaN : 100.0 * b02 / b07;
  double b04 = 100.0 * sd_pop_sel(tas, months, na_rm);
  double b15 = (std::isnan(pr_mean) || pr_mean <= 0.0)
                 ? R_NaN : 100.0 * sd_pop_sel(pr, months, na_rm) / pr_mean;

  bio[0]  = b01;
  bio[1]  = b02;
  bio[2]  = b03;
  bio[3]  = b04;
  bio[4]  = b05;
  bio[5]  = b06;
  bio[6]  = b07;
  bio[11] = b12;
  bio[12] = b13;
  bio[13] = b14;
  bio[14] = b15;

  if (months.size() >= window) {
    bool mask[12];
    build_month_mask(months, mask);

    int wet_start  = argmax_window_sum(pr, months, window, mask, na_rm);
    int dry_start  = argmin_window_sum(pr, months, window, mask, na_rm);
    int warm_start = argmax_window_mean(tas, months, window, mask, na_rm);
    int cold_start = argmin_window_mean(tas, months, window, mask, na_rm);

    if (wet_start > 0) {
      bio[7]  = window_mean_sel(tas, wet_start, window, na_rm);  // BIO08
      bio[15] = window_sum_sel(pr, wet_start, window, na_rm);     // BIO16
    }
    if (dry_start > 0) {
      bio[8]  = window_mean_sel(tas, dry_start, window, na_rm);  // BIO09
      bio[16] = window_sum_sel(pr, dry_start, window, na_rm);     // BIO17
    }
    if (warm_start > 0) {
      bio[9]  = window_mean_sel(tas, warm_start, window, na_rm); // BIO10
      bio[17] = window_sum_sel(pr, warm_start, window, na_rm);    // BIO18
    }
    if (cold_start > 0) {
      bio[10] = window_mean_sel(tas, cold_start, window, na_rm); // BIO11
      bio[18] = window_sum_sel(pr, cold_start, window, na_rm);    // BIO19
    }
  }

  return bio;
}

//' Compute bioclimatic variables over an arbitrary window of months
//'
//' @param tas    Numeric matrix (pixels x 12): monthly mean temperature.
//' @param tasmax Numeric matrix (pixels x 12): monthly maximum temperature.
//' @param tasmin Numeric matrix (pixels x 12): monthly minimum temperature.
//' @param pr     Numeric matrix (pixels x 12): monthly precipitation.
//' @param months Integer vector of 1-based month indices in the window.
//' @param window Integer: length (months) of the internal rolling sub-window
//'   used for the BIO08-BIO19 variables.  Must be >= 3 and <= length(months)
//'   for those variables to be non-NA.
//' @param na_rm  Logical: if `TRUE`, skip `NA` months.
//' @return Numeric matrix (pixels x 19) with columns `bio01`..`bio19`.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix bioclim_window_cpp(NumericMatrix tas,
                                 NumericMatrix tasmax,
                                 NumericMatrix tasmin,
                                 NumericMatrix pr,
                                 IntegerVector months,
                                 int window = 3,
                                 bool na_rm = false) {
  if (tas.nrow() != tasmax.nrow() || tas.nrow() != tasmin.nrow() ||
      tas.nrow() != pr.nrow()) {
    stop("All input matrices must have the same number of rows");
  }
  if (tas.ncol() != 12 || tasmax.ncol() != 12 || tasmin.ncol() != 12 ||
      pr.ncol() != 12) {
    stop("All input matrices must have 12 columns");
  }
  if (window < 2 || window > 12) {
    stop("'window' must be between 2 and 12");
  }
  for (int j = 0; j < months.size(); j++) {
    if (months[j] < 1 || months[j] > 12) {
      stop("month indices must be between 1 and 12");
    }
  }

  int n = tas.nrow();
  NumericMatrix result(n, 19);
  colnames(result) = CharacterVector::create(
    "bio01", "bio02", "bio03", "bio04", "bio05",
    "bio06", "bio07", "bio08", "bio09", "bio10",
    "bio11", "bio12", "bio13", "bio14", "bio15",
    "bio16", "bio17", "bio18", "bio19");

  for (int i = 0; i < n; i++) {
    std::array<double, 19> bio = compute_pixel_window(
      tas.row(i), tasmax.row(i), tasmin.row(i), pr.row(i),
      months, window, na_rm);
    for (int j = 0; j < 19; j++) {
      double v = bio[j];
      result(i, j) = std::isnan(v) ? NA_REAL : v;
    }
  }

  return result;
}

//' Compute bioclimatic variables using a rolling window of arbitrary length
//'
//' @param tas    Numeric matrix (pixels x 12): monthly mean temperature.
//' @param tasmax Numeric matrix (pixels x 12): monthly maximum temperature.
//' @param tasmin Numeric matrix (pixels x 12): monthly minimum temperature.
//' @param pr     Numeric matrix (pixels x 12): monthly precipitation.
//' @param window Integer: length (months) of the rolling window (2-11).
//' @param na_rm  Logical: if `TRUE`, skip `NA` months.
//' @return Numeric matrix (pixels x 19) with columns `bio01`..`bio19`.
//'   The base variables (bio01-bio07, bio12-bio15) are computed over the
//'   full 12 months; the rolling-window variables (bio08-bio11, bio16-bio19)
//'   are computed over the best `window`-month period.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix bioclim_rolling_cpp(NumericMatrix tas,
                                  NumericMatrix tasmax,
                                  NumericMatrix tasmin,
                                  NumericMatrix pr,
                                  int window = 3,
                                  bool na_rm = false) {
  if (tas.nrow() != tasmax.nrow() || tas.nrow() != tasmin.nrow() ||
      tas.nrow() != pr.nrow()) {
    stop("All input matrices must have the same number of rows");
  }
  if (tas.ncol() != 12 || tasmax.ncol() != 12 || tasmin.ncol() != 12 ||
      pr.ncol() != 12) {
    stop("All input matrices must have 12 columns");
  }
  if (window < 2 || window > 11) {
    stop("'window' must be between 2 and 11");
  }

  IntegerVector months(12);
  for (int m = 0; m < 12; m++) months[m] = m + 1;

  int n = tas.nrow();
  NumericMatrix result(n, 19);
  colnames(result) = CharacterVector::create(
    "bio01", "bio02", "bio03", "bio04", "bio05",
    "bio06", "bio07", "bio08", "bio09", "bio10",
    "bio11", "bio12", "bio13", "bio14", "bio15",
    "bio16", "bio17", "bio18", "bio19");

  for (int i = 0; i < n; i++) {
    // Base annual statistics.
    NumericVector r_tas    = tas.row(i);
    NumericVector r_tasmax = tasmax.row(i);
    NumericVector r_tasmin = tasmin.row(i);
    NumericVector r_pr     = pr.row(i);

    if (!na_rm && (has_na_sel(r_tas, months) || has_na_sel(r_tasmax, months) ||
                   has_na_sel(r_tasmin, months) || has_na_sel(r_pr, months))) {
      for (int j = 0; j < 19; j++) result(i, j) = NA_REAL;
      continue;
    }

    double b01 = mean_sel(r_tas, months, na_rm);
    double b05 = max_sel(r_tasmax, months, na_rm);
    double b06 = min_sel(r_tasmin, months, na_rm);
    double b12 = sum_sel(r_pr, months, na_rm);
    double b13 = max_sel(r_pr, months, na_rm);
    double b14 = min_sel(r_pr, months, na_rm);
    double pr_mean = mean_sel(r_pr, months, na_rm);

    double b02 = mean_sel(r_tasmax, months, na_rm) -
                 mean_sel(r_tasmin, months, na_rm);
    double b07 = (std::isnan(b05) || std::isnan(b06))
                   ? R_NaN : b05 - b06;
    double b03 = (std::isnan(b07) || b07 == 0.0) ? R_NaN : 100.0 * b02 / b07;
    double b04 = 100.0 * sd_pop_sel(r_tas, months, na_rm);
    double b15 = (std::isnan(pr_mean) || pr_mean <= 0.0)
                   ? R_NaN : 100.0 * sd_pop_sel(r_pr, months, na_rm) / pr_mean;

    result(i, 0)  = b01;
    result(i, 1)  = b02;
    result(i, 2)  = b03;
    result(i, 3)  = b04;
    result(i, 4)  = b05;
    result(i, 5)  = b06;
    result(i, 6)  = b07;
    result(i, 11) = b12;
    result(i, 12) = b13;
    result(i, 13) = b14;
    result(i, 14) = b15;

    // Rolling window (length `window`) over the full year.
    bool mask[12];
    build_month_mask(months, mask);

    int wet_start  = argmax_window_sum(r_pr, months, window, mask, na_rm);
    int dry_start  = argmin_window_sum(r_pr, months, window, mask, na_rm);
    int warm_start = argmax_window_mean(r_tas, months, window, mask, na_rm);
    int cold_start = argmin_window_mean(r_tas, months, window, mask, na_rm);

    if (wet_start > 0) {
      result(i, 7)  = window_mean_sel(r_tas, wet_start, window, na_rm);
      result(i, 15) = window_sum_sel(r_pr, wet_start, window, na_rm);
    } else {
      result(i, 7)  = NA_REAL;
      result(i, 15) = NA_REAL;
    }
    if (dry_start > 0) {
      result(i, 8)  = window_mean_sel(r_tas, dry_start, window, na_rm);
      result(i, 16) = window_sum_sel(r_pr, dry_start, window, na_rm);
    } else {
      result(i, 8)  = NA_REAL;
      result(i, 16) = NA_REAL;
    }
    if (warm_start > 0) {
      result(i, 9)  = window_mean_sel(r_tas, warm_start, window, na_rm);
      result(i, 17) = window_sum_sel(r_pr, warm_start, window, na_rm);
    } else {
      result(i, 9)  = NA_REAL;
      result(i, 17) = NA_REAL;
    }
    if (cold_start > 0) {
      result(i, 10) = window_mean_sel(r_tas, cold_start, window, na_rm);
      result(i, 18) = window_sum_sel(r_pr, cold_start, window, na_rm);
    } else {
      result(i, 10) = NA_REAL;
      result(i, 18) = NA_REAL;
    }
  }

  return result;
}
