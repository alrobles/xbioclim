#ifdef _OPENMP
#include <omp.h>
#endif
#include "xbioclim_omp.h"
#include <Rcpp.h>
#include "xbioclim_omp.h"
#include "xbioclim_core/bioclim.hpp"
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

// ── Batch: All 19 variables ───────────────────────────────────────────────────

//' Compute all 19 bioclimatic variables for a raster block
//'
//' @param tas    Numeric matrix (pixels x 12): monthly mean temperature.
//' @param tasmax Numeric matrix (pixels x 12): monthly max temperature.
//' @param tasmin Numeric matrix (pixels x 12): monthly min temperature.
//' @param pr     Numeric matrix (pixels x 12): monthly precipitation.
//' @param ncores Integer: number of OpenMP threads (default 1).
//' @return Numeric matrix (pixels x 19) with one column per variable
//'   (bio01..bio19), named accordingly. Rows with any NA input are returned
//'   as all-NA.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix bioclim_cpp(NumericMatrix tas,
                          NumericMatrix tasmax,
                          NumericMatrix tasmin,
                          NumericMatrix pr,
                          int ncores = 1) {
  int n = tas.nrow();
  NumericMatrix result(n, 19);

  CharacterVector cnames = CharacterVector::create(
    "bio01", "bio02", "bio03", "bio04", "bio05",
    "bio06", "bio07", "bio08", "bio09", "bio10",
    "bio11", "bio12", "bio13", "bio14", "bio15",
    "bio16", "bio17", "bio18", "bio19"
  );
  colnames(result) = cnames;

  xbioclim_core::ClimateBlock data;
  data.tas    = xbioclim_core::Array2D::from_shape({static_cast<std::size_t>(n), std::size_t(12)});
  data.tasmax = xbioclim_core::Array2D::from_shape({static_cast<std::size_t>(n), std::size_t(12)});
  data.tasmin = xbioclim_core::Array2D::from_shape({static_cast<std::size_t>(n), std::size_t(12)});
  data.pr     = xbioclim_core::Array2D::from_shape({static_cast<std::size_t>(n), std::size_t(12)});

  std::vector<int> na_rows(n, 0);
  for (int i = 0; i < n; i++) {
    bool has_na = false;
    for (int m = 0; m < 12; m++) {
      double t   = tas(i, m);
      double tmx = tasmax(i, m);
      double tmn = tasmin(i, m);
      double p   = pr(i, m);
      if (ISNA(t) || ISNA(tmx) || ISNA(tmn) || ISNA(p)) {
        has_na = true;
        data.tas(i, m)    = 0.0f;
        data.tasmax(i, m) = 0.0f;
        data.tasmin(i, m) = 0.0f;
        data.pr(i, m)     = 0.0f;
      } else {
        data.tas(i, m)    = static_cast<xbioclim_core::value_type>(t);
        data.tasmax(i, m) = static_cast<xbioclim_core::value_type>(tmx);
        data.tasmin(i, m) = static_cast<xbioclim_core::value_type>(tmn);
        data.pr(i, m)     = static_cast<xbioclim_core::value_type>(p);
      }
    }
    if (has_na) na_rows[i] = 1;
  }

  xbioclim_core::BioBlock bio = xbioclim_core::compute_bioclim(data);

  for (int i = 0; i < n; i++) {
    if (na_rows[i]) {
      for (int j = 0; j < 19; j++) result(i, j) = NA_REAL;
      continue;
    }
    result(i,  0) = static_cast<double>(bio.bio01(i));
    result(i,  1) = static_cast<double>(bio.bio02(i));
    result(i,  2) = static_cast<double>(bio.bio03(i));
    result(i,  3) = static_cast<double>(bio.bio04(i));
    result(i,  4) = static_cast<double>(bio.bio05(i));
    result(i,  5) = static_cast<double>(bio.bio06(i));
    result(i,  6) = static_cast<double>(bio.bio07(i));
    result(i,  7) = static_cast<double>(bio.bio08(i));
    result(i,  8) = static_cast<double>(bio.bio09(i));
    result(i,  9) = static_cast<double>(bio.bio10(i));
    result(i, 10) = static_cast<double>(bio.bio11(i));
    result(i, 11) = static_cast<double>(bio.bio12(i));
    result(i, 12) = static_cast<double>(bio.bio13(i));
    result(i, 13) = static_cast<double>(bio.bio14(i));
    result(i, 14) = static_cast<double>(bio.bio15(i));
    result(i, 15) = static_cast<double>(bio.bio16(i));
    result(i, 16) = static_cast<double>(bio.bio17(i));
    result(i, 17) = static_cast<double>(bio.bio18(i));
    result(i, 18) = static_cast<double>(bio.bio19(i));
  }

  return result;
}
