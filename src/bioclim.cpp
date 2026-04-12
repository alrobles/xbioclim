#include <Rcpp.h>
using namespace Rcpp;

// ── Primitive helpers ─────────────────────────────────────────────────────────

// Population standard deviation (denominator N, matching xbioclim convention)
static inline double sd_pop_c(const NumericVector& x) {
  int n = x.size();
  double m = mean(x);
  double ss = 0.0;
  for (int i = 0; i < n; i++) {
    double d = x[i] - m;
    ss += d * d;
  }
  return std::sqrt(ss / n);
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
//' @return Numeric matrix (pixels x 19) with one column per variable
//'   (bio01..bio19), named accordingly.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix bioclim_cpp(NumericMatrix tas,
                          NumericMatrix tasmax,
                          NumericMatrix tasmin,
                          NumericMatrix pr) {
  int n = tas.nrow();
  NumericMatrix result(n, 19);

  CharacterVector cnames = CharacterVector::create(
    "bio01", "bio02", "bio03", "bio04", "bio05",
    "bio06", "bio07", "bio08", "bio09", "bio10",
    "bio11", "bio12", "bio13", "bio14", "bio15",
    "bio16", "bio17", "bio18", "bio19"
  );
  colnames(result) = cnames;

  for (int i = 0; i < n; i++) {
    NumericVector t(tas.row(i));
    NumericVector tmx(tasmax.row(i));
    NumericVector tmn(tasmin.row(i));
    NumericVector p(pr.row(i));

    // Temperature basics
    double b01 = mean(t);
    double b02 = mean(tmx - tmn);
    double b05 = max(tmx);
    double b06 = min(tmn);
    double b07 = b05 - b06;
    double b03 = (b07 == 0.0) ? R_NaN : 100.0 * b02 / b07;
    double b04 = 100.0 * sd_pop_c(t);

    // Quarter indices
    int wet_start  = quarter_argmax_c(p);
    int dry_start  = quarter_argmin_c(p);
    int warm_start = quarter_argmax_c(t);
    int cold_start = quarter_argmin_c(t);

    // Temperature quarter means
    double b08 = quarter_mean_c(t, wet_start);
    double b09 = quarter_mean_c(t, dry_start);
    double b10 = quarter_mean_c(t, warm_start);
    double b11 = quarter_mean_c(t, cold_start);

    // Precipitation basics
    double b12 = sum(p);
    double b13 = max(p);
    double b14 = min(p);
    double pr_mean = mean(p);
    double b15 = (pr_mean == 0.0) ? R_NaN : 100.0 * sd_pop_c(p) / pr_mean;

    // Precipitation quarter sums
    double b16 = quarter_sum_c(p, wet_start);
    double b17 = quarter_sum_c(p, dry_start);
    double b18 = quarter_sum_c(p, warm_start);
    double b19 = quarter_sum_c(p, cold_start);

    result(i, 0)  = b01;
    result(i, 1)  = b02;
    result(i, 2)  = b03;
    result(i, 3)  = b04;
    result(i, 4)  = b05;
    result(i, 5)  = b06;
    result(i, 6)  = b07;
    result(i, 7)  = b08;
    result(i, 8)  = b09;
    result(i, 9)  = b10;
    result(i, 10) = b11;
    result(i, 11) = b12;
    result(i, 12) = b13;
    result(i, 13) = b14;
    result(i, 14) = b15;
    result(i, 15) = b16;
    result(i, 16) = b17;
    result(i, 17) = b18;
    result(i, 18) = b19;
  }

  return result;
}
