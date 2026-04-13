#ifdef _OPENMP
#include <omp.h>
#endif
#include <Rcpp.h>
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

#ifdef _OPENMP
  int prev_threads = omp_get_max_threads();
  omp_set_num_threads(ncores);
#endif

  // Use raw REAL() pointers instead of Rcpp proxy objects for OpenMP safety:
  // Rcpp proxy classes are not thread-safe; raw pointers allow concurrent
  // reads from input matrices and non-overlapping writes to the output matrix.
  // Column-major layout: element [i, m] is at ptr[i + m * n].
  const double* tas_ptr    = REAL(tas);
  const double* tasmax_ptr = REAL(tasmax);
  const double* tasmin_ptr = REAL(tasmin);
  const double* pr_ptr     = REAL(pr);
  double*       res_ptr    = REAL(result);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < n; i++) {
    // Copy row i into local stack buffers (private per thread)
    double t[12], tmx[12], tmn[12], p[12];
    for (int m = 0; m < 12; m++) {
      t[m]   = tas_ptr[i   + m * n];
      tmx[m] = tasmax_ptr[i + m * n];
      tmn[m] = tasmin_ptr[i + m * n];
      p[m]   = pr_ptr[i    + m * n];
    }

    // NA check: if any of the 48 input values is NA, write NA row and skip
    bool has_na = false;
    for (int m = 0; m < 12; m++) {
      if (ISNA(t[m]) || ISNA(tmx[m]) || ISNA(tmn[m]) || ISNA(p[m])) {
        has_na = true;
        break;
      }
    }
    if (has_na) {
      for (int j = 0; j < 19; j++) res_ptr[i + j * n] = NA_REAL;
      continue;
    }

    // Temperature basics
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
    double b01     = t_sum / 12.0;
    double b02     = diurnal_sum / 12.0;
    double b05     = tmx_max;
    double b06     = tmn_min;
    double b07     = b05 - b06;
    double b03     = (b07 == 0.0) ? R_NaN : 100.0 * b02 / b07;
    double t_sd    = std::sqrt(t_sumsq / 12.0 - b01 * b01);
    double b04     = 100.0 * t_sd;

    // Quarter indices (temperature)
    int warm_start = quarter_argmax_ptr(t);
    int cold_start = quarter_argmin_ptr(t);

    // Quarter indices (precipitation)
    int wet_start  = quarter_argmax_ptr(p);
    int dry_start  = quarter_argmin_ptr(p);

    // Temperature quarter means
    double b08 = quarter_mean_ptr(t, wet_start);
    double b09 = quarter_mean_ptr(t, dry_start);
    double b10 = quarter_mean_ptr(t, warm_start);
    double b11 = quarter_mean_ptr(t, cold_start);

    // Precipitation basics
    double p_sum = 0.0, p_sumsq = 0.0;
    double p_max = p[0], p_min = p[0];
    for (int m = 0; m < 12; m++) {
      p_sum   += p[m];
      p_sumsq += p[m] * p[m];
      if (p[m] > p_max) p_max = p[m];
      if (p[m] < p_min) p_min = p[m];
    }
    double b12     = p_sum;
    double b13     = p_max;
    double b14     = p_min;
    double pr_mean = p_sum / 12.0;
    double p_sd    = std::sqrt(p_sumsq / 12.0 - pr_mean * pr_mean);
    double b15     = (pr_mean == 0.0) ? R_NaN : 100.0 * p_sd / pr_mean;

    // Precipitation quarter sums
    double b16 = quarter_sum_ptr(p, wet_start);
    double b17 = quarter_sum_ptr(p, dry_start);
    double b18 = quarter_sum_ptr(p, warm_start);
    double b19 = quarter_sum_ptr(p, cold_start);

    // Write results (column-major: result[i, j] == res_ptr[i + j * n])
    res_ptr[i +  0 * n] = b01;
    res_ptr[i +  1 * n] = b02;
    res_ptr[i +  2 * n] = b03;
    res_ptr[i +  3 * n] = b04;
    res_ptr[i +  4 * n] = b05;
    res_ptr[i +  5 * n] = b06;
    res_ptr[i +  6 * n] = b07;
    res_ptr[i +  7 * n] = b08;
    res_ptr[i +  8 * n] = b09;
    res_ptr[i +  9 * n] = b10;
    res_ptr[i + 10 * n] = b11;
    res_ptr[i + 11 * n] = b12;
    res_ptr[i + 12 * n] = b13;
    res_ptr[i + 13 * n] = b14;
    res_ptr[i + 14 * n] = b15;
    res_ptr[i + 15 * n] = b16;
    res_ptr[i + 16 * n] = b17;
    res_ptr[i + 17 * n] = b18;
    res_ptr[i + 18 * n] = b19;
  }

#ifdef _OPENMP
  omp_set_num_threads(prev_threads);
#endif

  return result;
}
