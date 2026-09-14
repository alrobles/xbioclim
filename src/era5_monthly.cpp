// era5_monthly.cpp — Rcpp exports for ERA5-Land hourly-to-monthly aggregation.
//
// Provides two exported functions:
//   era5_t2m_to_monthly_cpp()  — hourly t2m → monthly tas, tasmax, tasmin
//   era5_tp_to_monthly_cpp()   — hourly tp  → monthly pr
//
// Both accept R NumericMatrix (n_pixels × n_hours, column-major) and return
// results using zero-copy writes into pre-allocated R memory.

#ifdef _OPENMP
#include <omp.h>
#endif
#include "xbioclim_omp.h"
#include <Rcpp.h>
#include "era5_monthly.h"
using namespace Rcpp;

//' Aggregate ERA5-Land hourly 2-m temperature to monthly statistics
//'
//' Converts an hourly temperature matrix to monthly mean temperature (tas),
//' monthly mean of daily maxima (tasmax), and monthly mean of daily minima
//' (tasmin), following the CHELSA variable convention.
//'
//' @param hourly Numeric matrix (n_pixels x n_hours): hourly 2-m temperature.
//'   Column-major layout.  n_hours must equal 24 * n_days.
//' @param n_days Integer: number of days in the month.
//' @param to_celsius Logical: if TRUE, convert Kelvin to Celsius (default
//'   FALSE, output in Kelvin matching CHELSA convention).
//' @param ncores Integer: number of OpenMP threads (default 1).
//' @return A named list with three numeric vectors of length n_pixels:
//'   \code{tas}, \code{tasmax}, \code{tasmin}.
//' @keywords internal
// [[Rcpp::export]]
List era5_t2m_to_monthly_cpp(NumericMatrix hourly,
                             int n_days,
                             bool to_celsius = false,
                             int ncores = 1) {
    const int n_pixels = hourly.nrow();
    const int n_hours  = hourly.ncol();

    if (n_hours != 24 * n_days) {
        stop("n_hours (%d) must equal 24 * n_days (%d)", n_hours, 24 * n_days);
    }

    NumericVector tas(n_pixels);
    NumericVector tasmax(n_pixels);
    NumericVector tasmin(n_pixels);

    const double* h_p   = REAL(hourly);
    double* tas_p    = REAL(tas);
    double* tasmax_p = REAL(tasmax);
    double* tasmin_p = REAL(tasmin);

    const double offset = to_celsius ? xbioclim::era5::K_TO_C : 0.0;

#ifdef _OPENMP
    int prev_threads = omp_get_max_threads();
    omp_set_num_threads(xbioclim_safe_threads(ncores));
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < n_pixels; ++i) {
        double sum_dmean = 0.0;
        double sum_dmax  = 0.0;
        double sum_dmin  = 0.0;

        for (int d = 0; d < n_days; ++d) {
            const int h_start = d * 24;
            double dsum = 0.0;
            double dmax = h_p[i + h_start * n_pixels];
            double dmin = dmax;
            int count = 0;

            for (int h = 0; h < 24; ++h) {
                double v = h_p[i + (h_start + h) * n_pixels];
                if (ISNA(v)) continue;
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

        tas_p[i]    = (sum_dmean / static_cast<double>(n_days)) + offset;
        tasmax_p[i] = (sum_dmax  / static_cast<double>(n_days)) + offset;
        tasmin_p[i] = (sum_dmin  / static_cast<double>(n_days)) + offset;
    }

#ifdef _OPENMP
    omp_set_num_threads(prev_threads);
#endif

    return List::create(
        Named("tas")    = tas,
        Named("tasmax") = tasmax,
        Named("tasmin") = tasmin
    );
}

//' Aggregate ERA5-Land hourly total precipitation to monthly total
//'
//' Sums hourly precipitation accumulations and converts from metres of water
//' to kg m-2 month-1 (equivalent to mm/month), matching the CHELSA \code{pr}
//' variable convention.
//'
//' @param hourly Numeric matrix (n_pixels x n_hours): hourly total
//'   precipitation in metres.
//' @param ncores Integer: number of OpenMP threads (default 1).
//' @return Numeric vector of length n_pixels: monthly total precipitation
//'   in kg m-2 (mm).
//' @keywords internal
// [[Rcpp::export]]
NumericVector era5_tp_to_monthly_cpp(NumericMatrix hourly,
                                     int ncores = 1) {
    const int n_pixels = hourly.nrow();
    const int n_hours  = hourly.ncol();

    NumericVector pr(n_pixels);
    const double* h_p = REAL(hourly);
    double* pr_p      = REAL(pr);

#ifdef _OPENMP
    int prev_threads = omp_get_max_threads();
    omp_set_num_threads(xbioclim_safe_threads(ncores));
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < n_pixels; ++i) {
        double total = 0.0;
        for (int h = 0; h < n_hours; ++h) {
            double v = h_p[i + h * n_pixels];
            if (!ISNA(v) && v > 0.0) {
                total += v;
            }
        }
        pr_p[i] = total * xbioclim::era5::M_TO_MM;
    }

#ifdef _OPENMP
    omp_set_num_threads(prev_threads);
#endif

    return pr;
}

//' Unified ERA5-Land hourly-to-monthly aggregation
//'
//' Converts hourly 2-m temperature and total precipitation to the four
//' CHELSA-compatible monthly climate variables in a single parallel pass.
//'
//' @param hourly_t2m Numeric matrix (n_pixels x n_hours_t2m): hourly 2-m
//'   temperature in Kelvin.  n_hours_t2m must equal 24 * n_days.
//' @param hourly_tp Numeric matrix (n_pixels x n_hours_tp): hourly total
//'   precipitation in metres.  n_pixels must match \code{hourly_t2m}.
//' @param n_days Integer: number of days in the month.
//' @param to_celsius Logical: convert temperatures from Kelvin to Celsius?
//'   Default FALSE.
//' @param ncores Integer: number of OpenMP threads (default 1).
//' @return A named list with four numeric vectors of length n_pixels:
//'   \code{tas}, \code{tasmax}, \code{tasmin}, \code{pr}.
//' @keywords internal
// [[Rcpp::export]]
List era5_to_monthly_cpp(NumericMatrix hourly_t2m,
                         NumericMatrix hourly_tp,
                         int n_days,
                         bool to_celsius = false,
                         int ncores = 1) {
    const int n_pixels   = hourly_t2m.nrow();
    const int n_hours_t  = hourly_t2m.ncol();
    const int n_hours_p  = hourly_tp.ncol();

    if (hourly_tp.nrow() != n_pixels) {
        stop("hourly_t2m (%d rows) and hourly_tp (%d rows) must have the same number of pixels",
             n_pixels, hourly_tp.nrow());
    }
    if (n_hours_t != 24 * n_days) {
        stop("hourly_t2m has %d columns but expected 24 * %d = %d",
             n_hours_t, n_days, 24 * n_days);
    }

    NumericVector tas(n_pixels);
    NumericVector tasmax(n_pixels);
    NumericVector tasmin(n_pixels);
    NumericVector pr(n_pixels);

    const double* t_p  = REAL(hourly_t2m);
    const double* p_p  = REAL(hourly_tp);
    double* tas_p      = REAL(tas);
    double* tasmax_p   = REAL(tasmax);
    double* tasmin_p   = REAL(tasmin);
    double* pr_p       = REAL(pr);

    const double offset = to_celsius ? xbioclim::era5::K_TO_C : 0.0;

#ifdef _OPENMP
    int prev_threads = omp_get_max_threads();
    omp_set_num_threads(xbioclim_safe_threads(ncores));
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < n_pixels; ++i) {
        // Temperature aggregation
        double sum_dmean = 0.0;
        double sum_dmax  = 0.0;
        double sum_dmin  = 0.0;

        for (int d = 0; d < n_days; ++d) {
            const int h_start = d * 24;
            double dsum = 0.0;
            double dmax = t_p[i + h_start * n_pixels];
            double dmin = dmax;
            int count = 0;

            for (int h = 0; h < 24; ++h) {
                double v = t_p[i + (h_start + h) * n_pixels];
                if (ISNA(v)) continue;
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

        tas_p[i]    = (sum_dmean / static_cast<double>(n_days)) + offset;
        tasmax_p[i] = (sum_dmax  / static_cast<double>(n_days)) + offset;
        tasmin_p[i] = (sum_dmin  / static_cast<double>(n_days)) + offset;

        // Precipitation aggregation
        double total = 0.0;
        for (int h = 0; h < n_hours_p; ++h) {
            double v = p_p[i + h * n_pixels];
            if (!ISNA(v) && v > 0.0) {
                total += v;
            }
        }
        pr_p[i] = total * xbioclim::era5::M_TO_MM;
    }

#ifdef _OPENMP
    omp_set_num_threads(prev_threads);
#endif

    return List::create(
        Named("tas")    = tas,
        Named("tasmax") = tasmax,
        Named("tasmin") = tasmin,
        Named("pr")     = pr
    );
}
