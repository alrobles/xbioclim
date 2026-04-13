#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

using namespace Rcpp;

// ── Internal C++ primitives ──────────────────────────────────────────────────

// Population standard deviation (denominator N, matching xbioclim convention)
static double sd_pop_impl(const std::vector<double>& x) {
  int n = static_cast<int>(x.size());
  if (n == 0) return std::numeric_limits<double>::quiet_NaN();
  double m = std::accumulate(x.begin(), x.end(), 0.0) / n;
  double s = 0.0;
  for (double v : x) s += (v - m) * (v - m);
  return std::sqrt(s / n);
}

// Rolling 3-month quarter sums with circular wrapping
static std::vector<double> rolling_quarter_sum_impl(
    const std::vector<double>& x) {
  std::vector<double> result(12);
  for (int i = 0; i < 12; ++i) {
    result[i] = x[i] + x[(i + 1) % 12] + x[(i + 2) % 12];
  }
  return result;
}

// 1-indexed starting month of the quarter with the maximum rolling sum
static int quarter_argmax_impl(const std::vector<double>& x) {
  std::vector<double> qsums = rolling_quarter_sum_impl(x);
  return static_cast<int>(
             std::max_element(qsums.begin(), qsums.end()) - qsums.begin()) +
         1;
}

// 1-indexed starting month of the quarter with the minimum rolling sum
static int quarter_argmin_impl(const std::vector<double>& x) {
  std::vector<double> qsums = rolling_quarter_sum_impl(x);
  return static_cast<int>(
             std::min_element(qsums.begin(), qsums.end()) - qsums.begin()) +
         1;
}

// Sum of 3 consecutive months starting at 1-indexed `start` (circular)
static double quarter_sum_impl(const std::vector<double>& x, int start) {
  double s = 0.0;
  for (int i = 0; i < 3; ++i) {
    s += x[static_cast<std::size_t>((start - 1 + i) % 12)];
  }
  return s;
}

// Mean of 3 consecutive months starting at 1-indexed `start` (circular)
static double quarter_mean_impl(const std::vector<double>& x, int start) {
  return quarter_sum_impl(x, start) / 3.0;
}

// ── C++ BioclimModel class ───────────────────────────────────────────────────

class BioclimModel {
 public:
  BioclimModel(const std::vector<double>& tas,
               const std::vector<double>& tasmax,
               const std::vector<double>& tasmin,
               const std::vector<double>& pr)
      : tas_(tas), tasmax_(tasmax), tasmin_(tasmin), pr_(pr) {}

  // BIO01: Mean Annual Temperature
  double bio01() const {
    return std::accumulate(tas_.begin(), tas_.end(), 0.0) / 12.0;
  }

  // BIO02: Mean Diurnal Range
  double bio02() const {
    double s = 0.0;
    for (int i = 0; i < 12; ++i) s += (tasmax_[i] - tasmin_[i]);
    return s / 12.0;
  }

  // BIO03: Isothermality (100 * BIO02 / BIO07)
  double bio03() const {
    double b02 = bio02();
    double b07 = bio07();
    if (b07 == 0.0) return std::numeric_limits<double>::quiet_NaN();
    return 100.0 * b02 / b07;
  }

  // BIO04: Temperature Seasonality
  double bio04() const { return 100.0 * sd_pop_impl(tas_); }

  // BIO05: Max Temperature of Warmest Month
  double bio05() const {
    return *std::max_element(tasmax_.begin(), tasmax_.end());
  }

  // BIO06: Min Temperature of Coldest Month
  double bio06() const {
    return *std::min_element(tasmin_.begin(), tasmin_.end());
  }

  // BIO07: Temperature Annual Range (BIO05 - BIO06)
  double bio07() const { return bio05() - bio06(); }

  // BIO08: Mean Temperature of Wettest Quarter
  double bio08() const {
    return quarter_mean_impl(tas_, quarter_argmax_impl(pr_));
  }

  // BIO09: Mean Temperature of Driest Quarter
  double bio09() const {
    return quarter_mean_impl(tas_, quarter_argmin_impl(pr_));
  }

  // BIO10: Mean Temperature of Warmest Quarter
  double bio10() const {
    return quarter_mean_impl(tas_, quarter_argmax_impl(tas_));
  }

  // BIO11: Mean Temperature of Coldest Quarter
  double bio11() const {
    return quarter_mean_impl(tas_, quarter_argmin_impl(tas_));
  }

  // BIO12: Annual Precipitation
  double bio12() const {
    return std::accumulate(pr_.begin(), pr_.end(), 0.0);
  }

  // BIO13: Precipitation of Wettest Month
  double bio13() const {
    return *std::max_element(pr_.begin(), pr_.end());
  }

  // BIO14: Precipitation of Driest Month
  double bio14() const {
    return *std::min_element(pr_.begin(), pr_.end());
  }

  // BIO15: Precipitation Seasonality (CV)
  double bio15() const {
    double pr_mean = std::accumulate(pr_.begin(), pr_.end(), 0.0) / 12.0;
    if (pr_mean == 0.0) return std::numeric_limits<double>::quiet_NaN();
    return 100.0 * sd_pop_impl(pr_) / pr_mean;
  }

  // BIO16: Precipitation of Wettest Quarter
  double bio16() const {
    return quarter_sum_impl(pr_, quarter_argmax_impl(pr_));
  }

  // BIO17: Precipitation of Driest Quarter
  double bio17() const {
    return quarter_sum_impl(pr_, quarter_argmin_impl(pr_));
  }

  // BIO18: Precipitation of Warmest Quarter
  double bio18() const {
    return quarter_sum_impl(pr_, quarter_argmax_impl(tas_));
  }

  // BIO19: Precipitation of Coldest Quarter
  double bio19() const {
    return quarter_sum_impl(pr_, quarter_argmin_impl(tas_));
  }

  // All 19 bioclimatic variables as a named numeric vector
  Rcpp::NumericVector compute_all() const {
    return Rcpp::NumericVector::create(
        Rcpp::Named("bio01") = bio01(), Rcpp::Named("bio02") = bio02(),
        Rcpp::Named("bio03") = bio03(), Rcpp::Named("bio04") = bio04(),
        Rcpp::Named("bio05") = bio05(), Rcpp::Named("bio06") = bio06(),
        Rcpp::Named("bio07") = bio07(), Rcpp::Named("bio08") = bio08(),
        Rcpp::Named("bio09") = bio09(), Rcpp::Named("bio10") = bio10(),
        Rcpp::Named("bio11") = bio11(), Rcpp::Named("bio12") = bio12(),
        Rcpp::Named("bio13") = bio13(), Rcpp::Named("bio14") = bio14(),
        Rcpp::Named("bio15") = bio15(), Rcpp::Named("bio16") = bio16(),
        Rcpp::Named("bio17") = bio17(), Rcpp::Named("bio18") = bio18(),
        Rcpp::Named("bio19") = bio19());
  }

 private:
  std::vector<double> tas_;
  std::vector<double> tasmax_;
  std::vector<double> tasmin_;
  std::vector<double> pr_;
};

// ── Rcpp-exported interface ──────────────────────────────────────────────────

//' Create a new C++ BioclimModel and return an external pointer
//'
//' @param tas    Numeric vector of length 12.
//' @param tasmax Numeric vector of length 12.
//' @param tasmin Numeric vector of length 12.
//' @param pr     Numeric vector of length 12.
//' @return An external pointer wrapping a \code{BioclimModel} C++ object.
//' @keywords internal
// [[Rcpp::export]]
SEXP bioclim_model_new(Rcpp::NumericVector tas, Rcpp::NumericVector tasmax,
                       Rcpp::NumericVector tasmin, Rcpp::NumericVector pr) {
  auto check_len = [](const Rcpp::NumericVector& v, const char* name) {
    if (v.size() != 12)
      Rcpp::stop("'%s' must have length 12 (one value per month), got %d",
                 name, (int)v.size());
  };
  check_len(tas,    "tas");
  check_len(tasmax, "tasmax");
  check_len(tasmin, "tasmin");
  check_len(pr,     "pr");
  if (Rcpp::is_true(Rcpp::any(Rcpp::is_na(tas)))    ||
      Rcpp::is_true(Rcpp::any(Rcpp::is_na(tasmax))) ||
      Rcpp::is_true(Rcpp::any(Rcpp::is_na(tasmin))) ||
      Rcpp::is_true(Rcpp::any(Rcpp::is_na(pr)))) {
    Rcpp::stop("Input vectors must not contain NA values.");
  }
  Rcpp::XPtr<BioclimModel> ptr(
      new BioclimModel(Rcpp::as<std::vector<double>>(tas),
                       Rcpp::as<std::vector<double>>(tasmax),
                       Rcpp::as<std::vector<double>>(tasmin),
                       Rcpp::as<std::vector<double>>(pr)),
      true);
  return ptr;
}

//' Test whether the C++ pointer is null
//' @param ptr An external pointer.
//' @return Logical scalar.
//' @keywords internal
// [[Rcpp::export]]
bool bioclim_model_is_null(SEXP ptr) {
  Rcpp::XPtr<BioclimModel> xptr(ptr);
  return !xptr;
}

// Individual bio-variable accessors -------------------------------------------

//' @keywords internal
// [[Rcpp::export]]
double bioclim_model_bio01(SEXP ptr) {
  return Rcpp::XPtr<BioclimModel>(ptr)->bio01();
}

//' @keywords internal
// [[Rcpp::export]]
double bioclim_model_bio02(SEXP ptr) {
  return Rcpp::XPtr<BioclimModel>(ptr)->bio02();
}

//' @keywords internal
// [[Rcpp::export]]
double bioclim_model_bio03(SEXP ptr) {
  return Rcpp::XPtr<BioclimModel>(ptr)->bio03();
}

//' @keywords internal
// [[Rcpp::export]]
double bioclim_model_bio04(SEXP ptr) {
  return Rcpp::XPtr<BioclimModel>(ptr)->bio04();
}

//' @keywords internal
// [[Rcpp::export]]
double bioclim_model_bio05(SEXP ptr) {
  return Rcpp::XPtr<BioclimModel>(ptr)->bio05();
}

//' @keywords internal
// [[Rcpp::export]]
double bioclim_model_bio06(SEXP ptr) {
  return Rcpp::XPtr<BioclimModel>(ptr)->bio06();
}

//' @keywords internal
// [[Rcpp::export]]
double bioclim_model_bio07(SEXP ptr) {
  return Rcpp::XPtr<BioclimModel>(ptr)->bio07();
}

//' @keywords internal
// [[Rcpp::export]]
double bioclim_model_bio08(SEXP ptr) {
  return Rcpp::XPtr<BioclimModel>(ptr)->bio08();
}

//' @keywords internal
// [[Rcpp::export]]
double bioclim_model_bio09(SEXP ptr) {
  return Rcpp::XPtr<BioclimModel>(ptr)->bio09();
}

//' @keywords internal
// [[Rcpp::export]]
double bioclim_model_bio10(SEXP ptr) {
  return Rcpp::XPtr<BioclimModel>(ptr)->bio10();
}

//' @keywords internal
// [[Rcpp::export]]
double bioclim_model_bio11(SEXP ptr) {
  return Rcpp::XPtr<BioclimModel>(ptr)->bio11();
}

//' @keywords internal
// [[Rcpp::export]]
double bioclim_model_bio12(SEXP ptr) {
  return Rcpp::XPtr<BioclimModel>(ptr)->bio12();
}

//' @keywords internal
// [[Rcpp::export]]
double bioclim_model_bio13(SEXP ptr) {
  return Rcpp::XPtr<BioclimModel>(ptr)->bio13();
}

//' @keywords internal
// [[Rcpp::export]]
double bioclim_model_bio14(SEXP ptr) {
  return Rcpp::XPtr<BioclimModel>(ptr)->bio14();
}

//' @keywords internal
// [[Rcpp::export]]
double bioclim_model_bio15(SEXP ptr) {
  return Rcpp::XPtr<BioclimModel>(ptr)->bio15();
}

//' @keywords internal
// [[Rcpp::export]]
double bioclim_model_bio16(SEXP ptr) {
  return Rcpp::XPtr<BioclimModel>(ptr)->bio16();
}

//' @keywords internal
// [[Rcpp::export]]
double bioclim_model_bio17(SEXP ptr) {
  return Rcpp::XPtr<BioclimModel>(ptr)->bio17();
}

//' @keywords internal
// [[Rcpp::export]]
double bioclim_model_bio18(SEXP ptr) {
  return Rcpp::XPtr<BioclimModel>(ptr)->bio18();
}

//' @keywords internal
// [[Rcpp::export]]
double bioclim_model_bio19(SEXP ptr) {
  return Rcpp::XPtr<BioclimModel>(ptr)->bio19();
}

//' Compute all 19 bioclimatic variables from the C++ object
//' @param ptr An external pointer to a \code{BioclimModel} C++ object.
//' @return Named numeric vector of length 19.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector bioclim_model_compute(SEXP ptr) {
  return Rcpp::XPtr<BioclimModel>(ptr)->compute_all();
}
