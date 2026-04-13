// rcpp_bioclim_mod.cpp — Rcpp Modules wrapper for the xbioclim C++ core.
//
// Defines a monolithic Rcpp module "bioclim_mod" that exposes the key
// xbioclim abstractions to R following the terra-style Rcpp Modules pattern:
//
//   ClimateBlock — input data class (four n_pixels × 12 matrices)
//                  Methods: n_pixels(), compute()
//
// Usage from R (after the package is loaded):
//
//   library(rxbioclim)
//   # Single-pixel (1 × 12 matrices)
//   tas    <- matrix(c(5,7,10,14,18,22,25,24,20,15,10,6), nrow=1)
//   tasmax <- matrix(c(8,10,14,18,23,28,32,31,26,19,13,9), nrow=1)
//   tasmin <- matrix(c(1,3,6,10,13,17,20,19,15,10,6,2), nrow=1)
//   pr     <- matrix(c(60,55,50,40,30,15,5,10,25,45,55,65), nrow=1)
//   block  <- new(ClimateBlock, tas, tasmax, tasmin, pr)
//   block$n_pixels()   # 1
//   block$compute()    # 1 × 19 matrix, columns bio01 … bio19

// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include "xbioclim.h"

using namespace Rcpp;

// ---------------------------------------------------------------------------
// Column names helper
// ---------------------------------------------------------------------------

static CharacterVector bio_colnames() {
    CharacterVector nms(19);
    char buf[8];
    for (int i = 0; i < 19; ++i) {
        std::snprintf(buf, sizeof(buf), "bio%02d", i + 1);
        nms[i] = buf;
    }
    return nms;
}

// ---------------------------------------------------------------------------
// ClimateBlock — R-facing wrapper for xbioclim::ClimateBlock
// ---------------------------------------------------------------------------

class ClimateBlock {
public:
    // ---- constructor -------------------------------------------------------

    /// Create a ClimateBlock from four R matrices of shape n_pixels × 12.
    ///
    /// @param tas    Monthly mean temperature matrix    [n_pixels × 12]
    /// @param tasmax Monthly maximum temperature matrix [n_pixels × 12]
    /// @param tasmin Monthly minimum temperature matrix [n_pixels × 12]
    /// @param pr     Monthly precipitation matrix       [n_pixels × 12]
    ClimateBlock(NumericMatrix tas,
                 NumericMatrix tasmax,
                 NumericMatrix tasmin,
                 NumericMatrix pr) {
        const int nrow = tas.nrow();
        const int ncol = tas.ncol();

        if (ncol != 12)
            stop("Each climate matrix must have exactly 12 columns (one per month)");

        if (tasmax.nrow() != nrow || tasmax.ncol() != 12 ||
            tasmin.nrow() != nrow || tasmin.ncol() != 12 ||
            pr.nrow()     != nrow || pr.ncol()     != 12)
            stop("All four matrices must have the same dimensions (n_pixels x 12)");

        const std::size_t N = static_cast<std::size_t>(nrow);

        // Convert R column-major matrices to C++ row-major vectors
        // (pixel × month layout: index = pixel * 12 + month)
        std::vector<double> tas_v(N * 12), tasmax_v(N * 12),
                            tasmin_v(N * 12), pr_v(N * 12);

        for (int p = 0; p < nrow; ++p) {
            for (int m = 0; m < 12; ++m) {
                tas_v[p * 12 + m]    = tas(p, m);
                tasmax_v[p * 12 + m] = tasmax(p, m);
                tasmin_v[p * 12 + m] = tasmin(p, m);
                pr_v[p * 12 + m]     = pr(p, m);
            }
        }

        data_ = xbioclim::ClimateBlock(
            std::move(tas_v), std::move(tasmax_v),
            std::move(tasmin_v), std::move(pr_v), N);
    }

    // ---- methods -----------------------------------------------------------

    /// Return the number of pixels stored in this block.
    int n_pixels() const {
        return static_cast<int>(data_.n_pixels);
    }

    /// Compute all 19 bioclimatic variables for every pixel.
    ///
    /// @return NumericMatrix of shape n_pixels × 19 with column names
    ///         bio01, bio02, …, bio19.
    NumericMatrix compute() const {
        auto bio = xbioclim::compute_bioclim(data_);
        const int N = static_cast<int>(data_.n_pixels);

        NumericMatrix result(N, 19);

        for (int p = 0; p < N; ++p) {
            result(p,  0) = bio.bio01[p];
            result(p,  1) = bio.bio02[p];
            result(p,  2) = bio.bio03[p];
            result(p,  3) = bio.bio04[p];
            result(p,  4) = bio.bio05[p];
            result(p,  5) = bio.bio06[p];
            result(p,  6) = bio.bio07[p];
            result(p,  7) = bio.bio08[p];
            result(p,  8) = bio.bio09[p];
            result(p,  9) = bio.bio10[p];
            result(p, 10) = bio.bio11[p];
            result(p, 11) = bio.bio12[p];
            result(p, 12) = bio.bio13[p];
            result(p, 13) = bio.bio14[p];
            result(p, 14) = bio.bio15[p];
            result(p, 15) = bio.bio16[p];
            result(p, 16) = bio.bio17[p];
            result(p, 17) = bio.bio18[p];
            result(p, 18) = bio.bio19[p];
        }

        colnames(result) = bio_colnames();
        return result;
    }

private:
    xbioclim::ClimateBlock data_;
};

// ---------------------------------------------------------------------------
// RCPP_MODULE definition — "bioclim_mod"
// ---------------------------------------------------------------------------

RCPP_MODULE(bioclim_mod) {
    class_<ClimateBlock>("ClimateBlock")

        .constructor<NumericMatrix, NumericMatrix, NumericMatrix, NumericMatrix>(
            "Create a ClimateBlock from four n_pixels x 12 matrices: "
            "tas (mean temperature), tasmax (maximum temperature), "
            "tasmin (minimum temperature), pr (precipitation).")

        .method("n_pixels", &ClimateBlock::n_pixels,
                "Return the number of pixels stored in this ClimateBlock.")

        .method("compute", &ClimateBlock::compute,
                "Compute all 19 bioclimatic variables (BIO01-BIO19). "
                "Returns a numeric matrix of shape n_pixels x 19 "
                "with columns named bio01 ... bio19.")
    ;
}
