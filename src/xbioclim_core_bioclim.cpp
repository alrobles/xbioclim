#include "xbioclim_core/bioclim.hpp"

#include <xtensor/xoperation.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace xbioclim_core {

static constexpr value_type NODATA = std::numeric_limits<value_type>::quiet_NaN();

BioBlock compute_bioclim(const ClimateBlock& data, bool na_rm) {
    BioBlock bio;

    // --- Validate ClimateBlock invariants ---
    const std::size_t N = data.tas.shape(0);
    if (data.tasmax.shape(0) != N ||
        data.tasmin.shape(0) != N ||
        data.pr.shape(0)     != N) {
        throw std::invalid_argument(
            "ClimateBlock: all fields must have the same n_pixels (shape(0))");
    }
    if (data.tas.shape(1)    != 12 ||
        data.tasmax.shape(1) != 12 ||
        data.tasmin.shape(1) != 12 ||
        data.pr.shape(1)     != 12) {
        throw std::invalid_argument(
            "ClimateBlock: all fields must have exactly 12 months (shape(1))");
    }

    // --- Diurnal range: [N_pixels, 12] ---
    Array2D diurnal = data.tasmax - data.tasmin;

    if (!na_rm) {
        // Classic / conservative path: NaN propagates (wrapper can enforce all-NA).
        bio.bio01 = row_mean(data.tas);
        bio.bio02 = row_mean(diurnal);
        bio.bio05 = row_max(data.tasmax);
        bio.bio06 = row_min(data.tasmin);
        bio.bio07 = bio.bio05 - bio.bio06;

        bio.bio03 = xt::where(xt::not_equal(bio.bio07, value_type(0)),
                              value_type(100) * bio.bio02 / bio.bio07,
                              NODATA);

        bio.bio04 = value_type(100) * row_std(data.tas);
        bio.bio12 = row_sum(data.pr);
        bio.bio13 = row_max(data.pr);
        bio.bio14 = row_min(data.pr);

        Array1D pr_mean = row_mean(data.pr);
        Array1D pr_std  = row_std(data.pr);
        bio.bio15 = xt::where(pr_mean > value_type(0),
                              value_type(100) * pr_std / pr_mean,
                              NODATA);

        IndexArray wet_q  = rolling_quarter_argmax(data.pr);
        IndexArray dry_q  = rolling_quarter_argmin(data.pr);
        IndexArray warm_q = rolling_quarter_argmax(data.tas);
        IndexArray cold_q = rolling_quarter_argmin(data.tas);

        bio.bio08 = quarter_mean(data.tas, wet_q);
        bio.bio09 = quarter_mean(data.tas, dry_q);
        bio.bio10 = quarter_mean(data.tas, warm_q);
        bio.bio11 = quarter_mean(data.tas, cold_q);

        bio.bio16 = quarter_sum(data.pr, wet_q);
        bio.bio17 = quarter_sum(data.pr, dry_q);
        bio.bio18 = quarter_sum(data.pr, warm_q);
        bio.bio19 = quarter_sum(data.pr, cold_q);
    } else {
        // na.rm = TRUE: skip NaNs; a quarter is valid if it has >= 1 valid month.
        bio.bio01 = row_nanmean(data.tas);
        bio.bio02 = row_nanmean(diurnal);
        bio.bio05 = row_nanmax(data.tasmax);
        bio.bio06 = row_nanmin(data.tasmin);
        bio.bio07 = bio.bio05 - bio.bio06;

        bio.bio03 = xt::where(xt::not_equal(bio.bio07, value_type(0)),
                              value_type(100) * bio.bio02 / bio.bio07,
                              NODATA);

        bio.bio04 = value_type(100) * row_nanstd(data.tas);
        bio.bio12 = row_nansum(data.pr);
        bio.bio13 = row_nanmax(data.pr);
        bio.bio14 = row_nanmin(data.pr);

        Array1D pr_mean = row_nanmean(data.pr);
        Array1D pr_std  = row_nanstd(data.pr);
        bio.bio15 = xt::where(pr_mean > value_type(0),
                              value_type(100) * pr_std / pr_mean,
                              NODATA);

        IndexArray wet_q  = nan_rolling_quarter_argmax(data.pr);
        IndexArray dry_q  = nan_rolling_quarter_argmin(data.pr);
        IndexArray warm_q = nan_rolling_quarter_argmax(data.tas);
        IndexArray cold_q = nan_rolling_quarter_argmin(data.tas);

        bio.bio08 = nan_quarter_mean(data.tas, wet_q);
        bio.bio09 = nan_quarter_mean(data.tas, dry_q);
        bio.bio10 = nan_quarter_mean(data.tas, warm_q);
        bio.bio11 = nan_quarter_mean(data.tas, cold_q);

        bio.bio16 = nan_quarter_sum(data.pr, wet_q);
        bio.bio17 = nan_quarter_sum(data.pr, dry_q);
        bio.bio18 = nan_quarter_sum(data.pr, warm_q);
        bio.bio19 = nan_quarter_sum(data.pr, cold_q);
    }

    return bio;
}

} // namespace xbioclim_core
