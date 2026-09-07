#pragma once

#include "xbioclim_core/primitives.hpp"

namespace xbioclim_core {

/// Compute all 19 bioclimatic variables from a block of monthly climate data.
/// @param data  input climate block
/// @param na_rm if true, treat NaN as missing and calculate each BIO from the
///              available months (quarter may be partial, min. 1 valid month).
///              if false, a single NaN in any input for a pixel yields all NaNs.
BioBlock compute_bioclim(const ClimateBlock& data, bool na_rm = false);

} // namespace xbioclim_core
