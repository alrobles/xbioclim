#pragma once

#include "xbioclim/primitives.hpp"

namespace xbioclim {

/// Compute all 19 bioclimatic variables from a block of monthly climate data.
BioBlock compute_bioclim(const ClimateBlock& data);

} // namespace xbioclim
