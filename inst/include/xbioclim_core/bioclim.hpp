#pragma once

#include "xbioclim_core/primitives.hpp"

namespace xbioclim_core {

/// Compute all 19 bioclimatic variables from a block of monthly climate data.
BioBlock compute_bioclim(const ClimateBlock& data);

} // namespace xbioclim_core
