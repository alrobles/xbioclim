#pragma once

#ifdef XBIOCLIM_USE_CUDA

#include "xbioclim/primitives.hpp"

namespace xbioclim {

/// Compute all 19 bioclimatic variables using a CUDA GPU kernel.
///
/// Each GPU thread handles one pixel independently. Input data is transferred
/// to device memory, the kernel is launched, and results are copied back.
///
/// Requires CUDA-capable GPU (target architecture ≥ SM_80).
/// Falls back gracefully when compiled without CUDA support.
BioBlock compute_bioclim_cuda(const ClimateBlock& data);

} // namespace xbioclim

#endif // XBIOCLIM_USE_CUDA
