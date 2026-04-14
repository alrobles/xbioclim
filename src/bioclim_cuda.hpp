// bioclim_cuda.hpp — host-side interface to the CUDA bioclimatic kernel.
//
// This header is included by BioclimEngine.cpp when HAVE_CUDA is defined.
// When the package is built without CUDA all of its contents are invisible.

#pragma once

#ifdef HAVE_CUDA

// Launch the CUDA bioclimatic kernel for a single tile.
//
// All pointer arguments point to host (CPU) memory.  The function handles
// device allocation, host-to-device transfer, kernel launch, and
// device-to-host copy internally.
//
// Data layout:
//   Input:  var[month * n_pix + i]   (12 months, n_pix pixels)
//   Mask:   mask[i]                  (nullptr means no mask)
//   Output: bio[bio_idx * n_pix + i] (19 variables, n_pix pixels)
//
// @param h_tas    Host pointer to mean temperature tile   [12 * n_pix].
// @param h_tasmax Host pointer to max temperature tile    [12 * n_pix].
// @param h_tasmin Host pointer to min temperature tile    [12 * n_pix].
// @param h_pr     Host pointer to precipitation tile      [12 * n_pix].
// @param h_mask   Host pointer to mask tile [n_pix], or nullptr for no mask.
// @param h_bio    Host pointer to output bio tile         [19 * n_pix].
// @param n_pix    Number of pixels in the tile.
void launch_bioclim_cuda(
    const double* h_tas,
    const double* h_tasmax,
    const double* h_tasmin,
    const double* h_pr,
    const double* h_mask,
    double*       h_bio,
    int           n_pix
);

#endif  // HAVE_CUDA
