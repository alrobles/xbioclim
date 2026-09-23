## R CMD check results

0 errors | 0 warnings | 1 note

The note is:

    * New submission

## Resubmission notes (1.0.3)

Version 1.0.2 received reviewer feedback requesting `\value` tags in .Rd
files for exported methods. `\value` sections were added to all
documentation pages that lacked them: `BioclimData-class`,
`BioclimModel-class`, `BioclimModel-methods`, `ClimateBlock`,
`bioclim-raster`, `primitives`, `quarterly_fixed`, `quarterly_rolling`,
`xbioclim-messages`, and `xbioclim-package`.

The reviewer also asked for references describing the methods in the
package. The computational method implemented here is new and has no
published reference with a DOI yet, so the 'xbioclim' C++ library that
implements it is cited in the Description field via URL:
Robles Fernandez (2026) <https://github.com/alrobles/xbioclimcpp>.

## Resubmission notes (1.0.2)

Version 1.0.1 was stopped at prescreen with two additional NOTEs:

* `Compilation used the following non-portable flag(s): '-march=native'` —
  removed from `src/Makevars.in` (now `-O2` only).
* `Possibly misspelled words in DESCRIPTION` (`Bioclimatic`, `WorldClim`,
  `bioclimatic`) — these are technical terms now covered by the package
  `.aspell` dictionary shipped in the tarball.

## Test environments

* **Linux** (Ubuntu 22.04, x86_64), R 4.4.0 — CPU-only build (no GDAL, no CUDA)
* **Linux** (Ubuntu 22.04, x86_64), R 4.4.0 — GDAL 3.8 build, CPU-only
* **Linux** (Ubuntu 22.04, x86_64), R devel — CPU-only build
* **macOS** (macOS 14 Sonoma, arm64), R 4.4.0 — CPU-only build
* **Windows** (Windows Server 2022, x86_64), R 4.4.0 — CPU-only build
  (checked via `devtools::check_win_devel()`)
* **win-builder** (R-devel) — CPU-only build

GPU environments (informational; not required for CRAN acceptance):

* **Linux** (Ubuntu 22.04, x86_64), R 4.4.0 — GDAL 3.8 + CUDA 12.2 build,
  tested on an NVIDIA A100 GPU

## Optional system dependencies

Two optional compile-time features are detected automatically by `configure.ac`
and have **no impact** on package functionality when absent:

* **GDAL ≥ 2.0.1** — enables the native tiled I/O pipeline (`BioclimEngine`).
  When `gdal-config` is not found the package compiles without GDAL; the two
  diagnostic functions (`gdal_can_open()`, `gdal_info()`) stop with an
  informative error message. All GDAL-dependent tests are skipped automatically
  via `skip_without_gdal()`.

* **CUDA Toolkit ≥ 11.0** (`nvcc`) — enables GPU-accelerated computation.
  When `nvcc` is not found the package compiles in CPU-only mode; `has_cuda()`
  returns `FALSE` and `cuda_info()` returns an empty list. All CUDA-dependent
  tests are skipped automatically via `skip_if(!has_cuda())`.

Both features can be disabled explicitly at configure time:
`--without-gdal` / `--without-cuda`.

## Downstream dependencies

There are currently **no downstream dependencies** on CRAN or Bioconductor.
