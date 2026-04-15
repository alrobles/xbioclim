## R CMD check results

0 errors | 0 warnings | 1 note

The note is:

    * New submission

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
