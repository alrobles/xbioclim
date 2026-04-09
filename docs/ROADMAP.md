# xbioclim — Project Roadmap

> **Goal:** Translate the R/Rcpp proof-of-concept ([alrobles/fastBioClim](https://github.com/alrobles/fastBioClim)) into a high-performance, pure C++17 library that computes the 19 CHELSA/WorldClim bioclimatic variables (BIO01–BIO19) from monthly climate rasters using GDAL for raster I/O and xtensor for array algebra, targeting HPC clusters with GPU hardware.

---

## Table of Contents

1. [Phase 0 — Project Setup](#phase-0--project-setup)
2. [Phase 1 — Test Data & Validation Framework](#phase-1--test-data--validation-framework)
3. [Phase 2 — Core Library (xtensor Primitives)](#phase-2--core-library-xtensor-primitives)
4. [Phase 3 — GDAL I/O Layer](#phase-3--gdal-io-layer)
5. [Phase 4 — BIO Variable Computation](#phase-4--bio-variable-computation)
6. [Phase 5 — CLI Tool & Integration](#phase-5--cli-tool--integration)
7. [Phase 6 — GPU Acceleration](#phase-6--gpu-acceleration)
8. [Phase 7 — Performance Benchmarking & Release](#phase-7--performance-benchmarking--release)

---

## Phase 0 — Project Setup

### 0.1 Repository Directory Layout

```
xbioclim/
├── CMakeLists.txt               # Top-level build configuration
├── cmake/                       # CMake helper modules
│   ├── FindGDAL.cmake
│   └── FindXtensor.cmake
├── include/
│   └── xbioclim/
│       ├── primitives.hpp       # Core xtensor primitives declarations
│       ├── bioclim.hpp          # BIO01–BIO19 computation API
│       ├── gdal_io.hpp          # GDAL raster reader/writer
│       └── version.hpp          # Library version macros
├── src/
│   ├── primitives.cpp           # xtensor primitive implementations
│   ├── bioclim.cpp              # BIO variable implementations
│   └── gdal_io.cpp              # GDAL block I/O implementation
├── tests/
│   ├── CMakeLists.txt
│   ├── data/                    # Synthetic GeoTIFF test fixtures (generated)
│   ├── test_primitives.cpp      # Unit tests for xtensor primitives
│   ├── test_bioclim.cpp         # Integration tests for all 19 BIO variables
│   └── test_gdal_io.cpp         # GDAL reader/writer roundtrip tests
├── tools/
│   └── generate_test_data.py    # Script to create synthetic 10×10 GeoTIFFs
├── scripts/
│   └── run_slurm_array.sh       # Slurm batch submission script
├── docs/
│   ├── ROADMAP.md               # This file
│   └── PROTOCOL.md              # Detailed algorithm & data-layout protocol
├── .github/
│   └── workflows/
│       ├── ci.yml               # Build + test on Ubuntu
│       └── release.yml          # Tag-triggered release pipeline
├── LICENSE
└── README.md
```

### 0.2 CMakeLists.txt Skeleton

The top-level `CMakeLists.txt` must:

- Require **CMake ≥ 3.15** and **C++17**.
- Find **GDAL ≥ 3.0** (mandatory).
- Find **xtensor** (mandatory; optionally via vcpkg or Conan).
- Find **OpenMP** (optional; enables CPU multi-threading).
- Find **CUDA** (optional; enables GPU kernels).
- Export an `xbioclim::xbioclim` CMake target for downstream consumers.

Key CMake options:

| Option | Default | Description |
|---|---|---|
| `XBIOCLIM_BUILD_TESTS` | `ON` | Build the test suite |
| `XBIOCLIM_USE_OPENMP` | `ON` | Enable OpenMP CPU parallelism |
| `XBIOCLIM_USE_CUDA` | `OFF` | Enable CUDA GPU kernels |
| `XBIOCLIM_USE_OPENMP_OFFLOAD` | `OFF` | Enable OpenMP 4.5 GPU offload |

Minimum `CMakeLists.txt` excerpt:

```cmake
cmake_minimum_required(VERSION 3.15)
project(xbioclim VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Dependencies ---
find_package(GDAL 3.0 REQUIRED)
find_package(xtensor REQUIRED)

option(XBIOCLIM_USE_OPENMP "Enable OpenMP" ON)
if(XBIOCLIM_USE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

# --- Library target ---
add_library(xbioclim
  src/primitives.cpp
  src/bioclim.cpp
  src/gdal_io.cpp
)
target_include_directories(xbioclim PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
target_link_libraries(xbioclim PUBLIC
  GDAL::GDAL
  xtensor
  $<$<BOOL:${XBIOCLIM_USE_OPENMP}>:OpenMP::OpenMP_CXX>
)

# --- Tests ---
option(XBIOCLIM_BUILD_TESTS "Build tests" ON)
if(XBIOCLIM_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
```

### 0.3 Dependency Management

#### Option A — System Packages (Ubuntu/Debian, recommended for CI)

```bash
sudo apt-get install -y \
  libgdal-dev \
  libxtensor-dev \
  catch2              # or libgtest-dev for GoogleTest
```

#### Option B — vcpkg

```bash
vcpkg install gdal xtensor catch2
cmake -DCMAKE_TOOLCHAIN_FILE=$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake ..
```

#### Option C — Conan

```ini
# conanfile.txt
[requires]
gdal/3.6.2
xtensor/0.24.7
catch2/3.4.0

[generators]
CMakeDeps
CMakeToolchain
```

```bash
conan install . --output-folder=build --build=missing
cmake -B build -DCMAKE_TOOLCHAIN_FILE=build/conan_toolchain.cmake
```

### 0.4 CI/CD Skeleton (GitHub Actions)

#### `.github/workflows/ci.yml`

```yaml
name: CI

on:
  push:
    branches: [main, develop]
  pull_request:
    branches: [main]

jobs:
  build-and-test:
    runs-on: ubuntu-22.04
    steps:
      - uses: actions/checkout@v4

      - name: Install system dependencies
        run: |
          sudo apt-get update
          sudo apt-get install -y \
            build-essential cmake \
            libgdal-dev gdal-bin \
            libxtensor-dev \
            catch2 python3-gdal python3-numpy

      - name: Generate test data
        run: python3 tools/generate_test_data.py --outdir tests/data

      - name: Configure
        run: cmake -B build -DCMAKE_BUILD_TYPE=Release -DXBIOCLIM_BUILD_TESTS=ON

      - name: Build
        run: cmake --build build --parallel $(nproc)

      - name: Test
        run: ctest --test-dir build --output-on-failure

  static-analysis:
    runs-on: ubuntu-22.04
    steps:
      - uses: actions/checkout@v4
      - name: Run clang-tidy
        run: |
          sudo apt-get install -y clang-tidy libgdal-dev libxtensor-dev
          cmake -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
          clang-tidy src/*.cpp include/xbioclim/*.hpp \
            -p build --warnings-as-errors='*'
```

---

## Phase 1 — Test Data & Validation Framework

### 1.1 Synthetic GeoTIFF Specification

Generate **48 compressed 10×10 integer GeoTIFFs** (4 variables × 12 months) using `tools/generate_test_data.py` (Python + GDAL/rasterio).

All pixels within a layer share the **same uniform integer value** — this makes analytical calculation of all 19 BIO values straightforward.

| Variable | Month *m* (1–12) | Layer value | Unit (CHELSA raw) |
|---|---|---|---|
| `tas` | *m* | *m* | K |
| `tasmax` | *m* | *m* + 1 | K |
| `tasmin` | *m* | *m* − 1 | K |
| `pr` | *m* | *m* | kg m⁻² month⁻¹ |

So the 12 per-month values for each variable are:

| Month | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `tas` | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
| `tasmax` | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 |
| `tasmin` | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
| `pr` | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |

Output files follow the naming convention `<var>_<MM>.tif`, e.g. `tas_01.tif`, `pr_12.tif`.

#### File metadata:
- **Size:** 10 × 10 pixels
- **Bands:** 1
- **Data type:** Int16 (GDT_Int16)
- **Compression:** LZW (COMPRESS=LZW)
- **CRS:** EPSG:4326
- **Geotransform:** `(-180.0, 36.0, 0.0, 90.0, 0.0, -18.0)` (10° cells)
- **NoData:** -9999

#### `tools/generate_test_data.py` outline:

```python
#!/usr/bin/env python3
"""Generate synthetic 10×10 GeoTIFF test fixtures for xbioclim."""

import argparse, pathlib, numpy as np
from osgeo import gdal, osr

VARS = {
    "tas":    lambda m: m,
    "tasmax": lambda m: m + 1,
    "tasmin": lambda m: m - 1,
    "pr":     lambda m: m,
}
GT = (-180.0, 36.0, 0.0, 90.0, 0.0, -18.0)

def write_tif(path, value, nodata=-9999):
    drv = gdal.GetDriverByName("GTiff")
    ds = drv.Create(str(path), 10, 10, 1, gdal.GDT_Int16,
                    ["COMPRESS=LZW"])
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    ds.SetProjection(srs.ExportToWkt())
    ds.SetGeoTransform(GT)
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    band.WriteArray(np.full((10, 10), value, dtype=np.int16))
    ds.FlushCache()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="tests/data")
    args = ap.parse_args()
    out = pathlib.Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    for var, fn in VARS.items():
        for m in range(1, 13):
            write_tif(out / f"{var}_{m:02d}.tif", fn(m))
    print(f"Generated {4 * 12} GeoTIFFs in {out}/")

if __name__ == "__main__":
    main()
```

### 1.2 Test Framework

Use **Catch2 v3** (header-only option available; integrates cleanly with CMake `FetchContent`).

```cmake
# tests/CMakeLists.txt
find_package(Catch2 3 REQUIRED)

add_executable(test_xbioclim
  test_primitives.cpp
  test_bioclim.cpp
  test_gdal_io.cpp
)
target_link_libraries(test_xbioclim PRIVATE xbioclim Catch2::Catch2WithMain)

include(Catch)
catch_discover_tests(test_xbioclim)
```

Alternatively, **GoogleTest** is an equally valid choice — the main criterion is that the framework is available as a system package in the CI environment (`libgtest-dev` on Ubuntu).

### 1.3 Expected Analytical Values for All 19 BIO Variables

Given the mock data defined in §1.1, all pixels are spatially uniform so each BIO variable reduces to a scalar computation. The derivations below use 0-based month indices for array notation; the rolling quarter windows wrap circularly (month 12 is followed by month 1).

**Intermediate quantities:**

```
tas    = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]   (months 1–12)
tasmax = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
tasmin = [0, 1, 2, 3, 4, 5, 6, 7,  8,  9, 10, 11]
pr     = [1, 2, 3, 4, 5, 6, 7, 8,  9, 10, 11, 12]

mean(tas)  = 6.5
mean(pr)   = 6.5
diurnal    = tasmax - tasmin = [2, 2, 2, ..., 2]    (all 2)
mean(diurnal) = 2.0

population variance(tas):
  Σ(xᵢ - 6.5)² / 12 = 2×(0.5²+1.5²+2.5²+3.5²+4.5²+5.5²) / 12
                     = 2×(0.25+2.25+6.25+12.25+20.25+30.25) / 12
                     = 2×71.5 / 12 = 143/12 ≈ 11.9167
sd(tas) = √(143/12) ≈ 3.4521   [population SD; same for pr since data identical]

Rolling 3-month pr sums (starting month index 1..12, circular):
  start=1:  1+ 2+ 3 =  6
  start=2:  2+ 3+ 4 =  9
  start=3:  3+ 4+ 5 = 12
  start=4:  4+ 5+ 6 = 15
  start=5:  5+ 6+ 7 = 18
  start=6:  6+ 7+ 8 = 21
  start=7:  7+ 8+ 9 = 24
  start=8:  8+ 9+10 = 27
  start=9:  9+10+11 = 30
  start=10: 10+11+12 = 33   ← wettest
  start=11: 11+12+ 1 = 24
  start=12: 12+ 1+ 2 = 15

Rolling 3-month tas means (starting month index 1..12, circular):
  start=1:  (1+ 2+ 3)/3  = 2.0   ← coldest
  start=2:  (2+ 3+ 4)/3  = 3.0
  ...
  start=10: (10+11+12)/3 = 11.0  ← warmest
  start=11: (11+12+ 1)/3 = 8.0
  start=12: (12+ 1+ 2)/3 = 5.0
```

**Complete BIO variable expected values:**

| Variable | Formula | Expected Value |
|---|---|---|
| **BIO01** | mean(tas) | **6.5** |
| **BIO02** | mean(tasmax − tasmin) | **2.0** |
| **BIO03** | 100 × BIO02 / BIO07 | **≈ 15.385** (= 200/13) |
| **BIO04** | 100 × sd(tas) [population] | **≈ 345.21** |
| **BIO05** | max(tasmax) | **13.0** |
| **BIO06** | min(tasmin) | **0.0** |
| **BIO07** | BIO05 − BIO06 | **13.0** |
| **BIO08** | mean tas, wettest quarter (start=10, months 10,11,12) | **11.0** |
| **BIO09** | mean tas, driest quarter (start=1, months 1,2,3) | **2.0** |
| **BIO10** | mean tas, warmest quarter (start=10) | **11.0** |
| **BIO11** | mean tas, coldest quarter (start=1) | **2.0** |
| **BIO12** | sum(pr) | **78.0** |
| **BIO13** | max(pr) | **12.0** |
| **BIO14** | min(pr) | **1.0** |
| **BIO15** | 100 × sd(pr) / mean(pr) | **≈ 53.11** (= 100 × 3.4521/6.5) |
| **BIO16** | mean pr, wettest quarter (start=10, months 10,11,12) | **11.0** |
| **BIO17** | mean pr, driest quarter (start=1, months 1,2,3) | **2.0** |
| **BIO18** | mean pr, warmest quarter by tas (start=10) | **11.0** |
| **BIO19** | mean pr, coldest quarter by tas (start=1) | **2.0** |

> **Notes on BIO03 and BIO15:**
> - BIO03 = 100 × (2/13) = 200/13 ≈ 15.3846; test with tolerance ε = 1×10⁻³.
> - BIO15 = 100 × √(143/12) / 6.5; test with tolerance ε = 1×10⁻³.
> - BIO04 = 100 × √(143/12) ≈ 345.21; test with tolerance ε = 1×10⁻³.

#### Sample Catch2 test structure:

```cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "xbioclim/bioclim.hpp"

using Catch::Matchers::WithinAbs;

TEST_CASE("BIO01 - Mean Annual Temperature", "[bioclim]") {
    // tas = [1..12], all pixels uniform
    auto bio = compute_bioclim(load_test_data("tests/data"));
    CHECK_THAT(bio.bio01[0], WithinAbs(6.5, 1e-3));
}

TEST_CASE("BIO07 - Annual Temperature Range", "[bioclim]") {
    auto bio = compute_bioclim(load_test_data("tests/data"));
    CHECK_THAT(bio.bio07[0], WithinAbs(13.0, 1e-3));
}

// ... (one section per BIO variable)
```

---

## Phase 2 — Core Library (xtensor Primitives)

### 2.1 Data Layout

All input data is loaded into a 2-D xtensor array of shape **[N_pixels, 12]**, where `N_pixels = n_rows × n_cols` of the raster tile.

```cpp
using Array2D = xt::xtensor<float, 2>;  // shape: [N_pixels, 12]
using Array1D = xt::xtensor<float, 1>;  // shape: [N_pixels]
```

A tile-processing struct holds the four input arrays:

```cpp
struct ClimateBlock {
    Array2D tas;    // [N_pixels, 12]  monthly mean temperature
    Array2D tasmax; // [N_pixels, 12]  monthly max temperature
    Array2D tasmin; // [N_pixels, 12]  monthly min temperature
    Array2D pr;     // [N_pixels, 12]  monthly precipitation
};

struct BioBlock {
    Array1D bio01, bio02, bio03, bio04, bio05,
            bio06, bio07, bio08, bio09, bio10,
            bio11, bio12, bio13, bio14, bio15,
            bio16, bio17, bio18, bio19;  // each [N_pixels]
};
```

### 2.2 Rcpp → xtensor Primitive Mapping

| Rcpp Primitive | xtensor Equivalent | Notes |
|---|---|---|
| `parallel_average(mat)` | `xt::mean(mat, 1)` | Row-wise mean over 12 months; axis=1 |
| `parallel_sum(mat)` | `xt::sum(mat, 1)` | Row-wise sum; NaN handling via `xt::nansum` if needed |
| `parallel_difference(A, B)` | `A - B` | Element-wise; xtensor overloads `operator-` |
| `parallel_which_max_row(mat)` | `xt::argmax(mat, 1)` | Row-wise argmax → index of max column |
| `parallel_which_min_row(mat)` | `xt::argmin(mat, 1)` | Row-wise argmin |
| `parallel_sd(mat)` | `xt::stddev(mat, 1)` | Population SD (ddof=0); use `xt::stddev(mat, {1})` |
| `parallel_variance(mat)` | `xt::variance(mat, 1)` | Population variance (ddof=0) |
| `parallel_which_max_rolling_quarter(mat)` | Custom `rolling_quarter_argmax(mat)` | 12-position circular rolling window, returns start index |
| `parallel_which_min_rolling_quarter(mat)` | Custom `rolling_quarter_argmin(mat)` | Circular rolling window, returns start index |
| `parallel_average_quarter(mat, idx)` | Custom `quarter_mean(mat, idx)` | Mean of columns `[idx, idx+1, idx+2]` with circular wrap |

### 2.3 Custom Primitive Implementations

#### `rolling_quarter_argmax` / `rolling_quarter_argmin`

```cpp
/// Returns the starting month index (0-based) of the 3-month rolling window
/// whose column sum is maximum (or minimum), with circular wrap.
inline Array1D rolling_quarter_argmax(const Array2D& mat) {
    const std::size_t N = mat.shape(0);
    Array1D result = xt::zeros<float>({N});

    for (std::size_t p = 0; p < N; ++p) {
        float best = -std::numeric_limits<float>::infinity();
        std::size_t best_idx = 0;
        for (std::size_t i = 0; i < 12; ++i) {
            float s = mat(p, i)
                    + mat(p, (i + 1) % 12)
                    + mat(p, (i + 2) % 12);
            if (s > best) { best = s; best_idx = i; }
        }
        result(p) = static_cast<float>(best_idx);
    }
    return result;
}
```

> **OpenMP acceleration:** Add `#pragma omp parallel for` before the outer loop when `XBIOCLIM_USE_OPENMP` is enabled.

#### `quarter_mean`

```cpp
/// Returns row-wise mean of the 3-month window starting at each pixel's index.
inline Array1D quarter_mean(const Array2D& mat, const Array1D& start_idx) {
    const std::size_t N = mat.shape(0);
    Array1D result = xt::zeros<float>({N});

    for (std::size_t p = 0; p < N; ++p) {
        std::size_t i = static_cast<std::size_t>(start_idx(p));
        result(p) = (mat(p, i % 12)
                   + mat(p, (i + 1) % 12)
                   + mat(p, (i + 2) % 12)) / 3.0f;
    }
    return result;
}
```

### 2.4 Unit Test Coverage for Primitives

Each primitive must have an independent unit test in `tests/test_primitives.cpp` covering:

- **Normal case** — known analytical inputs/outputs.
- **Edge case** — all-equal values (ties in argmax/argmin), zeros in `parallel_sd`.
- **NoData handling** — verify that NaN propagation or masking is consistent.

---

## Phase 3 — GDAL I/O Layer

### 3.1 Responsibilities

- **Reader (`GdalReader`):** Opens 48 input GeoTIFFs, reads tiles/blocks into `ClimateBlock` structs, converts integer raw values to float (applying CHELSA scale/offset if needed).
- **Writer (`GdalWriter`):** Creates 19 output GeoTIFFs matching the input spatial metadata (projection, geotransform, NoData).

### 3.2 Tile-Based Processing

To support arbitrarily large global rasters without exhausting RAM, processing uses GDAL's block I/O:

```
for each tile_y in range(n_tiles_y):
    for each tile_x in range(n_tiles_x):
        block = reader.read_block(tile_x, tile_y)   // ClimateBlock
        bio   = compute_bioclim(block)               // BioBlock
        writer.write_block(tile_x, tile_y, bio)
```

Tile size is configurable (default 512×512 pixels). Cloud Optimized GeoTIFFs (COG) are produced by default.

### 3.3 CHELSA Raw-Value Conventions

CHELSA distributes integer-coded rasters with linear scale/offset. The library must support per-variable configuration:

| Variable | CHELSA scale | CHELSA offset | Resulting unit |
|---|---|---|---|
| `tas` | 0.1 | −273.15 (after ×0.1) | °C |
| `tasmax` | 0.1 | −273.15 | °C |
| `tasmin` | 0.1 | −273.15 | °C |
| `pr` | 0.1 | 0 | kg m⁻² month⁻¹ |

> For synthetic test data the raw values are used directly (integer = float), simplifying validation.

---

## Phase 4 — BIO Variable Computation

### 4.1 Compute Function

```cpp
BioBlock compute_bioclim(const ClimateBlock& data);
```

All 19 BIO variables are computed from a single pass over the `ClimateBlock` to maximize cache efficiency. The implementation should minimise temporary allocations by reusing intermediate results (e.g., `diurnal_range`, `rolling pr sums`, `rolling tas means`).

### 4.2 Computation Order (dependency graph)

```
tas, tasmax, tasmin, pr  →  diurnal = tasmax - tasmin
diurnal  →  BIO02
tas      →  BIO01, BIO04
tas      →  BIO05 (via tasmax max), BIO06 (via tasmin min)
BIO05, BIO06  →  BIO07
BIO02, BIO07  →  BIO03
pr       →  BIO12, BIO13, BIO14, BIO15
rolling_argmax(pr) → BIO08 (wettest quarter index) → BIO16
rolling_argmin(pr) → BIO09 (driest quarter index)  → BIO17
rolling_argmax(tas)→ BIO10 (warmest quarter index) → BIO18
rolling_argmin(tas)→ BIO11 (coldest quarter index) → BIO19
BIO08 index + tas → BIO08 (mean tas wettest quarter)
BIO09 index + tas → BIO09 (mean tas driest quarter)
BIO16 index + pr  → BIO16 (mean pr wettest quarter)
...
```

### 4.3 Output GeoTIFF Naming

| BIO | Output filename |
|---|---|
| BIO01 | `bio01.tif` |
| BIO02 | `bio02.tif` |
| … | … |
| BIO19 | `bio19.tif` |

---

## Phase 5 — CLI Tool & Integration

### 5.1 Command-Line Interface

A command-line executable `xbioclim` wraps the library for standalone use:

```
Usage: xbioclim [OPTIONS]

Options:
  --tas    GLOB       Pattern for tas GeoTIFFs   (e.g. "data/tas_??.tif")
  --tasmax GLOB       Pattern for tasmax GeoTIFFs
  --tasmin GLOB       Pattern for tasmin GeoTIFFs
  --pr     GLOB       Pattern for pr GeoTIFFs
  --outdir DIR        Output directory            (default: ".")
  --prefix STR        Output filename prefix      (default: "bio")
  --tile-size N       Tile width/height in pixels (default: 512)
  --threads N         OpenMP thread count         (default: all cores)
  --help              Show this message
```

### 5.2 Slurm HPC Integration

`scripts/run_slurm_array.sh` allows distributing work across tiles:

```bash
#!/bin/bash
#SBATCH --job-name=xbioclim
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --array=0-9  # 10 longitude strips

xbioclim \
  --tas   /data/chelsa/tas/CHELSA_tas_{01..12}_*.tif \
  --tasmax /data/chelsa/tasmax/... \
  --tasmin /data/chelsa/tasmin/... \
  --pr    /data/chelsa/pr/...     \
  --outdir /scratch/$SLURM_JOB_ID \
  --tile-size 1024 \
  --threads $SLURM_CPUS_PER_TASK
```

---

## Phase 6 — GPU Acceleration

### 6.1 Strategy A — CUDA Custom Kernels

Implement a CUDA kernel for the core per-pixel computation loop. Each thread handles one pixel (or a small warp-coalesced block of pixels).

Key considerations:
- Memory transfer: copy the `[N_pixels, 48]` float input matrix to device in one `cudaMemcpy`.
- Kernel maps pixel index → thread; 12-month loops fully unrolled by the compiler.
- Output: copy `[N_pixels, 19]` BIO result matrix back to host.
- Target: **5–6 billion pixels/second** on NVIDIA A100 (limited by HBM2e bandwidth).

### 6.2 Strategy B — OpenMP 4.5 Target Offload with xtensor

```cpp
#pragma omp target teams distribute parallel for \
  map(to: tas_data, tasmax_data, tasmin_data, pr_data) \
  map(from: bio_data)
for (std::size_t p = 0; p < N_pixels; ++p) {
    // per-pixel xtensor views / inline computation
}
```

Target: **2–4 billion pixels/second** — portable across AMD (ROCm) and NVIDIA (NVCC + nvhpc) GPUs.

### 6.3 CMake GPU Options

```cmake
option(XBIOCLIM_USE_CUDA "Enable CUDA kernels" OFF)
option(XBIOCLIM_USE_OPENMP_OFFLOAD "Enable OpenMP GPU offload" OFF)

if(XBIOCLIM_USE_CUDA)
  enable_language(CUDA)
  target_sources(xbioclim PRIVATE src/bioclim_cuda.cu)
  set_target_properties(xbioclim PROPERTIES CUDA_ARCHITECTURES "80;86;90")
endif()
```

---

## Phase 7 — Performance Benchmarking & Release

### 7.1 Benchmarks

Use Google Benchmark or Catch2 `BENCHMARK` blocks to measure:

- **Primitive throughput:** pixels/second for each primitive at N = 10⁶.
- **End-to-end throughput:** pixels/second for all 19 BIO variables.
- **I/O vs compute ratio:** confirm compute dominates for large tiles.

### 7.2 Release Checklist

- [ ] All 19 BIO variables pass analytical tests (Phase 1 values, ε < 1×10⁻³).
- [ ] CI green on Ubuntu 22.04 (GCC 11 + Clang 14).
- [ ] `clang-tidy` and `cppcheck` pass with zero warnings.
- [ ] Doxygen API documentation generated from header comments.
- [ ] COG output validated with `gdalinfo --checksum`.
- [ ] Performance benchmarks documented in `docs/BENCHMARKS.md`.
- [ ] Semantic versioning tag `v0.1.0` created.
- [ ] `vcpkg` and `Conan` recipes contributed upstream.

---

*Last updated: April 2026*
