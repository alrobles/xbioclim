# xbioclim

[![CI](https://github.com/alrobles/xbioclim/actions/workflows/ci.yml/badge.svg)](https://github.com/alrobles/xbioclim/actions/workflows/ci.yml)
MIT License

## Overview

**xbioclim** is a high-performance C++17 library for computing WorldClim-style
bioclimatic variables (BIO1–BIO19) from monthly temperature and precipitation
rasters. It uses [xtensor](https://github.com/xtensor-stack/xtensor) for
vectorised array operations, [GDAL](https://gdal.org/) for geospatial I/O,
and optionally OpenMP for CPU parallelism.

The computation is embarrassingly parallel — each pixel is processed
independently — making it well suited for HPC clusters and GPU offload.

## Dependencies

| Dependency | Minimum version | Notes |
|------------|-----------------|-------|
| C++17 compiler | GCC 9+, Clang 10+ | Required |
| CMake | ≥ 3.15 | Build system |
| GDAL | ≥ 3.0 | Raster I/O (C++ bindings) |
| xtensor | 0.25 | Array primitives (header-only) |
| xtl | 0.7 | Required by xtensor |
| OpenMP | any | Optional — CPU parallelism |

## Build Instructions

```bash
# Clone the repository
git clone https://github.com/alrobles/xbioclim.git
cd xbioclim

# Configure and build
cmake -B build -DCMAKE_BUILD_TYPE=Release -DXBIOCLIM_BUILD_TESTS=ON
cmake --build build --parallel $(nproc)
```

### CMake Options

| Option | Default | Description |
|--------|---------|-------------|
| `XBIOCLIM_BUILD_TESTS` | `ON` | Build the Catch2 test suite |
| `XBIOCLIM_USE_OPENMP` | `ON` | Enable OpenMP CPU parallelism |
| `XBIOCLIM_USE_OPENMP_OFFLOAD` | `OFF` | Enable OpenMP 4.5 GPU offload |
| `XBIOCLIM_USE_CUDA` | `OFF` | Enable CUDA GPU kernels |
| `XBIOCLIM_OPENMP_OFFLOAD_FLAGS` | _(empty)_ | Extra compiler/linker flags for GPU offload target |

### GPU Build (OpenMP 4.5 Target Offload)

xbioclim supports GPU-accelerated computation via OpenMP 4.5 `target` offload.
When `XBIOCLIM_USE_OPENMP_OFFLOAD=ON`, the rolling-quarter and quarter-mean
primitives use `#pragma omp target teams distribute parallel for` to offload
work to an available device (GPU or host fallback).

**NVIDIA GPU (GCC + nvptx):**

```bash
cmake -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DXBIOCLIM_USE_OPENMP_OFFLOAD=ON \
  -DXBIOCLIM_OPENMP_OFFLOAD_FLAGS="-foffload=nvptx-none;-foffload-options=nvptx-none=-misa=sm_80"
cmake --build build --parallel $(nproc)
```

**NVIDIA GPU (Clang/LLVM + nvptx64):**

```bash
cmake -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=clang++ \
  -DXBIOCLIM_USE_OPENMP_OFFLOAD=ON \
  -DXBIOCLIM_OPENMP_OFFLOAD_FLAGS="-fopenmp-targets=nvptx64-nvidia-cuda"
cmake --build build --parallel $(nproc)
```

**AMD GPU (Clang/LLVM + ROCm):**

```bash
cmake -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=clang++ \
  -DXBIOCLIM_USE_OPENMP_OFFLOAD=ON \
  -DXBIOCLIM_OPENMP_OFFLOAD_FLAGS="-fopenmp-targets=amdgcn-amd-amdhsa;-Xopenmp-target=amdgcn-amd-amdhsa;-march=gfx90a"
cmake --build build --parallel $(nproc)
```

**Host-only testing (no GPU required):**

```bash
cmake -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DXBIOCLIM_BUILD_TESTS=ON \
  -DXBIOCLIM_USE_OPENMP_OFFLOAD=ON
cmake --build build --parallel $(nproc)
OMP_TARGET_OFFLOAD=DISABLED ctest --test-dir build --output-on-failure
```

When no GPU is available, set `OMP_TARGET_OFFLOAD=DISABLED` to run target
regions on the host CPU. The `XBIOCLIM_USE_OPENMP_OFFLOAD` tests in the
Catch2 suite use this path in CI.

## Running Tests

```bash
# Generate synthetic test fixtures (built with XBIOCLIM_BUILD_TESTS=ON)
./build/xbioclim_generate_test_data --outdir tests/data

# Run all tests
XBIOCLIM_TEST_DATA=tests/data ctest --test-dir build --output-on-failure
```

## Usage Example

```bash
# CLI tool (when built)
./build/xbioclim_cli \
  --tas tas_{01..12}.tif \
  --tasmax tasmax_{01..12}.tif \
  --tasmin tasmin_{01..12}.tif \
  --pr pr_{01..12}.tif \
  --outdir output/ \
  --prefix bio_ \
  --tile-size 512
```

For Slurm clusters:

```bash
sbatch scripts/run_slurm_array.sh
```

## Dependency Management

### vcpkg

```bash
vcpkg install  # reads vcpkg.json manifest automatically
cmake -B build -DCMAKE_TOOLCHAIN_FILE=$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake
```

### Conan

```bash
conan install . --output-folder=build --build=missing
cmake -B build -DCMAKE_TOOLCHAIN_FILE=build/conan_toolchain.cmake
```

## Documentation

- [Protocol](docs/PROTOCOL.md) — Theoretical background and implementation details
- [Roadmap](docs/ROADMAP.md) — Development phases and release checklist
- [Benchmarks](docs/BENCHMARKS.md) — CPU performance benchmarks for v0.1.0

## Project Structure

```
xbioclim/
├── include/xbioclim/   # Public headers (primitives, bioclim, gdal_io, version)
├── src/                # Library implementation
├── tests/              # Catch2 test suite
├── tools/              # Helper scripts (test data generation)
├── docs/               # Protocol and roadmap documentation
├── cmake/              # CMake helper modules
└── scripts/            # HPC / Slurm job scripts
```

## License

MIT © Angel Luis Robles Fernández — see [LICENSE](LICENSE) for details.
