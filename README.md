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

## Running Tests

```bash
# Generate synthetic test fixtures (requires Python 3 + GDAL)
python3 tools/generate_test_data.py --outdir tests/data

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

## Documentation

- [Protocol](docs/PROTOCOL.md) — Theoretical background and implementation details
- [Roadmap](docs/ROADMAP.md) — Development phases and release checklist

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
