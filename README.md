# `bioclim-gpu`

[GitHub Repository](https://github.com/yourusername/bioclim-gpu)  
MIT License  

## High‑Performance Bioclimatic Variables

Fast, parallel computation of WorldClim-style bioclimatic variables (BIO1--BIO19) from monthly temperature and precipitation rasters. Designed for HPC clusters and GPU acceleration -- processes billions of pixels per second on an NVIDIA A100.

### Problem & Approach

Bioclimatic variables are derived per pixel from 12 monthly inputs. The computation is:

- **Embarrassingly parallel** – each pixel independent.
- **Memory‑bound** – ~170 bytes/pixel I/O, ~150 FLOPs/pixel.
- **Ideal for GPUs** – custom CUDA kernels or OpenMP offload achieve 5–6 billion pixels/s on A100.

Two implementation strategies are provided:

1. **Custom CUDA kernels** – maximum performance (near peak memory bandwidth).
2. **`xtensor` + OpenMP target offload** – portable, maintainable, still 2–4 billion pixels/s.

Both use GDAL for tile‑based I/O and support Cloud Optimized GeoTIFFs (COGs).

### Performance Expectations (NVIDIA A100)

| Metric                     | Value             |
|----------------------------|-------------------|
| Peak memory bandwidth      | 1.6 TB/s         |
| Theoretical pixels/s       | 9.4 billion      |
| Realistic (CUDA)           | 5–6 billion pixels/s |
| Realistic (OpenMP offload) | 2–4 billion pixels/s |
| Global 1 km land grid      | ~0.03 seconds (CUDA) |

**Note:** I/O becomes the bottleneck – COGs and fast parallel file systems are essential.

### Dependencies

- C++17 compiler (GCC 9+, Clang 10+, or NVIDIA HPC SDK)
- GDAL (≥3.0) with C++ bindings
- OpenMP (≥4.5) for offload strategy
- CUDA Toolkit (≥11.0) for custom kernels
- `xtensor` (optional, for strategy 2)
- CMake (≥3.15)

### Build Instructions

```bash
git clone https://github.com/yourusername/bioclim-gpu.git
cd bioclim-gpu
mkdir build && cd build
cmake .. -DUSE_CUDA=ON      # for custom CUDA kernels
# or
cmake .. -DUSE_OPENMP_OFFLOAD=ON   # for xtensor+OpenMP
make -j
```

### Usage Example

```bash
# Input: 24 files (tmean_01.tif ... tmean_12.tif, prec_01.tif ... prec_12.tif)
# Output: 19 files (bio1.tif ... bio19.tif)

./bioclim_gpu \
  --tmean monthly_tmean_{01..12}.tif \
  --prec monthly_prec_{01..12}.tif \
  --out_prefix bioclim_ \
  --tile_size 512
```

For Slurm clusters:

```bash
sbatch scripts/run_slurm_array.sh
```

### Full Protocol

The detailed theoretical and implementation protocol is available in `docs/protocol.md`.

### Contributing

Issues and pull requests welcome. See `CONTRIBUTING.md`.

### License

MIT © Angel Luis Robles Fernández
