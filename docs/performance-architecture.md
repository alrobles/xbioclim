# rxbioclim Performance Architecture

This document captures the full performance analysis and optimization vision for
the `rxbioclim` package. It serves as the authoritative reference the team
consults throughout the project.

---

## 1 Current Architecture Analysis

### 1.1 Data Flow

The current data path for `bioclim_raster()` traverses **six serialisation
boundaries** before a single bioclimatic variable reaches the output file:

```
Disk (GeoTIFF / NetCDF)
  → terra::readValues()          # GDAL → R matrix (column-major, double)
    → R for-loop per pixel       # .compute_bioclim_block()
      → S4 dispatch bioclim()    # method resolution per call
        → R primitives           # sd_pop, quarter_argmax, …
    → R result matrix            # 19 columns
  → terra::writeValues()         # R matrix → GDAL → Disk
```

When `ncores > 1`, an additional serialisation step is inserted: each chunk of
the block is **deep-copied and serialised** to a socket-based PSOCK cluster via
`parallel::parLapply`.

The C++ batch path (`bioclim_cpp()`, `ClimateBlock$compute()`) exists but is
**not wired into** the raster processing loop. The R-level
`.compute_bioclim_block()` function in `R/bioclim_raster.R` (lines 46-60) calls
the **R** `bioclim()` generic in a `for` loop instead.

### 1.2 Critical Bottlenecks

| # | Location | Issue | Impact |
|---|----------|-------|--------|
| 1 | `R/bioclim_raster.R` `.compute_bioclim_block()` | R-level `for` loop calls S4 generic `bioclim()` per pixel instead of compiled `bioclim_cpp()` | Orders of magnitude slower; S4 dispatch + R interpreter overhead per pixel |
| 2 | `src/rcpp_bioclim_mod.cpp` lines 72-84 | `ClimateBlock` constructor copies every element one-by-one from R column-major to C++ row-major layout | O(n × 12 × 4) copy operations; prevents zero-copy adaptor use |
| 3 | `src/bioclim.cpp` lines 424-486 | `bioclim_cpp()` creates temporary `NumericVector` objects on every iteration via `NumericVector(tasmax.row(i))` | Heap allocation per pixel per variable; GC pressure |
| 4 | `R/bioclim_raster.R` lines 194-216, 244-256 | `parallel::parLapply` serialises entire matrix chunks to PSOCK worker processes via sockets | Serialisation + deserialisation cost dominates for small chunks; no shared memory |

### 1.3 The Disconnect with xbioclim

`rxbioclim` ships a **self-contained** `src/xbioclim.h` header that reimplements
the algorithm core using `std::vector<double>`. The parent C++ project
[alrobles/xbioclim](https://github.com/alrobles/xbioclim) uses a fundamentally
different technology stack:

| Feature | `rxbioclim/src/xbioclim.h` | `alrobles/xbioclim` (C++ library) |
|---------|---------------------------|-----------------------------------|
| Data type | `std::vector<double>` | `xt::xtensor<float, 2>` |
| File I/O | None (R / terra handles it) | GDAL via `gdal_io.cpp` |
| Vectorisation | Scalar per-pixel loop | Whole-array xtensor expressions |
| Parallelism | None | OpenMP (`primitives.cpp`) |
| GPU | None | CUDA kernels (`bioclim_cuda.cu`) |
| Memory layout | Row-major std::vector | Configurable xtensor layout |

The real xbioclim library already solves every performance problem rxbioclim
has — rxbioclim just doesn't use any of it.

---

## 2 The Zero-R-Memory Vision

### 2.1 Target Architecture

R becomes a **pure declarative layer**. No pixel data ever enters R's heap.

```
┌─────────────────────────────────────────────────────────┐
│                       R (user-facing)                    │
│                                                          │
│  bioclim_engine("tmin_*.tif", "tmax_*.tif",             │
│                 "prec_*.tif",                             │
│                  mask = "region.shp",                     │
│                  output = "bioclim_out.tif",              │
│                  threads = 8)                             │
│                                                          │
│  ← returns: terra::rast() metadata only (no pixel data)  │
└──────────────┬──────────────────────────────────────────┘
               │ XPtr / file paths + parameters
               ▼
┌─────────────────────────────────────────────────────────┐
│              C++ Engine (all data lives here)             │
│                                                          │
│  1. GDAL opens input rasters → windowed tile reads       │
│  2. Optional polygon mask applied at read time           │
│  3. xtensor computes all 19 BIO vars per tile            │
│  4. GDAL writes output tiles directly to disk            │
│                                                          │
│  Memory: 1-2 tiles in RAM at a time                      │
│  Speed: no R serialisation, no GC pressure               │
└─────────────────────────────────────────────────────────┘
```

### 2.2 Comparison Table

| Metric | Current rxbioclim | fastbioclim_check | Proposed engine |
|--------|-------------------|-------------------|-----------------|
| Data path | Disk → Terra → R → Rcpp → R → Terra → Disk | Disk → Terra → R → Rfast → R → Terra → Disk | Disk → GDAL → xtensor → GDAL → Disk |
| R-side copies | 4+ full-raster copies | 2 full-raster copies | **Zero** |
| Parallelism | PSOCK cluster (serialisation) | R-level vectorisation | OpenMP (shared memory) |
| Peak memory | ~6× raster size | ~3× raster size | **1-2 tiles** |
| GPU support | No | No | CUDA (optional) |
| Bottleneck | R interpreter + S4 dispatch | R→C++ boundary | **Disk I/O** (compute-bound eliminated) |

### 2.3 Key Insight

The `xbioclim` C++ project **already has** every component needed:

- **GDAL I/O**: `gdal_io.cpp` — `GdalReader` / `GdalWriter` with tiled windowed reads/writes
- **Vectorised compute**: `bioclim.cpp` — `xt::mean`, `xt::amax` over full arrays
- **OpenMP parallelism**: `primitives.cpp` — `rolling_quarter_argmax`, etc.
- **CUDA GPU**: `bioclim_cuda.cu` — GPU kernels for massive grids
- **CLI pipeline**: `main.cpp` — complete tiled processing pipeline

The task is to **bridge** these into the R package, not to rewrite them.

---

## 3 Components Already Built (across repositories)

### 3.1 alrobles/xbioclim (C++ library)

| File | Purpose | Key types / functions |
|------|---------|-----------------------|
| `src/gdal_io.cpp` / `gdal_io.hpp` | GDAL raster I/O | `GdalReader`, `GdalWriter`, tiled windowed reads, scale/offset decoding, COG output |
| `src/bioclim.cpp` | Vectorised bioclim | Whole-array computation via `xt::mean`, `xt::amax` |
| `src/primitives.cpp` | Parallel primitives | OpenMP `rolling_quarter_argmax`, `rolling_quarter_argmin`, etc. |
| `src/bioclim_cuda.cu` | CUDA GPU kernels | GPU-parallel bioclim computation |
| `src/main.cpp` | CLI entry point | Complete tiled processing pipeline |

### 3.2 alrobles/xtensor (forked tensor library)

Header-only C++ tensor library providing:

- `xt::adapt()` — zero-copy pointer adaptors onto existing memory
- `xt::xtensor_fixed` — stack-allocated tensors (no heap for 12-element arrays)
- Expression templates — lazy evaluation, fused operations
- `xt::partition()` — O(n) partial sort via `std::nth_element`
- Configurable memory layout (row-major / column-major)

### 3.3 alrobles/sf `configure.ac`

Reference implementation for detecting and linking GDAL from an R package's
build system. The pattern can be adapted for rxbioclim's `configure.ac`.

### 3.4 rxbioclim (this package)

| File | Purpose | Notes |
|------|---------|-------|
| `src/xbioclim.h` | Self-contained C++ core | `compute_pixel()`, `compute_bioclim()`, `ClimateBlock`, `BioBlock` |
| `src/bioclim.cpp` | Per-variable + batch Rcpp exports | `bio01_cpp()` … `bio19_cpp()`, `bioclim_cpp()` |
| `src/rcpp_bioclim_mod.cpp` | Rcpp Module `ClimateBlock` class | Column-to-row-major copy in constructor |
| `src/bioclim_model.cpp` | XPtr-based `BioclimModel` | Single-pixel XPtr pattern (proof of concept for engine approach) |
| `R/bioclim_raster.R` | Block-loop raster processing | `.compute_bioclim_block()` — the R for-loop bottleneck |
| `R/BioclimData.R` | S4 class + methods | Dispatches to `*_cpp()` functions for `BioclimData` objects |

---

## 4 Benchmark Targets

| Grid size | Pixels | Target time | Where tested |
|-----------|--------|-------------|--------------|
| Small (100 × 100) | 10,000 | < 0.1 s | CI (every PR) |
| Medium (1000 × 1000) | 1,000,000 | < 1 s | CI (every PR) |
| Large (global 0.5°) | ~260,000+ | Manual comparison vs `fastbioclim_check` | End of project (NOT in CI) |

The small and medium benchmarks use **synthetic rasters** with realistic climate
value ranges (temperature: −30 to 45 °C, precipitation: 0 to 500 mm). The large
grid benchmark requires real-world data and is reserved for manual end-of-project
validation.

See [`benchmarks/README.md`](../benchmarks/README.md) for instructions on
running benchmarks locally.

---

## 5 xtensor Capabilities for Exploitation

### 5.1 Zero-Copy `xt::adapt()`

```cpp
// Non-owning view over R's REAL() memory — zero copies
auto xt_view = xt::adapt(REAL(sexp), size, xt::no_ownership(), shape);
```

### 5.2 Stack-Allocated Fixed Tensors

```cpp
// 96 bytes on the stack for 12-element monthly data — no heap
xt::xtensor_fixed<double, xt::xshape<12>> monthly;
```

### 5.3 O(n) Partial Sort

```cpp
// Find kth element without full sort (for wettest/driest quarter)
xt::partition(quarterly_sums, kth);
```

### 5.4 Expression Templates

All arithmetic produces lazy expression trees evaluated in a single pass:

```cpp
auto bio02 = xt::mean(tasmax - tasmin, {1});  // fused, no intermediate array
```
