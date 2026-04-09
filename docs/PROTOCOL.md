# xbioclim — Technical Protocol

> **Scope:** Complete algorithmic specification for computing the 19 CHELSA/WorldClim bioclimatic variables (BIO01–BIO19) using a pure C++17 pipeline with GDAL for raster I/O and xtensor for array algebra.

---

## Table of Contents

1. [Input Variables](#1-input-variables)
2. [Output Variables — Definitions & Formulas](#2-output-variables--definitions--formulas)
3. [Data Layout & Memory Model](#3-data-layout--memory-model)
4. [Core Algorithmic Primitives](#4-core-algorithmic-primitives)
5. [BIO Variable Computation Pipeline](#5-bio-variable-computation-pipeline)
6. [GDAL I/O Protocol](#6-gdal-io-protocol)
7. [NoData & Edge-Case Handling](#7-nodata--edge-case-handling)
8. [Numerical Precision](#8-numerical-precision)
9. [Performance Model](#9-performance-model)
10. [Reference Implementation Mapping (Rcpp → xtensor)](#10-reference-implementation-mapping-rcpp--xtensor)

---

## 1. Input Variables

### 1.1 Source

Monthly climate rasters are sourced from **CHELSA v2.1**  
(https://www.chelsa-climate.org/datasets/chelsa_monthly).

### 1.2 Variables

| Symbol | Long name | Unit | Temporal coverage |
|---|---|---|---|
| `tas` | Daily mean near-surface air temperature | K | 12 monthly layers |
| `tasmax` | Daily maximum near-surface air temperature | K | 12 monthly layers |
| `tasmin` | Daily minimum near-surface air temperature | K | 12 monthly layers |
| `pr` | Precipitation | kg m⁻² month⁻¹ | 12 monthly layers |

> **Note:** `pr` includes both liquid (rain) and solid (snow) phases. `tas`, `tasmax`, and `tasmin` are measured at approximately 2 m above the surface.

### 1.3 CHELSA Integer Encoding

CHELSA distributes data as integer GeoTIFFs with a linear scale/offset encoding:

```
physical_value = raw_integer × scale + offset
```

| Variable | Scale | Offset | Physical unit after conversion |
|---|---|---|---|
| `tas` | 0.1 | −273.15 (applied after scaling) | °C |
| `tasmax` | 0.1 | −273.15 | °C |
| `tasmin` | 0.1 | −273.15 | °C |
| `pr` | 0.1 | 0 | kg m⁻² month⁻¹ |

Conversion for temperature in practice:

```
tas_celsius[m] = tas_raw[m] × 0.1 − 273.15
```

### 1.4 Spatial Reference

- **CRS:** EPSG:4326 (WGS 84 geographic)
- **Resolution:** 30 arc-seconds (~1 km at the equator) for CHELSA V2.1
- **Extent:** Global land cover

---

## 2. Output Variables — Definitions & Formulas

All indices below are 1-based month indices (1 = January … 12 = December). Quarter indices refer to the **start** month of a rolling 3-month window; the quarter window wraps circularly (e.g., start=12 covers months 12, 1, 2).

Let:
- **T[m]** = `tas` (monthly mean temperature in °C) for month *m*
- **Tx[m]** = `tasmax` for month *m*  
- **Tn[m]** = `tasmin` for month *m*  
- **P[m]** = `pr` (monthly precipitation) for month *m*  
- **DR[m]** = Tx[m] − Tn[m] (diurnal range for month *m*)

### BIO01 — Mean Annual Near-Surface Air Temperature (°C)

$$\text{BIO01} = \frac{1}{12} \sum_{m=1}^{12} T[m]$$

### BIO02 — Mean Diurnal Range (°C)

$$\text{BIO02} = \frac{1}{12} \sum_{m=1}^{12} \bigl( Tx[m] - Tn[m] \bigr)$$

### BIO03 — Isothermality (dimensionless, expressed as %)

$$\text{BIO03} = 100 \times \frac{\text{BIO02}}{\text{BIO07}}$$

Isothermality measures how large the mean daily temperature oscillation is compared to the annual temperature range. A value of 100 means the daily range equals the annual range.

### BIO04 — Temperature Seasonality (°C × 100)

$$\text{BIO04} = 100 \times \sigma_{\text{pop}}(T[1], \ldots, T[12])$$

where $\sigma_{\text{pop}}$ is the **population standard deviation** (denominator *N* = 12, not *N*−1).

$$\sigma_{\text{pop}} = \sqrt{ \frac{1}{12} \sum_{m=1}^{12} \bigl(T[m] - \overline{T}\bigr)^2 }$$

### BIO05 — Maximum Temperature of Warmest Month (°C)

$$\text{BIO05} = \max_{m} \bigl( Tx[m] \bigr)$$

### BIO06 — Minimum Temperature of Coldest Month (°C)

$$\text{BIO06} = \min_{m} \bigl( Tn[m] \bigr)$$

### BIO07 — Annual Temperature Range (°C)

$$\text{BIO07} = \text{BIO05} - \text{BIO06}$$

### BIO08 — Mean Temperature of Wettest Quarter (°C)

Let $q^*$ be the start month of the 3-month window maximising total precipitation:

$$q^* = \arg\max_{q \in \{1,\ldots,12\}} \; \sum_{k=0}^{2} P[(q + k - 1) \bmod 12 + 1]$$

$$\text{BIO08} = \frac{1}{3} \sum_{k=0}^{2} T[(q^* + k - 1) \bmod 12 + 1]$$

### BIO09 — Mean Temperature of Driest Quarter (°C)

$$q^* = \arg\min_{q} \sum_{k=0}^{2} P[\ldots]$$

$$\text{BIO09} = \frac{1}{3} \sum_{k=0}^{2} T[\ldots]$$

### BIO10 — Mean Temperature of Warmest Quarter (°C)

$$q^* = \arg\max_{q} \frac{1}{3} \sum_{k=0}^{2} T[\ldots]$$

$$\text{BIO10} = \frac{1}{3} \sum_{k=0}^{2} T[(q^* + k - 1) \bmod 12 + 1]$$

### BIO11 — Mean Temperature of Coldest Quarter (°C)

$$q^* = \arg\min_{q} \frac{1}{3} \sum_{k=0}^{2} T[\ldots]$$

$$\text{BIO11} = \frac{1}{3} \sum_{k=0}^{2} T[\ldots]$$

### BIO12 — Annual Precipitation (kg m⁻² year⁻¹)

$$\text{BIO12} = \sum_{m=1}^{12} P[m]$$

### BIO13 — Precipitation of Wettest Month (kg m⁻² month⁻¹)

$$\text{BIO13} = \max_{m} P[m]$$

### BIO14 — Precipitation of Driest Month (kg m⁻² month⁻¹)

$$\text{BIO14} = \min_{m} P[m]$$

### BIO15 — Precipitation Seasonality (coefficient of variation, %)

$$\text{BIO15} = 100 \times \frac{\sigma_{\text{pop}}(P[1],\ldots,P[12])}{\overline{P}}$$

> **Warning:** When $\overline{P} = 0$ (permanently dry grid cells), BIO15 is undefined. Assign NoData.

### BIO16 — Precipitation of Wettest Quarter (kg m⁻² month⁻¹)

Uses the same wettest quarter start index $q^*$ determined from rolling pr sums (as in BIO08):

$$\text{BIO16} = \frac{1}{3} \sum_{k=0}^{2} P[(q^* + k - 1) \bmod 12 + 1]$$

> **Note:** BIO16 is the *mean* monthly precipitation of the wettest quarter, not the sum.

### BIO17 — Precipitation of Driest Quarter (kg m⁻² month⁻¹)

$$\text{BIO17} = \frac{1}{3} \sum_{k=0}^{2} P[(q^* + k - 1) \bmod 12 + 1]$$

where $q^*$ is the driest-quarter start index (same as BIO09).

### BIO18 — Precipitation of Warmest Quarter (kg m⁻² month⁻¹)

Uses the warmest-quarter start index from BIO10:

$$\text{BIO18} = \frac{1}{3} \sum_{k=0}^{2} P[(q^* + k - 1) \bmod 12 + 1]$$

### BIO19 — Precipitation of Coldest Quarter (kg m⁻² month⁻¹)

Uses the coldest-quarter start index from BIO11:

$$\text{BIO19} = \frac{1}{3} \sum_{k=0}^{2} P[(q^* + k - 1) \bmod 12 + 1]$$

---

## 3. Data Layout & Memory Model

### 3.1 Tile-Based Flattening

Processing is tile-based. A tile of *W* × *H* pixels is flattened in **row-major order** to a 1-D pixel index:

```
pixel_idx = row × W + col      (0-based)
N_pixels  = W × H
```

### 3.2 Input Array Shape

Each of the four variables (`tas`, `tasmax`, `tasmin`, `pr`) is represented as a 2-D array:

```
shape: [N_pixels, 12]     (float32)
```

So the 12 monthly values for pixel *p* are contiguous in memory at `arr[p, 0..11]`.  
This layout is **optimal for the per-pixel computation** because each BIO variable requires all 12 months of a single pixel — no cross-pixel communication.

### 3.3 Memory Footprint per Tile

For a tile of size 512 × 512 = 262,144 pixels, with float32 storage:

| Array | Shape | Size |
|---|---|---|
| `tas` | [262144, 12] | 12.6 MB |
| `tasmax` | [262144, 12] | 12.6 MB |
| `tasmin` | [262144, 12] | 12.6 MB |
| `pr` | [262144, 12] | 12.6 MB |
| **Total input** | | **~50.3 MB** |
| BIO01–BIO19 output | [262144, 19] | ~20 MB |
| **Total per tile** | | **~70 MB** |

This easily fits in L3 cache on modern CPUs when processing smaller tiles (e.g., 128 × 128 = 16 MB input).

### 3.4 C++ Type Aliases

```cpp
// include/xbioclim/primitives.hpp

#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xmath.hpp>

namespace xbioclim {

/// 2-D array: shape [N_pixels, 12]
using Array2D = xt::xtensor<float, 2>;

/// 1-D array: shape [N_pixels]
using Array1D = xt::xtensor<float, 1>;

/// Integer index array: shape [N_pixels]
using IndexArray = xt::xtensor<std::size_t, 1>;

/// Bundle of four input variable blocks
struct ClimateBlock {
    Array2D tas;    ///< Monthly mean temperature  [N_pixels, 12]
    Array2D tasmax; ///< Monthly max temperature   [N_pixels, 12]
    Array2D tasmin; ///< Monthly min temperature   [N_pixels, 12]
    Array2D pr;     ///< Monthly precipitation     [N_pixels, 12]
    std::size_t n_pixels() const { return tas.shape(0); }
};

/// Bundle of 19 output BIO variable arrays
struct BioBlock {
    Array1D bio01, bio02, bio03, bio04, bio05,
            bio06, bio07, bio08, bio09, bio10,
            bio11, bio12, bio13, bio14, bio15,
            bio16, bio17, bio18, bio19;
};

} // namespace xbioclim
```

---

## 4. Core Algorithmic Primitives

All primitives operate on `Array2D` (shape `[N_pixels, 12]`) and return `Array1D` (shape `[N_pixels]`), implementing the row-wise operation for all pixels simultaneously.

### 4.1 `row_mean` — Row-wise Mean

**Formula:** $\bar{x}_p = \frac{1}{12}\sum_{m=0}^{11} A[p, m]$

```cpp
Array1D row_mean(const Array2D& A) {
    return xt::mean(A, /*axis=*/1);
}
```

**xtensor note:** `xt::mean(A, {1})` reduces over axis 1 (months), producing shape `[N_pixels]`.

### 4.2 `row_sum` — Row-wise Sum

**Formula:** $s_p = \sum_{m=0}^{11} A[p, m]$

```cpp
Array1D row_sum(const Array2D& A) {
    return xt::sum(A, /*axis=*/1);
}
```

For NoData-aware summation: replace NaN values with 0 before summing, or use a masked approach.

### 4.3 `row_max` / `row_min` — Row-wise Maximum / Minimum

```cpp
Array1D row_max(const Array2D& A) { return xt::amax(A, /*axis=*/1); }
Array1D row_min(const Array2D& A) { return xt::amin(A, /*axis=*/1); }
```

### 4.4 `row_std` — Row-wise Population Standard Deviation

**Formula:** $\sigma_p = \sqrt{\frac{1}{12} \sum_{m=0}^{11}(A[p,m] - \bar{A}_p)^2}$

```cpp
Array1D row_std(const Array2D& A) {
    // ddof=0 for population SD (default in xtensor)
    return xt::stddev(A, /*axis=*/1);
}
```

**Welford's online algorithm** (numerically stable, O(1) extra memory) is the preferred implementation for the manual scalar path (e.g., GPU kernels) — see §4.7.

### 4.5 `elementwise_diff` — Element-wise Difference

```cpp
Array2D elementwise_diff(const Array2D& A, const Array2D& B) {
    return A - B;  // xtensor overloads operator-
}
```

### 4.6 `rolling_quarter_argmax` / `rolling_quarter_argmin`

Computes the 0-based start index of the 3-month circular rolling window maximising (or minimising) the sum across all 12 possible windows.

**Algorithm:**

```
for each pixel p:
    best_sum  = -∞
    best_idx  = 0
    for i in 0..11:
        s = A[p, i] + A[p, (i+1)%12] + A[p, (i+2)%12]
        if s > best_sum:
            best_sum = s
            best_idx = i
    result[p] = best_idx
```

```cpp
IndexArray rolling_quarter_argmax(const Array2D& A) {
    const std::size_t N = A.shape(0);
    IndexArray result = xt::zeros<std::size_t>({N});

    #ifdef XBIOCLIM_USE_OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (std::size_t p = 0; p < N; ++p) {
        float best = std::numeric_limits<float>::lowest();
        std::size_t bi = 0;
        for (std::size_t i = 0; i < 12; ++i) {
            float s = A(p, i) + A(p, (i+1)%12) + A(p, (i+2)%12);
            if (s > best) { best = s; bi = i; }
        }
        result(p) = bi;
    }
    return result;
}

IndexArray rolling_quarter_argmin(const Array2D& A) {
    const std::size_t N = A.shape(0);
    IndexArray result = xt::zeros<std::size_t>({N});

    #ifdef XBIOCLIM_USE_OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (std::size_t p = 0; p < N; ++p) {
        float best = std::numeric_limits<float>::max();
        std::size_t bi = 0;
        for (std::size_t i = 0; i < 12; ++i) {
            float s = A(p, i) + A(p, (i+1)%12) + A(p, (i+2)%12);
            if (s < best) { best = s; bi = i; }
        }
        result(p) = bi;
    }
    return result;
}
```

> **Tie-breaking:** When two windows have equal sums, the **first** (lowest month index) is chosen. This matches the Rcpp PoC behaviour.

### 4.7 `quarter_mean` — Mean of a Rolling 3-Month Window

Given a start-index array (output of `rolling_quarter_argmax/min`):

```cpp
Array1D quarter_mean(const Array2D& A, const IndexArray& starts) {
    const std::size_t N = A.shape(0);
    Array1D result = xt::zeros<float>({N});

    #ifdef XBIOCLIM_USE_OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (std::size_t p = 0; p < N; ++p) {
        const std::size_t i = starts(p);
        result(p) = (A(p, i % 12)
                   + A(p, (i+1) % 12)
                   + A(p, (i+2) % 12)) / 3.0f;
    }
    return result;
}
```

### 4.8 Welford's Online Algorithm (for GPU scalar path)

For environments where `xt::stddev` is not available (e.g., hand-written CUDA kernels):

```cuda
__device__ float welford_std(const float* arr, int n) {
    float mean = 0.0f, M2 = 0.0f;
    for (int k = 0; k < n; ++k) {
        float delta = arr[k] - mean;
        mean += delta / (k + 1);
        M2   += delta * (arr[k] - mean);
    }
    return sqrtf(M2 / n);  // population SD (ddof=0)
}
```

---

## 5. BIO Variable Computation Pipeline

### 5.1 Full Compute Function

```cpp
// include/xbioclim/bioclim.hpp
namespace xbioclim {
    BioBlock compute_bioclim(const ClimateBlock& data);
}
```

### 5.2 Step-by-Step Implementation

```cpp
BioBlock xbioclim::compute_bioclim(const ClimateBlock& data) {
    BioBlock bio;

    // --- Intermediates ---
    // Diurnal range: [N_pixels, 12]
    Array2D diurnal = data.tasmax - data.tasmin;

    // --- Simple statistics ---
    bio.bio01 = row_mean(data.tas);            // mean(T)
    bio.bio02 = row_mean(diurnal);             // mean(Tx-Tn)
    bio.bio05 = row_max(data.tasmax);          // max(Tx)
    bio.bio06 = row_min(data.tasmin);          // min(Tn)
    bio.bio07 = bio.bio05 - bio.bio06;         // BIO05 - BIO06
    bio.bio03 = 100.0f * bio.bio02 / bio.bio07;// 100 * BIO02/BIO07
    bio.bio04 = 100.0f * row_std(data.tas);   // 100 * sd(T)
    bio.bio12 = row_sum(data.pr);             // sum(P)
    bio.bio13 = row_max(data.pr);             // max(P)
    bio.bio14 = row_min(data.pr);             // min(P)

    // BIO15: precipitation CV
    Array1D pr_mean = row_mean(data.pr);
    Array1D pr_std  = row_std(data.pr);
    bio.bio15 = 100.0f * pr_std / pr_mean;    // 100 * sd(P)/mean(P)

    // --- Rolling quarter indices ---
    IndexArray wet_q  = rolling_quarter_argmax(data.pr);   // wettest by pr
    IndexArray dry_q  = rolling_quarter_argmin(data.pr);   // driest  by pr
    IndexArray warm_q = rolling_quarter_argmax(data.tas);  // warmest by tas
    IndexArray cold_q = rolling_quarter_argmin(data.tas);  // coldest by tas

    // --- BIO08/09: mean temperature of wettest/driest quarter ---
    bio.bio08 = quarter_mean(data.tas, wet_q);
    bio.bio09 = quarter_mean(data.tas, dry_q);

    // --- BIO10/11: mean temperature of warmest/coldest quarter ---
    bio.bio10 = quarter_mean(data.tas, warm_q);
    bio.bio11 = quarter_mean(data.tas, cold_q);

    // --- BIO16/17: mean precipitation of wettest/driest quarter ---
    bio.bio16 = quarter_mean(data.pr, wet_q);
    bio.bio17 = quarter_mean(data.pr, dry_q);

    // --- BIO18/19: mean precipitation of warmest/coldest quarter ---
    bio.bio18 = quarter_mean(data.pr, warm_q);
    bio.bio19 = quarter_mean(data.pr, cold_q);

    return bio;
}
```

### 5.3 Dependency Graph

```
tasmax, tasmin ──► diurnal ──────────────────────► BIO02
                                                      │
tas ──────────────────────────────────────────────► BIO01
                                                      │
tas ──► row_std ──────────────────────────────────► BIO04
                                                      │
tasmax ──► row_max ───────────────────────────────► BIO05
                                                      │
tasmin ──► row_min ───────────────────────────────► BIO06
                                                      │
BIO05, BIO06 ────────────────────────────────────► BIO07
                                                      │
BIO02, BIO07 ────────────────────────────────────► BIO03
                                                      │
pr ──► row_sum ───────────────────────────────────► BIO12
pr ──► row_max ───────────────────────────────────► BIO13
pr ──► row_min ───────────────────────────────────► BIO14
pr ──► row_std, row_mean ────────────────────────► BIO15

pr ──► rolling_quarter_argmax ──► wet_q ──► quarter_mean(tas) ──► BIO08
                                        └──► quarter_mean(pr)  ──► BIO16

pr ──► rolling_quarter_argmin ──► dry_q ──► quarter_mean(tas) ──► BIO09
                                        └──► quarter_mean(pr)  ──► BIO17

tas ──► rolling_quarter_argmax ──► warm_q ──► quarter_mean(tas) ──► BIO10
                                          └──► quarter_mean(pr)  ──► BIO18

tas ──► rolling_quarter_argmin ──► cold_q ──► quarter_mean(tas) ──► BIO11
                                          └──► quarter_mean(pr)  ──► BIO19
```

---

## 6. GDAL I/O Protocol

### 6.1 Opening Input Files

All 48 input GeoTIFFs must be opened in read-only mode. Before processing begins, verify:

1. All files share the same **spatial reference** (CRS), **geotransform**, and **raster dimensions**.
2. All files have exactly **1 band**.
3. The **NoData** value is consistent across all layers of the same variable.

```cpp
GDALDataset* ds = (GDALDataset*) GDALOpen(filename.c_str(), GA_ReadOnly);
if (!ds) throw std::runtime_error("Cannot open: " + filename);
```

### 6.2 Block Reading Strategy

GDAL block-aligned I/O maximises throughput. Read blocks using `RasterIO`:

```cpp
CPLErr err = band->RasterIO(
    GF_Read,
    x_off, y_off,    // pixel offset within the raster
    block_w, block_h, // number of pixels to read
    buffer.data(),   // output buffer
    block_w, block_h,
    GDT_Float32,     // convert to float32 on read
    0, 0             // pixel/line space (0 = default tightly packed)
);
```

Use `GDALDataset::GetRasterXSize()` / `GetRasterYSize()` to determine full extent; iterate blocks with configurable tile size.

### 6.3 Writing Output Files

Create output GeoTIFFs using the **GeoTIFF** driver with COG-compatible options:

```cpp
GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("GTiff");
char** opts = nullptr;
opts = CSLSetNameValue(opts, "COMPRESS", "LZW");
opts = CSLSetNameValue(opts, "TILED", "YES");
opts = CSLSetNameValue(opts, "BLOCKXSIZE", "512");
opts = CSLSetNameValue(opts, "BLOCKYSIZE", "512");
opts = CSLSetNameValue(opts, "BIGTIFF", "IF_SAFER");

GDALDataset* out_ds = drv->Create(
    output_path.c_str(),
    raster_x_size, raster_y_size,
    1,              // single band
    GDT_Float32,
    opts
);
out_ds->SetProjection(srs_wkt.c_str());
out_ds->SetGeoTransform(gt);
out_ds->GetRasterBand(1)->SetNoDataValue(NODATA_VALUE);
```

### 6.4 Class Interface

```cpp
// include/xbioclim/gdal_io.hpp
namespace xbioclim {

class GdalReader {
public:
    explicit GdalReader(const std::vector<std::string>& tas_files,
                        const std::vector<std::string>& tasmax_files,
                        const std::vector<std::string>& tasmin_files,
                        const std::vector<std::string>& pr_files);

    /// Read a tile into a ClimateBlock (pixel values as float32).
    ClimateBlock read_block(int x_off, int y_off,
                            int block_w, int block_h) const;

    int raster_x_size() const;
    int raster_y_size() const;
    std::string projection_wkt() const;
    std::array<double, 6> geotransform() const;
};

class GdalWriter {
public:
    GdalWriter(const std::string& output_dir,
               const std::string& prefix,
               int x_size, int y_size,
               const std::string& projection_wkt,
               const std::array<double, 6>& geotransform);

    /// Write a BioBlock tile to the output GeoTIFFs.
    void write_block(int x_off, int y_off, const BioBlock& bio);

    void finalize();  ///< Flush and close all open datasets
};

} // namespace xbioclim
```

---

## 7. NoData & Edge-Case Handling

### 7.1 NoData Propagation

If **any** of the 12 monthly inputs for a pixel is NoData (typically −9999 as integer or NaN as float):

- All 19 BIO outputs for that pixel are set to the output NoData value (`std::numeric_limits<float>::quiet_NaN()` or a configurable sentinel like −9999.0f).
- A per-pixel validity mask is computed once at block load time:

```cpp
// Mark pixel as invalid if any of the 48 input values is NoData
Array1D valid = xt::ones<float>({N_pixels});
for (int v = 0; v < 4; ++v)
    for (int m = 0; m < 12; ++m)
        valid = valid * xt::cast<float>(
            xt::col(input_array[v], m) != NODATA_FLOAT);
```

- After computing all BIO values, apply the mask:

```cpp
bio.bio01 = xt::where(valid, bio.bio01, NODATA_FLOAT);
// ... repeat for bio02–bio19
```

### 7.2 BIO15 — Division by Zero

When mean precipitation $\overline{P} = 0$, BIO15 is undefined. The implementation must guard against division by zero:

```cpp
bio.bio15 = xt::where(pr_mean > 0.0f,
                      100.0f * pr_std / pr_mean,
                      NODATA_FLOAT);
```

### 7.3 BIO03 — Division by Zero

When `BIO07 = 0` (no annual temperature range), BIO03 is undefined:

```cpp
bio.bio03 = xt::where(bio.bio07 != 0.0f,
                      100.0f * bio.bio02 / bio.bio07,
                      NODATA_FLOAT);
```

### 7.4 Rolling Quarter Tie-Breaking

When multiple quarters share the same minimum or maximum value, always select the **lowest start-month index** (0-based). This ensures deterministic output identical to the Rcpp PoC.

---

## 8. Numerical Precision

### 8.1 Floating-Point Format

All internal computations use **float32** (single precision). This is consistent with CHELSA raster storage and with the performance goals (float32 is twice as fast as float64 on modern SIMD and GPU hardware).

Final output GeoTIFFs are also stored as **GDT_Float32**.

### 8.2 Precision Expectations

| Variable | Relative error target |
|---|---|
| BIO01, BIO02, BIO05, BIO06, BIO07, BIO12, BIO13, BIO14 | < 1 × 10⁻⁵ (exact arithmetic) |
| BIO03, BIO04, BIO10, BIO11, BIO15, BIO08, BIO09, BIO16–BIO19 | < 1 × 10⁻³ (division / sqrt) |

Test tolerances use absolute difference (not relative), as some values are near zero:

```cpp
CHECK_THAT(bio.bio01[0], WithinAbs(expected, 1e-3f));
```

### 8.3 Integer Overflow

CHELSA raw values are Int16 (range −32768 to 32767). When computing sums across 12 months, overflow can occur in integer arithmetic (e.g., 12 × 32767 > 32767). Always convert to float32 before accumulation:

```cpp
float val = static_cast<float>(raw_int16) * scale + offset;
```

---

## 9. Performance Model

### 9.1 Arithmetic Intensity

For a single pixel, the computation performs approximately:

| Stage | Operations |
|---|---|
| Load 48 values (4 vars × 12 months) | 48 reads |
| Diurnal range: Tx − Tn | 12 subtractions |
| BIO01,02,04,12,15: mean/sum/std ×3 | ~84 FLOPs |
| BIO05,06,13,14: max/min ×4 | 44 comparisons |
| Rolling quarters (4×): 12×3 comparisons + additions | 144 FLOPs |
| Quarter means (8×): 3 adds + 1 div | 32 FLOPs |
| NoData mask application | 19 selects |
| Store 19 values | 19 writes |
| **Total** | **~150 FLOPs, ~67 bytes I/O** |

**Arithmetic intensity:** ~150 FLOPs / 67 bytes ≈ **2.2 FLOP/byte**.  
This is below the ridge point of all modern GPUs (typically 10–50 FLOP/byte), confirming the computation is **memory-bandwidth bound**.

### 9.2 Throughput Targets

| Hardware | Bandwidth | Expected pixel/s |
|---|---|---|
| Laptop CPU (DDR4-3200, 4-ch) | 50 GB/s | ~750 M pixels/s |
| HPC node CPU (DDR5-6400, 8-ch) | 100 GB/s | ~1.5 B pixels/s |
| NVIDIA A100 (HBM2e) | 2 TB/s | ~30 B pixels/s (theoretical) |
| NVIDIA A100 (realistic, with kernel overhead) | — | 5–6 B pixels/s |
| AMD MI250X (HBM2e) | 3.2 TB/s | ~48 B pixels/s (theoretical) |

### 9.3 Global 30-arc-second Land Grid

The CHELSA 30-arc-second land mask has approximately **1.57 billion land pixels**. At 5 B pixels/s (A100 CUDA), total processing time ≈ **0.31 seconds** (compute) + I/O overhead.

---

## 10. Reference Implementation Mapping (Rcpp → xtensor)

The original proof-of-concept ([alrobles/fastBioClim](https://github.com/alrobles/fastBioClim)) used Rcpp + RcppParallel with Terra block I/O. The following table maps each PoC primitive to its xtensor equivalent:

| PoC Primitive (Rcpp) | xbioclim Primitive (xtensor) | Notes |
|---|---|---|
| `parallel_average(mat)` | `xt::mean(mat, 1)` | Axis-1 (months) reduction |
| `parallel_sum(mat)` | `xt::sum(mat, 1)` | Use `xt::nansum` for NaN tolerance |
| `parallel_difference(A, B)` | `A - B` | Element-wise; no explicit call needed |
| `parallel_which_max_row(mat)` | `xt::amax(mat, 1)` | Returns values; use `xt::argmax` for indices |
| `parallel_which_min_row(mat)` | `xt::amin(mat, 1)` | Returns values |
| `parallel_sd(mat)` | `xt::stddev(mat, 1)` | Population SD, ddof=0 (xtensor default) |
| `parallel_variance(mat)` | `xt::variance(mat, 1)` | Population variance, ddof=0 |
| `parallel_which_max_rolling_quarter(mat)` | `rolling_quarter_argmax(mat)` | Custom loop; returns 0-based start index |
| `parallel_which_min_rolling_quarter(mat)` | `rolling_quarter_argmin(mat)` | Custom loop |
| `parallel_average_quarter(mat, idx)` | `quarter_mean(mat, idx)` | Circular 3-month window mean |

### Key Differences from Rcpp PoC

1. **Thread model:** Rcpp used `tbb::parallel_for` (Intel TBB) via RcppParallel. The xtensor version uses `#pragma omp parallel for` (OpenMP), which is more portable on HPC clusters and directly supported by GPU offload pragmas.

2. **Memory layout:** R matrices are column-major (Fortran order). The xtensor implementation uses **row-major** (C order) with shape `[N_pixels, 12]`, so each pixel's 12 values are contiguous — this is better for sequential per-pixel access patterns and GPU coalescing.

3. **I/O:** Terra (R package) handled raster I/O using R's native GDAL bindings. The C++ version uses GDAL's C++ API directly, enabling fine-grained control over tile sizes, compression, and Cloud Optimized GeoTIFF output.

4. **Scale/offset:** The Rcpp PoC operated on pre-converted floating-point values passed from R. The C++ library handles CHELSA integer encoding internally via the `GdalReader`, applying the documented scale and offset at load time.

5. **NoData handling:** RcppParallel treated `NA_real_` as IEEE NaN. The xtensor implementation uses explicit masking to propagate NoData efficiently without relying on NaN arithmetic (which can be slower on some hardware).

---

## Appendix A — BIO Variable Quick Reference

| Code | Name | Input | Unit |
|---|---|---|---|
| BIO01 | Mean Annual Near-Surface Air Temperature | tas | °C |
| BIO02 | Mean Diurnal Range | tasmax, tasmin | °C |
| BIO03 | Isothermality | BIO02, BIO07 | % |
| BIO04 | Temperature Seasonality | tas | °C×100 |
| BIO05 | Max Temperature of Warmest Month | tasmax | °C |
| BIO06 | Min Temperature of Coldest Month | tasmin | °C |
| BIO07 | Annual Temperature Range | BIO05, BIO06 | °C |
| BIO08 | Mean Temperature of Wettest Quarter | tas, pr | °C |
| BIO09 | Mean Temperature of Driest Quarter | tas, pr | °C |
| BIO10 | Mean Temperature of Warmest Quarter | tas | °C |
| BIO11 | Mean Temperature of Coldest Quarter | tas | °C |
| BIO12 | Annual Precipitation | pr | kg m⁻² yr⁻¹ |
| BIO13 | Precipitation of Wettest Month | pr | kg m⁻² month⁻¹ |
| BIO14 | Precipitation of Driest Month | pr | kg m⁻² month⁻¹ |
| BIO15 | Precipitation Seasonality (CV) | pr | % |
| BIO16 | Precipitation of Wettest Quarter | pr | kg m⁻² month⁻¹ |
| BIO17 | Precipitation of Driest Quarter | pr | kg m⁻² month⁻¹ |
| BIO18 | Precipitation of Warmest Quarter | pr, tas | kg m⁻² month⁻¹ |
| BIO19 | Precipitation of Coldest Quarter | pr, tas | kg m⁻² month⁻¹ |

---

## Appendix B — CHELSA References

- CHELSA V2.1 dataset: https://www.chelsa-climate.org/downloads/
- Monthly climate variables: https://www.chelsa-climate.org/datasets/chelsa_monthly
- Bioclimatic variables: https://www.chelsa-climate.org/datasets/chelsa_bioclim
- Technical specification: Karger et al. (2017) *Scientific Data* 4, 170122. https://doi.org/10.1038/sdata.2017.122
- WorldClim bioclimatic variable definitions: https://www.worldclim.org/data/bioclim.html

---

*Last updated: April 2026*
