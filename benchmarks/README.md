# xbioclim Benchmarks

Standardised benchmark suite for tracking `xbioclim` performance across
changes.

---

## Quick Start

```bash
# From the repository root (package must be installed first)
Rscript benchmarks/run_benchmarks.R
```

For accurate timings, install the package from an optimized build. If you ran `roxygen2::roxygenise()` first, run `make clean` (or `Rscript tools/clean-obj.R`) and reinstall with `R CMD INSTALL --configure-args='--without-cuda' .` before benchmarking. Better still, use `R CMD build` and install from the tarball.

The script exits with code **0** if all targets are met and **1** otherwise,
making it suitable for CI integration.

---

## What is Benchmarked

| Benchmark | What it measures |
|-----------|------------------|
| `bioclim_raster()` | Full raster pipeline: Terra I/O → computation → output |
| `bioclim_cpp()` | C++ batch function (matrix → matrix, no Terra overhead) |
| `ClimateBlock$compute()` | Rcpp Module path (includes column→row-major copy) |
| Individual `bio*_cpp()` | All 19 individual C++ functions called sequentially |

Each benchmark is run on two synthetic grids:

| Grid | Dimensions | Pixels | Target (median) |
|------|-----------|--------|-----------------|
| Small | 100 × 100 | 10,000 | < 0.1 s |
| Medium | 1000 × 1000 | 1,000,000 | < 1.0 s |

---

## Large Grid (Global 0.5°)

Large grid benchmarks (comparison against `fastbioclim_check`) are **not run in
CI**. They require real-world climate data and significant time/memory. These are
reserved for manual end-of-project validation.

To run manually:

1. Obtain global 0.5° monthly climate rasters (tas, tasmax, tasmin, pr × 12
   months)
2. Use `bioclim_raster()` and compare against `fastbioclim_check::bioclim()`
3. Record wall-clock time, peak memory, and output agreement

---

## Synthetic Data

The benchmark script generates synthetic rasters with realistic climate value
ranges:

- **Temperature (tas):** Seasonal cycle from −5 °C to 25 °C with Gaussian noise
  (sd = 3)
- **Diurnal range:** 5–15 °C (uniform)
- **Max/Min temperature:** Derived from mean ± half the diurnal range
- **Precipitation:** Seasonal pattern 5–80 mm base with Gaussian noise (sd = 20),
  clamped ≥ 0

All random generation uses `set.seed(42)` for reproducibility.

---

## CI Integration

Add the following step to a GitHub Actions workflow to run benchmarks on every
PR:

```yaml
- name: Run benchmarks
  run: Rscript benchmarks/run_benchmarks.R
```

The script will fail the CI step if:

- Small grid (100 × 100) median time > 0.1 s
- Medium grid (1000 × 1000) median time > 1.0 s

---

## Output

The script prints a summary table:

```
=== Summary ===

Benchmark                       Small (s)   Medium (s)   Status
---                                   ---          ---      ---
bioclim_raster()                  0.0XXX       0.XXXX     PASS
bioclim_cpp()                     0.00XX       0.0XXX     PASS
ClimateBlock$compute()            0.00XX       0.0XXX     PASS
Individual bio*_cpp()             0.00XX       0.0XXX     PASS
```

It also reports approximate memory usage for the medium grid allocation.
