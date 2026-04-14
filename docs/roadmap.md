# rxbioclim Optimisation Roadmap

Structured implementation plan with milestones, phases, and checklists.

---

## 🗓️ Milestone Dashboard

> Last updated: 2026-04-14  
> All 13 actionable issues are **closed**. All phases complete. ✅

| Phase | Milestone | Issues | Status |
|-------|-----------|--------|--------|
| **Phase 0** | Baseline Performance | [#16 (Issue 1)][i16], [#17 (Issue 2)][i17] | ✅ Complete |
| **Phase 1** | xtensor Integration | [#18 (Issue 3)][i18], [#19 (Issue 4)][i19], [#20 (Issue 5)][i20], [#21 (Issue 6)][i21] | ✅ Complete |
| **Phase 2** | Native Engine | [#22 (Issue 7)][i22], [#23 (Issue 8)][i23], [#24 (Issue 9)][i24], [#25 (Issue 10)][i25] | ✅ Complete |
| **Phase 3** | GPU & CI | [#26 (Issue 11)][i26], [#27 (Issue 12)][i27], [#28 (Issue 13)][i28] | ✅ Complete |

[i16]: https://github.com/alrobles/rxbioclim/issues/16
[i17]: https://github.com/alrobles/rxbioclim/issues/17
[i18]: https://github.com/alrobles/rxbioclim/issues/18
[i19]: https://github.com/alrobles/rxbioclim/issues/19
[i20]: https://github.com/alrobles/rxbioclim/issues/20
[i21]: https://github.com/alrobles/rxbioclim/issues/21
[i22]: https://github.com/alrobles/rxbioclim/issues/22
[i23]: https://github.com/alrobles/rxbioclim/issues/23
[i24]: https://github.com/alrobles/rxbioclim/issues/24
[i25]: https://github.com/alrobles/rxbioclim/issues/25
[i26]: https://github.com/alrobles/rxbioclim/issues/26
[i27]: https://github.com/alrobles/rxbioclim/issues/27
[i28]: https://github.com/alrobles/rxbioclim/issues/28

---

## Phase 0: Immediate Performance Fix (Tier 1) — Match fastbioclim_check

**Milestone: "Baseline Performance"**

### Issue 1: Replace R pixel loop with existing C++ batch call ✅

Replace the R-level `for` loop in `.compute_bioclim_block()` with the compiled
`bioclim_cpp()` function that already exists in `src/bioclim.cpp`.

- [x] In `R/bioclim_raster.R`, replace `.compute_bioclim_block()` body with a
  call to `bioclim_xt(v_tas, v_tasmax, v_tasmin, v_pr, ncores)` (xtensor-accelerated)
- [x] Add NA-row handling at the C++ level: detect rows where any of the four
  input columns contain `NA` and write `NA` to all 19 output columns
- [x] Remove the R-level `for` loop and per-pixel `bioclim()` call
- [x] Verify `bioclim_raster()` produces identical output with the new path
- [x] **Expected: 10–100× speedup** (eliminate R interpreter + S4 dispatch per
  pixel)

### Issue 2: Replace `parallel::parLapply` with OpenMP ✅

Replace R-level PSOCK parallelism with C++-level OpenMP shared-memory
parallelism.

- [x] Add `-fopenmp` to `src/Makevars` and `src/Makevars.win` (via `$(SHLIB_OPENMP_CXXFLAGS)`)
- [x] Add `#pragma omp parallel for` to the pixel loop in `bioclim_xt()`
  (`src/bioclim_xt.cpp`)
- [x] Remove R-level `parallel::makeCluster()` / `parLapply()` code from
  `bioclim_raster()`
- [x] Add `ncores` parameter to `bioclim_xt()` and call
  `omp_set_num_threads()` accordingly
- [x] **Expected: linear speedup with cores** (no serialisation overhead)

---

## Phase 1: Link the Real xbioclim Engine (Tier 2) — Beat fastbioclim_check

**Milestone: "xtensor Integration"**

### Issue 3: Add `configure` / `configure.ac` for GDAL detection ✅

- [x] Write simplified `configure.ac` (< 100 lines) using `gdal-config
  --cflags` and `gdal-config --libs`
- [x] Create `src/Makevars.in` template that substitutes `@GDAL_CFLAGS@` and
  `@GDAL_LIBS@`
- [x] Create `src/Makevars.win` / `src/Makevars.ucrt` for Windows (reference:
  `alrobles/sf` `configure.ac`)
- [x] Pass `R CMD check --as-cran` on Linux and macOS

### Issue 4: Vendor xtensor headers into `inst/include/` ✅

- [x] Copy `xtensor` + `xtl` headers (header-only, no compilation needed) into
  `inst/include/xtensor/` and `inst/include/xtl/`
- [x] Update `src/Makevars.in` to include `-I../inst/include`
- [x] Verify compilation with `R CMD INSTALL`

> **Depends on:** Issue 3 (Makevars.in must exist)

### Issue 5: Replace `xbioclim.h` with real xbioclim primitives ✅

- [x] Use `xt::xtensor<float, 2>` instead of `std::vector<double>` in the
  computation core
- [x] Port `primitives.cpp` (OpenMP-parallelised `rolling_quarter_argmax`,
  etc.) from `alrobles/xbioclim`
- [x] Port `bioclim.cpp` (vectorised `xt::mean`, `xt::amax` over full arrays)
  from `alrobles/xbioclim`
- [x] All operations become whole-array vectorised (no per-pixel loop)
- [x] **Expected: 2–5× over Tier 1**

> **Depends on:** Issue 4 (xtensor headers must be vendored)

### Issue 6: Create zero-copy R↔xtensor bridge ✅

- [x] Use `xt::adapt(REAL(x), shape, xt::no_ownership())` to map R matrices
  directly onto xtensor views — zero copy on input
- [x] For output, write directly into pre-allocated R `NumericMatrix` memory
- [x] Consider `float32` computation (50% memory savings) with a single
  double→float transpose copy at the boundary
- [x] Benchmark to confirm copy elimination

> **Depends on:** Issue 5 (xtensor types must be in use)

---

## Phase 2: Zero R Memory Engine (Tier 3) — Unreachable by R-based solutions

**Milestone: "Native Engine"**

### Issue 7: Integrate GdalReader/GdalWriter from xbioclim ✅

- [x] Add `src/gdal_io.hpp` and `src/gdal_io.cpp` — GdalReader / GdalWriter
  C++ classes with full `#ifdef HAVE_GDAL` guards
- [x] Tiled windowed reads — only 1–2 tiles in memory at a time
  (`GdalReader::read_window` via GDAL `RasterIO`)
- [x] Scale/offset decoding for packed integer rasters (applied in
  `read_window` when GDAL band metadata is present)
- [x] COG-compatible GTiff output — `GdalWriter(cog_compatible = true)` sets
  `TILED=YES COMPRESS=LZW BLOCKXSIZE=256 BLOCKYSIZE=256 BIGTIFF=IF_SAFER`
- [x] Optional GDAL — all code compiles without GDAL; functions stop with a
  clear message at runtime when GDAL is absent
- [x] Rcpp-exported diagnostics: `gdal_can_open(path)` and `gdal_info(path)`
- [x] Unit tests in `tests/testthat/test-gdal-io.R` that skip when GDAL is
  absent; tiny test raster in `inst/extdata/tiny.tif`

> **Depends on:** Issue 3 (GDAL must be linked)

### Issue 8: Create `BioclimEngine` C++ class ✅

- [x] XPtr-based engine exposed to R: `open()` → `set_mask()` →
  `set_threads()` → `compute()`
- [x] R sends file paths + parameters only; C++ does all data handling
- [x] Tiled processing: read tile → compute → write tile → next
- [x] Peak memory proportional to tile size, not raster size

> **Depends on:** Issues 5 + 7 (xtensor primitives + GDAL I/O)

### Issue 9: Create R user-facing API `bioclim_engine()` ✅

- [x] `bioclim_engine(tas, tasmax, tasmin, pr, output, mask, threads,
  tile_size, overwrite, device)` — single function call
- [x] Returns `terra::rast()` pointing at the output file (metadata only, no
  pixel data in R)
- [x] Document with roxygen2 and add examples

> **Depends on:** Issue 8

### Issue 10: Add polygon masking support ✅

- [x] Use GDAL to rasterise the polygon mask → binary mask raster
- [x] Skip masked pixels during tile reads (no wasted computation)
- [x] Accept both `sf` objects and file paths for the mask argument

> **Depends on:** Issue 7 (GDAL must be available)

---

## Phase 3: Advanced Optimisation

**Milestone: "GPU & CI"**

### Issue 11: Add CUDA backend option (conditional compilation) ✅

- [x] Detect CUDA toolkit in `configure.ac`
- [x] Conditionally compile `bioclim_cuda.cu` when CUDA is available
- [x] Fall back to CPU (OpenMP) path when CUDA is absent
- [x] Expose `device = "auto"` / `device = "cpu"` / `device = "gpu"` argument in
  `bioclim_engine()`

### Issue 12: Benchmark suite ✅

- [x] `benchmarks/run_benchmarks.R` — synthetic small (100×100) and medium
  (1000×1000) grid benchmarks
- [x] `benchmarks/README.md` — instructions and target times
- [x] CI integration: GitHub Actions runs benchmarks on every PR
- [x] Fail the check if small > 0.1 s or medium > 1 s
- [x] Large grid (global 0.5°) comparison reserved for manual end-of-project
  testing

### Issue 13: Keep backward compatibility ✅

- [x] `bioclim_raster()` continues to work with the same signature, just
  faster
- [x] `bioclim_engine()` is the new high-performance API for advanced users
- [x] All existing tests pass without modification

---

## Dependency Graph

```
Issue 1 (R loop → C++ batch)     ← standalone, do first
Issue 2 (OpenMP)                 ← standalone, can parallel with Issue 1
Issue 3 (configure.ac / GDAL)   ← standalone
Issue 4 (vendor xtensor)         ← depends on Issue 3
Issue 5 (xtensor primitives)     ← depends on Issue 4
Issue 6 (zero-copy bridge)       ← depends on Issue 5
Issue 7 (GDAL I/O)              ← depends on Issue 3
Issue 8 (BioclimEngine)          ← depends on Issues 5 + 7
Issue 9 (R API bioclim_engine)   ← depends on Issue 8
Issue 10 (polygon masking)       ← depends on Issue 7
Issue 11 (CUDA)                  ← depends on Issue 8
Issue 12 (benchmarks)            ← standalone (this issue)
Issue 13 (backward compat)       ← ongoing throughout
```

**Recommended execution order:**

1. Issues 1 + 2 + 12 in parallel (immediate wins + benchmarks)
2. Issue 3 (build system)
3. Issues 4 → 5 → 6 (xtensor integration, sequential)
4. Issue 7 (GDAL I/O, can start after Issue 3)
5. Issues 8 → 9 (engine, sequential)
6. Issues 10, 11 (advanced, can be deferred)
7. Issue 13 (ongoing validation)

---

## CI Integration

A GitHub Actions workflow should run the benchmark suite on every PR:

- **Runs:** `Rscript benchmarks/run_benchmarks.R`
- **Grids tested:** small (100×100) and medium (1000×1000) only
- **Failure condition:** exit code 1 if small > 0.1 s or medium > 1 s
- **No large grid in CI** — large grid benchmarks use real-world data and are
  too slow / data-heavy for CI runners
