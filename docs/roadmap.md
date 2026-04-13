# rxbioclim Optimisation Roadmap

Structured implementation plan with milestones, phases, and checklists.

---

## Phase 0: Immediate Performance Fix (Tier 1) — Match fastbioclim_check

**Milestone: "Baseline Performance"**

### Issue 1: Replace R pixel loop with existing C++ batch call

Replace the R-level `for` loop in `.compute_bioclim_block()` with the compiled
`bioclim_cpp()` function that already exists in `src/bioclim.cpp`.

- [ ] In `R/bioclim_raster.R`, replace `.compute_bioclim_block()` body with a
  call to `bioclim_cpp(v_tas, v_tasmax, v_tasmin, v_pr)`
- [ ] Add NA-row handling at the C++ level: detect rows where any of the four
  input columns contain `NA` and write `NA` to all 19 output columns
- [ ] Remove the R-level `for` loop and per-pixel `bioclim()` call
- [ ] Verify `bioclim_raster()` produces identical output with the new path
- [ ] **Expected: 10–100× speedup** (eliminate R interpreter + S4 dispatch per
  pixel)

### Issue 2: Replace `parallel::parLapply` with OpenMP

Replace R-level PSOCK parallelism with C++-level OpenMP shared-memory
parallelism.

- [ ] Add `-fopenmp` to `src/Makevars` and `src/Makevars.win`
- [ ] Add `#pragma omp parallel for` to the pixel loop in `bioclim_cpp()`
  (`src/bioclim.cpp` line 424)
- [ ] Remove R-level `parallel::makeCluster()` / `parLapply()` code from
  `bioclim_raster()`
- [ ] Add `ncores` parameter to `bioclim_cpp()` and call
  `omp_set_num_threads()` accordingly
- [ ] **Expected: linear speedup with cores** (no serialisation overhead)

---

## Phase 1: Link the Real xbioclim Engine (Tier 2) — Beat fastbioclim_check

**Milestone: "xtensor Integration"**

### Issue 3: Add `configure` / `configure.ac` for GDAL detection

- [ ] Write simplified `configure.ac` (< 100 lines) using `gdal-config
  --cflags` and `gdal-config --libs`
- [ ] Create `src/Makevars.in` template that substitutes `@GDAL_CFLAGS@` and
  `@GDAL_LIBS@`
- [ ] Create `src/Makevars.win` / `src/Makevars.ucrt` for Windows (reference:
  `alrobles/sf` `configure.ac`)
- [ ] Pass `R CMD check --as-cran` on Linux and macOS

### Issue 4: Vendor xtensor headers into `inst/include/`

- [ ] Copy `xtensor` + `xtl` headers (header-only, no compilation needed) into
  `inst/include/xtensor/` and `inst/include/xtl/`
- [ ] Update `src/Makevars.in` to include `-I../inst/include`
- [ ] Verify compilation with `R CMD INSTALL`

> **Depends on:** Issue 3 (Makevars.in must exist)

### Issue 5: Replace `xbioclim.h` with real xbioclim primitives

- [ ] Use `xt::xtensor<float, 2>` instead of `std::vector<double>` in the
  computation core
- [ ] Port `primitives.cpp` (OpenMP-parallelised `rolling_quarter_argmax`,
  etc.) from `alrobles/xbioclim`
- [ ] Port `bioclim.cpp` (vectorised `xt::mean`, `xt::amax` over full arrays)
  from `alrobles/xbioclim`
- [ ] All operations become whole-array vectorised (no per-pixel loop)
- [ ] **Expected: 2–5× over Tier 1**

> **Depends on:** Issue 4 (xtensor headers must be vendored)

### Issue 6: Create zero-copy R↔xtensor bridge

- [ ] Use `xt::adapt(REAL(x), shape, xt::no_ownership())` to map R matrices
  directly onto xtensor views — zero copy on input
- [ ] For output, write directly into pre-allocated R `NumericMatrix` memory
- [ ] Consider `float32` computation (50% memory savings) with a single
  double→float transpose copy at the boundary
- [ ] Benchmark to confirm copy elimination

> **Depends on:** Issue 5 (xtensor types must be in use)

---

## Phase 2: Zero R Memory Engine (Tier 3) — Unreachable by R-based solutions

**Milestone: "Native Engine"**

### Issue 7: Integrate GdalReader/GdalWriter from xbioclim

- [ ] Port `gdal_io.cpp` and `gdal_io.hpp` from `alrobles/xbioclim`
- [ ] Tiled windowed reads — only 1-2 tiles in memory at a time
- [ ] Scale/offset decoding for packed integer rasters
- [ ] COG (Cloud-Optimised GeoTIFF) output support
- [ ] GDAL linked via `configure.ac` from Issue 3

> **Depends on:** Issue 3 (GDAL must be linked)

### Issue 8: Create `BioclimEngine` C++ class

- [ ] XPtr-based engine exposed to R: `open()` → `set_mask()` →
  `set_threads()` → `compute()`
- [ ] R sends file paths + parameters only; C++ does all data handling
- [ ] Tiled processing: read tile → compute → write tile → next
- [ ] Peak memory proportional to tile size, not raster size

> **Depends on:** Issues 5 + 7 (xtensor primitives + GDAL I/O)

### Issue 9: Create R user-facing API `bioclim_engine()`

- [ ] `bioclim_engine(tas_files, tasmax_files, tasmin_files, pr_files,
  output, mask, threads)` — single function call
- [ ] Returns `terra::rast()` pointing at the output file (metadata only, no
  pixel data in R)
- [ ] Document with roxygen2 and add examples

> **Depends on:** Issue 8

### Issue 10: Add polygon masking support

- [ ] Use GDAL to rasterise the polygon mask → binary mask raster
- [ ] Skip masked pixels during tile reads (no wasted computation)
- [ ] Accept both `sf` objects and file paths for the mask argument

> **Depends on:** Issue 7 (GDAL must be available)

---

## Phase 3: Advanced Optimisation

**Milestone: "GPU & CI"**

### Issue 11: Add CUDA backend option (conditional compilation)

- [ ] Detect CUDA toolkit in `configure.ac`
- [ ] Conditionally compile `bioclim_cuda.cu` when CUDA is available
- [ ] Fall back to CPU (OpenMP) path when CUDA is absent
- [ ] Expose `device = "cpu"` / `device = "gpu"` argument in
  `bioclim_engine()`

### Issue 12: Benchmark suite (this issue)

- [ ] `benchmarks/run_benchmarks.R` — synthetic small (100×100) and medium
  (1000×1000) grid benchmarks
- [ ] `benchmarks/README.md` — instructions and target times
- [ ] CI integration: GitHub Actions runs benchmarks on every PR
- [ ] Fail the check if small > 0.1 s or medium > 1 s
- [ ] Large grid (global 0.5°) comparison reserved for manual end-of-project
  testing

### Issue 13: Keep backward compatibility

- [ ] `bioclim_raster()` continues to work with the same signature, just
  faster
- [ ] `bioclim_engine()` is the new high-performance API for advanced users
- [ ] All existing tests pass without modification

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
