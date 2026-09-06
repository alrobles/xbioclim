# xbioclim: CRAN Submission Roadmap

> Last updated: 2026-04-14
> Status: **In Progress** — All prior optimization phases complete. CRAN submission prep begins now.

---

## Overview

This roadmap covers the engineering work needed to get `xbioclim` accepted on
CRAN while preserving its optional CUDA GPU acceleration. The central challenge
is that **CRAN's check farm has no GPUs** — the package must build, install, and
pass `R CMD check --as-cran` cleanly on CPU-only machines across Linux, macOS,
and Windows. GPU functionality must be **strictly optional**.

**Strategy:** Compile-time optional CUDA (precedented by [MPCR](https://cran.r-project.org/package=MPCR)
on CRAN). The existing `#ifdef HAVE_CUDA` guard architecture is already ~80%
CRAN-ready.

---

## Phase 0: Fix CRAN Blockers (P0)

These issues will cause **immediate rejection** if not addressed.

| # | Issue | File(s) affected | Tracking |
|---|-------|-------------------|----------|
| 1 | Rename `.cu` → `.cpp` + conditional nvcc rule | `src/bioclim_cuda.cu`, `src/Makevars.in` | [#29][i29] |
| 2 | Remove GNU Make extensions from Makevars | `src/Makevars.in`, `DESCRIPTION` | [#30][i30] |
| 3 | Add `SystemRequirements` and `NeedsCompilation` to DESCRIPTION | `DESCRIPTION` | [#31][i31] |
| 4 | Drop `CXX_STD = CXX17` (default since R 4.3) | `src/Makevars.in`, `src/Makevars.win`, `src/Makevars.ucrt` | [#32][i32] |
| 5 | Remove generated `src/Makevars` from git | `src/Makevars`, `.gitignore` | [#33][i33] |

[i29]: https://github.com/alrobles/xbioclim/issues/29
[i30]: https://github.com/alrobles/xbioclim/issues/30
[i31]: https://github.com/alrobles/xbioclim/issues/31
[i32]: https://github.com/alrobles/xbioclim/issues/32
[i33]: https://github.com/alrobles/xbioclim/issues/33

---

## Phase 1: Harden Build System (P1)

These will likely cause rejection during CRAN human review.

| # | Issue | File(s) affected | Tracking |
|---|-------|-------------------|----------|
| 6 | Expand `cleanup` script | `cleanup` | [#34][i34] |
| 7 | Fix `tools/` directory conflict with `.Rbuildignore` | `configure.ac`, `.Rbuildignore` | [#35][i35] |
| 8 | Add `--without-cuda` configure flag | `configure.ac` | [#36][i36] |
| 9 | Add nvcc compile test in configure | `configure.ac` | [#36][i36] |
| 10 | Cap OpenMP threads for CRAN compliance | `src/BioclimEngine.cpp`, `src/bioclim.cpp`, `src/bioclim_xt.cpp` | [#37][i37] |

[i34]: https://github.com/alrobles/xbioclim/issues/34
[i35]: https://github.com/alrobles/xbioclim/issues/35
[i36]: https://github.com/alrobles/xbioclim/issues/36
[i37]: https://github.com/alrobles/xbioclim/issues/37

---

## Phase 2: Polish and Pre-Submission (P2)

| # | Issue | File(s) affected | Tracking |
|---|-------|-------------------|----------|
| 11 | Add CUDA architecture flags for portability | `configure.ac` | [#36][i36] |
| 12 | Wrap engine examples in `\donttest{}` | `R/bioclim_engine.R` | [#38][i38] |
| 13 | Write `cran-comments.md` + test on R-hub / win-builder | `cran-comments.md` | [#39][i39] |

[i38]: https://github.com/alrobles/xbioclim/issues/38
[i39]: https://github.com/alrobles/xbioclim/issues/39

---

## CRAN Policy Quick Reference

| Requirement | Impact on xbioclim |
|-------------|---------------------|
| `configure` must be POSIX `/bin/sh` | Current `configure.ac` is POSIX-clean ✓ |
| No GNU Make extensions without declaration | **VIOLATION**: `$(wildcard)` / `$(patsubst)` in `Makevars.in` |
| `R CMD check` no ERRORs or WARNINGs | `.cu` extension triggers WARNING |
| Max 2 threads during checks | `omp_set_num_threads(ncores)` uncapped |
| Examples run in a few seconds | Engine examples need `\donttest{}` |
| Source for all components provided | CUDA `.cu` sources included ✓ |

---

## Cross-Platform Compilation Matrix

| Platform | C++ | CUDA/nvcc | GDAL | OpenMP | Expected |
|----------|-----|-----------|------|--------|----------|
| Linux x86_64 (Debian/Fedora) | gcc/clang | ✗ on check farm | ✓ | ✓ | CPU-only build |
| macOS arm64 (Apple Silicon) | clang | ✗ (no NVIDIA) | ✓ | Limited | CPU-only, always |
| macOS x86_64 | clang | ✗ | ✓ | Limited | CPU-only |
| Windows x86_64 (Rtools) | gcc | ✗ (no nvcc) | ✓ | ✓ | CPU-only via `Makevars.win` |

---

## CRAN Precedent for GPU Packages

| Package | Strategy | Status |
|---------|----------|--------|
| **MPCR** | CMake + optional CUDA at configure time | ✓ Active on CRAN (2026) |
| **torch** | Runtime download of GPU libs | ✓ Active |
| **gpuR** | OpenCL compile-time detection | ✓ Active |
| **gputools** | Hard nvcc requirement | ✗ Archived (2017) |

**xbioclim follows the MPCR approach** — the safest compile-time strategy with
proven CRAN acceptance.

---

## Testing Matrix Before Submission

| Scenario | Expected | How to test |
|----------|----------|-------------|
| Linux, no CUDA, no GDAL | Installs; CPU tests pass; `has_cuda()` → FALSE | `--configure-args="--without-cuda"` |
| Linux, no CUDA, with GDAL | Engine works in CPU mode | Standard install |
| Linux, with CUDA + GDAL | Full GPU acceleration | GPU machine |
| macOS arm64 | CPU-only; all tests pass | R-hub / mac-builder |
| Windows Rtools | CPU-only via `Makevars.win` | win-builder |
| `R CMD check --as-cran` | 0 errors, 0 warnings, ≤1 note | Local |

---

## Estimated Timeline

| Phase | Effort | Calendar |
|-------|--------|----------|
| Phase 0 (blockers) | 2–3 days | Week 1 |
| Phase 1 (build hardening) | 1–2 days | Week 1–2 |
| Phase 2 (polish) | 2–3 days | Week 2 |
| R-hub / win-builder testing | 2–3 days wait | Week 2–3 |
| CRAN submission + reviewer rounds | 3–5 days | Week 3 |

**Total: ~3 weeks to CRAN acceptance.**

---

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| `.cu` → `.cpp` rename rejected by reviewer | Low | High | Backup: move `.cu` to `inst/cuda/` + runtime `dlopen()` |
| nvcc/g++ flag incompatibility | Medium | Medium | `--without-cuda` escape hatch |
| CRAN compiler update breaks CUDA | Medium | Low | GPU is optional; CPU always works |
| macOS OpenMP issues | Medium | Medium | Test on mac-builder first |
| Reviewer requests removal of all CUDA | Low | High | Cite MPCR, gpuR precedent |

---

## Dependency Graph

```
Issue #29 (rename .cu → .cpp)        ← P0, do first
Issue #30 (GNU Make extensions)       ← P0, depends on #29 (Makevars rewrite)
Issue #31 (DESCRIPTION fields)        ← P0, standalone
Issue #32 (drop CXX_STD)             ← P0, standalone
Issue #33 (remove src/Makevars)       ← P0, standalone
Issue #34 (cleanup script)            ← P1, standalone
Issue #35 (tools/ directory)          ← P1, standalone
Issue #36 (configure hardening)       ← P1, depends on #29
Issue #37 (OpenMP thread cap)         ← P1, standalone
Issue #38 (example guards)            ← P2, standalone
Issue #39 (cran-comments + testing)   ← P2, depends on all above
```

**Execution order:**
1. #29, #31, #32, #33 in parallel (independent P0 fixes)
2. #30 after #29 (Makevars rewrite depends on .cu rename)
3. #34, #35, #37 in parallel (independent P1 fixes)
4. #36 after #29 (configure hardening includes nvcc compile test for renamed file)
5. #38 (example guards)
6. #39 last (submission prep depends on everything)
