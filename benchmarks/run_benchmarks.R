#!/usr/bin/env Rscript
#
# benchmarks/run_benchmarks.R
#
# Standalone benchmark script for rxbioclim.
#
# Usage:
#   Rscript benchmarks/run_benchmarks.R
#
# Exit codes:
#   0 — all targets met
#   1 — at least one target exceeded
#
# Targets (median wall-clock time):
#   Small  grid (100 x 100):   < 0.1 s
#   Medium grid (1000 x 1000): < 1.0 s
#
# NOTE: Large grid (global 0.5 deg) comparisons against fastbioclim_check are
#       done manually at the end of the project, NOT in CI.

# ── Dependencies ──────────────────────────────────────────────────────────────

stopifnot(
  requireNamespace("terra",     quietly = TRUE),
  requireNamespace("rxbioclim", quietly = TRUE)
)

library(rxbioclim)

# ── Helpers ───────────────────────────────────────────────────────────────────

#' Create synthetic monthly climate SpatRasters with realistic values
#'
#' @param nrow,ncol Grid dimensions (pixels).
#' @return A named list with elements `tas`, `tasmax`, `tasmin`, `pr` — each a
#'   SpatRaster with 12 layers.
make_climate_rasters <- function(nrow, ncol) {
  n <- nrow * ncol
  set.seed(42L)

  # Monthly mean temperatures: seasonal cycle -5 to 25 deg C with noise
  base_cycle <- seq(-5, 25, length.out = 7)
  base_cycle <- c(base_cycle, rev(base_cycle[-c(1, 7)]))[1:12]  # symmetric
  tas_vals <- matrix(
    rep(base_cycle, each = n) + stats::rnorm(n * 12L, sd = 3),
    nrow = n, ncol = 12L
  )

  # Diurnal range: 5-15 deg C
  diurnal <- matrix(
    stats::runif(n * 12L, min = 5, max = 15),
    nrow = n, ncol = 12L
  )

  tasmax_vals <- tas_vals + diurnal / 2
  tasmin_vals <- tas_vals - diurnal / 2

  # Precipitation: 0-500 mm with seasonal pattern
  pr_base <- c(80, 70, 60, 40, 25, 10, 5, 10, 30, 50, 65, 80)
  pr_vals <- matrix(
    pmax(0, rep(pr_base, each = n) + stats::rnorm(n * 12L, sd = 20)),
    nrow = n, ncol = 12L
  )

  make_rast <- function(vals, nr, nc) {
    r <- terra::rast(nrows = nr, ncols = nc, nlyr = 12L)
    terra::values(r) <- vals
    r
  }

  list(
    tas    = make_rast(tas_vals,    nrow, ncol),
    tasmax = make_rast(tasmax_vals, nrow, ncol),
    tasmin = make_rast(tasmin_vals, nrow, ncol),
    pr     = make_rast(pr_vals,    nrow, ncol)
  )
}

#' Time an expression, returning median seconds over `n` iterations
#'
#' @param expr An expression to benchmark (quoted or brace-wrapped).
#' @param n Number of iterations.
#' @return Median elapsed time in seconds.
bench_median <- function(expr, n = 3L) {
  timings <- numeric(n)
  for (i in seq_len(n)) {
    gc(verbose = FALSE)
    t0 <- proc.time()[["elapsed"]]
    eval(expr, envir = parent.frame())
    timings[i] <- proc.time()[["elapsed"]] - t0
  }
  stats::median(timings)
}

# ── Benchmark: bioclim_cpp() batch path (matrix interface) ────────────────────

bench_bioclim_cpp <- function(nrow, ncol, n_iter = 3L) {
  n <- nrow * ncol
  set.seed(42L)

  base_cycle <- seq(-5, 25, length.out = 7)
  base_cycle <- c(base_cycle, rev(base_cycle[-c(1, 7)]))[1:12]
  tas <- matrix(
    rep(base_cycle, each = n) + stats::rnorm(n * 12L, sd = 3),
    nrow = n, ncol = 12L
  )
  diurnal <- matrix(stats::runif(n * 12L, min = 5, max = 15), nrow = n, ncol = 12L)
  tasmax <- tas + diurnal / 2
  tasmin <- tas - diurnal / 2
  pr_base <- c(80, 70, 60, 40, 25, 10, 5, 10, 30, 50, 65, 80)
  pr <- matrix(
    pmax(0, rep(pr_base, each = n) + stats::rnorm(n * 12L, sd = 20)),
    nrow = n, ncol = 12L
  )

  bench_median(quote(bioclim_cpp(tas, tasmax, tasmin, pr)), n = n_iter)
}

# ── Benchmark: ClimateBlock$compute() (Rcpp module path) ──────────────────────

bench_climate_block <- function(nrow, ncol, n_iter = 3L) {
  n <- nrow * ncol
  set.seed(42L)

  base_cycle <- seq(-5, 25, length.out = 7)
  base_cycle <- c(base_cycle, rev(base_cycle[-c(1, 7)]))[1:12]
  tas <- matrix(
    rep(base_cycle, each = n) + stats::rnorm(n * 12L, sd = 3),
    nrow = n, ncol = 12L
  )
  diurnal <- matrix(stats::runif(n * 12L, min = 5, max = 15), nrow = n, ncol = 12L)
  tasmax <- tas + diurnal / 2
  tasmin <- tas - diurnal / 2
  pr_base <- c(80, 70, 60, 40, 25, 10, 5, 10, 30, 50, 65, 80)
  pr <- matrix(
    pmax(0, rep(pr_base, each = n) + stats::rnorm(n * 12L, sd = 20)),
    nrow = n, ncol = 12L
  )

  bench_median(quote({
    block <- new(ClimateBlock, tas, tasmax, tasmin, pr)
    block$compute()
  }), n = n_iter)
}

# ── Benchmark: individual bio*_cpp() functions ────────────────────────────────

bench_individual_bio <- function(nrow, ncol, n_iter = 3L) {
  n <- nrow * ncol
  set.seed(42L)

  base_cycle <- seq(-5, 25, length.out = 7)
  base_cycle <- c(base_cycle, rev(base_cycle[-c(1, 7)]))[1:12]
  tas <- matrix(
    rep(base_cycle, each = n) + stats::rnorm(n * 12L, sd = 3),
    nrow = n, ncol = 12L
  )
  diurnal <- matrix(stats::runif(n * 12L, min = 5, max = 15), nrow = n, ncol = 12L)
  tasmax <- tas + diurnal / 2
  tasmin <- tas - diurnal / 2
  pr_base <- c(80, 70, 60, 40, 25, 10, 5, 10, 30, 50, 65, 80)
  pr <- matrix(
    pmax(0, rep(pr_base, each = n) + stats::rnorm(n * 12L, sd = 20)),
    nrow = n, ncol = 12L
  )

  bench_median(quote({
    bio01_cpp(tas)
    bio02_cpp(tasmax, tasmin)
    bio03_cpp(tasmax, tasmin)
    bio04_cpp(tas)
    bio05_cpp(tasmax)
    bio06_cpp(tasmin)
    bio07_cpp(tasmax, tasmin)
    bio08_cpp(tas, pr)
    bio09_cpp(tas, pr)
    bio10_cpp(tas)
    bio11_cpp(tas)
    bio12_cpp(pr)
    bio13_cpp(pr)
    bio14_cpp(pr)
    bio15_cpp(pr)
    bio16_cpp(pr)
    bio17_cpp(pr)
    bio18_cpp(tas, pr)
    bio19_cpp(tas, pr)
  }), n = n_iter)
}

# ── Benchmark: bioclim_raster() full pipeline ─────────────────────────────────

bench_bioclim_raster <- function(nrow, ncol, n_iter = 3L) {
  clim <- make_climate_rasters(nrow, ncol)
  bench_median(quote(
    bioclim_raster(clim$tas, clim$tasmax, clim$tasmin, clim$pr)
  ), n = n_iter)
}

# ── Run benchmarks ────────────────────────────────────────────────────────────

cat("=== rxbioclim benchmark suite ===\n\n")

# Small grid: 100 x 100 = 10,000 pixels
cat("--- Small grid (100 x 100 = 10,000 pixels) ---\n")
small_raster  <- bench_bioclim_raster(100L, 100L)
small_cpp     <- bench_bioclim_cpp(100L, 100L)
small_block   <- bench_climate_block(100L, 100L)
small_indiv   <- bench_individual_bio(100L, 100L)

cat(sprintf("  bioclim_raster():        %.4f s\n", small_raster))
cat(sprintf("  bioclim_cpp():           %.4f s\n", small_cpp))
cat(sprintf("  ClimateBlock$compute():  %.4f s\n", small_block))
cat(sprintf("  Individual bio*_cpp():   %.4f s\n", small_indiv))
cat(sprintf("  Target:                  < 0.1 s\n\n"))

# Medium grid: 1000 x 1000 = 1,000,000 pixels
cat("--- Medium grid (1000 x 1000 = 1,000,000 pixels) ---\n")
medium_raster <- bench_bioclim_raster(1000L, 1000L)
medium_cpp    <- bench_bioclim_cpp(1000L, 1000L)
medium_block  <- bench_climate_block(1000L, 1000L)
medium_indiv  <- bench_individual_bio(1000L, 1000L)

cat(sprintf("  bioclim_raster():        %.4f s\n", medium_raster))
cat(sprintf("  bioclim_cpp():           %.4f s\n", medium_cpp))
cat(sprintf("  ClimateBlock$compute():  %.4f s\n", medium_block))
cat(sprintf("  Individual bio*_cpp():   %.4f s\n", medium_indiv))
cat(sprintf("  Target:                  < 1.0 s\n\n"))

# ── Summary table ─────────────────────────────────────────────────────────────

cat("=== Summary ===\n\n")
cat(sprintf("%-28s %12s %12s %8s\n", "Benchmark", "Small (s)", "Medium (s)", "Status"))
cat(sprintf("%-28s %12s %12s %8s\n", "---", "---", "---", "---"))

report_row <- function(label, small_t, medium_t, small_lim, medium_lim) {
  s_ok <- small_t < small_lim
  m_ok <- medium_t < medium_lim
  status <- if (s_ok && m_ok) "PASS" else "FAIL"
  cat(sprintf("%-28s %12.4f %12.4f %8s\n", label, small_t, medium_t, status))
  s_ok && m_ok
}

pass <- TRUE
pass <- report_row("bioclim_raster()",       small_raster, medium_raster, 0.1, 1.0) && pass
pass <- report_row("bioclim_cpp()",          small_cpp,    medium_cpp,    0.1, 1.0) && pass
pass <- report_row("ClimateBlock$compute()", small_block,  medium_block,  0.1, 1.0) && pass
pass <- report_row("Individual bio*_cpp()",  small_indiv,  medium_indiv,  0.1, 1.0) && pass

cat("\n")

# ── Memory usage snapshot ─────────────────────────────────────────────────────

cat("=== Memory usage (approximate) ===\n\n")
gc_before <- gc(verbose = FALSE)
clim_medium <- make_climate_rasters(1000L, 1000L)
gc_after  <- gc(verbose = FALSE)
mem_mb <- sum(gc_after[, 2] - gc_before[, 2])
cat(sprintf("  Raster creation (1000x1000, 4 vars x 12 months): ~%.1f MB\n", mem_mb))
rm(clim_medium)
gc(verbose = FALSE)

cat("\n")

# ── Large grid note ───────────────────────────────────────────────────────────

cat("NOTE: Large grid (global 0.5 deg) comparison against fastbioclim_check\n")
cat("      is done manually at the end of the project, NOT in CI.\n")
cat("      See docs/performance-architecture.md for details.\n\n")

# ── Exit code ─────────────────────────────────────────────────────────────────

if (pass) {
  cat("All benchmark targets met.\n")
  quit(save = "no", status = 0L)
} else {
  cat("WARNING: One or more benchmark targets exceeded.\n")
  quit(save = "no", status = 1L)
}
