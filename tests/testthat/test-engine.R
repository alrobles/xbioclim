# Tests for BioclimEngine (Issue #23) — XPtr-based tiled computation engine.
#
# Test structure:
#   1. Engine creation — works on all platforms (no GDAL needed).
#   2. Setter methods — work without GDAL.
#   3. compute() without GDAL — stops with informative error.
#   4. Full round-trip with GDAL — skipped when GDAL is unavailable.

# ── Helpers ──────────────────────────────────────────────────────────────────

skip_without_gdal <- function() {
  if (!has_gdal()) testthat::skip("Package built without GDAL support")
}

skip_without_terra <- function() {
  testthat::skip_if_not_installed("terra")
}

# Create a tiny (3×3, 12-band) GeoTIFF in a temp directory filled with a
# constant value.  Returns the file path.
make_tiny_raster <- function(value, path, nbands = 12L) {
  r <- terra::rast(
    nrows = 3L, ncols = 3L, nlyrs = nbands,
    xmin = 0, xmax = 1, ymin = 0, ymax = 1,
    crs = "EPSG:4326"
  )
  terra::values(r) <- value
  terra::writeRaster(r, path, overwrite = TRUE)
  path
}

# Build a list of 12 single-band file paths for one climate variable, all
# filled with `value`, stored under `dir`.
make_monthly_files <- function(value, dir, prefix) {
  paths <- character(12L)
  for (m in seq_len(12L)) {
    paths[m] <- file.path(dir, sprintf("%s_%02d.tif", prefix, m))
    make_tiny_raster(value, paths[m], nbands = 1L)
  }
  paths
}

# ── 1. Engine creation ────────────────────────────────────────────────────────

test_that("engine_create() returns an externalptr", {
  ptr <- engine_create()
  expect_true(is(ptr, "externalptr"))
})

test_that("engine_create() produces a non-null pointer", {
  ptr <- engine_create()
  expect_false(is.null(ptr))
})

# ── 2. Setter methods ─────────────────────────────────────────────────────────

test_that("engine setter methods run without error", {
  ptr <- engine_create()
  expect_no_error(engine_set_output(ptr, tempfile(fileext = ".tif")))
  expect_no_error(engine_set_mask(ptr, ""))
  expect_no_error(engine_set_threads(ptr, 2L))
  expect_no_error(engine_set_tile_size(ptr, 64L))
})

test_that("engine_set_threads clamps values < 1 to 1 without error", {
  ptr <- engine_create()
  expect_no_error(engine_set_threads(ptr, 0L))
  expect_no_error(engine_set_threads(ptr, -5L))
})

test_that("engine_set_tile_size clamps values < 1 to 1 without error", {
  ptr <- engine_create()
  expect_no_error(engine_set_tile_size(ptr, 0L))
})

# ── 3. compute() without GDAL ─────────────────────────────────────────────────

test_that("engine_compute stops with clear message when GDAL absent", {
  if (has_gdal()) testthat::skip("GDAL is present — testing absence path N/A")

  ptr <- engine_create()
  expect_error(
    engine_compute(ptr),
    regexp = "GDAL is required for BioclimEngine"
  )
})

# ── 4. has_gdal() ─────────────────────────────────────────────────────────────

test_that("has_gdal() returns a logical scalar", {
  result <- has_gdal()
  expect_type(result, "logical")
  expect_length(result, 1L)
})

# ── 5. Full round-trip (GDAL required) ───────────────────────────────────────

test_that("engine_compute produces a 19-band GeoTIFF (GDAL + terra)", {
  skip_without_gdal()
  skip_without_terra()

  tmpdir <- tempfile("engine_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  # Monthly temperature files (12 files, each 1 band, 3×3 pixels).
  # Values chosen to produce valid (non-NA) bioclimatic variables:
  #   tas    ≈ realistic monthly mean temps (°C × 10)
  #   tasmax > tas, tasmin < tas, pr ≥ 0
  tas_vals    <- c(5, 7, 10, 14, 18, 22, 25, 24, 20, 15, 10,  6)
  tasmax_vals <- c(8, 10, 14, 18, 23, 28, 32, 31, 26, 19, 13,  9)
  tasmin_vals <- c(1,  3,  6, 10, 13, 17, 20, 19, 15, 10,  6,  2)
  pr_vals     <- c(60, 55, 48, 35, 28, 22, 18, 20, 35, 55, 65, 68)

  tas_files    <- make_monthly_files(tas_vals,    tmpdir, "tas")
  tasmax_files <- make_monthly_files(tasmax_vals, tmpdir, "tasmax")
  tasmin_files <- make_monthly_files(tasmin_vals, tmpdir, "tasmin")
  pr_files     <- make_monthly_files(pr_vals,     tmpdir, "pr")

  # But each make_monthly_files creates files with a single repeated value.
  # Re-create them so month m has the correct scalar value.
  for (m in seq_len(12L)) {
    make_tiny_raster(tas_vals[m],    tas_files[m],    nbands = 1L)
    make_tiny_raster(tasmax_vals[m], tasmax_files[m], nbands = 1L)
    make_tiny_raster(tasmin_vals[m], tasmin_files[m], nbands = 1L)
    make_tiny_raster(pr_vals[m],     pr_files[m],     nbands = 1L)
  }

  output_path <- file.path(tmpdir, "bioclim_output.tif")

  ptr <- engine_create()
  engine_open(ptr, tas_files, tasmax_files, tasmin_files, pr_files)
  engine_set_output(ptr, output_path)
  engine_set_tile_size(ptr, 2L)   # tiny tiles to exercise edge-tile code
  engine_set_threads(ptr, 1L)

  result_path <- engine_compute(ptr)

  # Output file must exist and match the configured path.
  expect_equal(result_path, output_path)
  expect_true(file.exists(output_path))

  # Output must have exactly 19 bands.
  out_rast <- terra::rast(output_path)
  expect_equal(terra::nlyr(out_rast), 19L)

  # Output dimensions must match input (3 rows × 3 cols).
  expect_equal(terra::nrow(out_rast), 3L)
  expect_equal(terra::ncol(out_rast), 3L)

  # BIO01 (mean annual temperature) should be finite for all pixels.
  bio01_vals <- terra::values(out_rast[[1L]])
  expect_true(all(is.finite(bio01_vals)),
              info = "BIO01 should be finite for valid inputs")

  # BIO12 (annual precipitation) should equal sum of monthly pr values.
  expected_pr_sum <- sum(pr_vals)
  bio12_vals <- terra::values(out_rast[[12L]])
  expect_true(all(abs(bio12_vals - expected_pr_sum) < 1e-6),
              info = "BIO12 should equal the sum of monthly precipitation")
})

test_that("engine_compute handles 1-multi-band-file input (GDAL + terra)", {
  skip_without_gdal()
  skip_without_terra()

  tmpdir <- tempfile("engine_multiband_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  tas_vals    <- c(5, 7, 10, 14, 18, 22, 25, 24, 20, 15, 10,  6)
  tasmax_vals <- c(8, 10, 14, 18, 23, 28, 32, 31, 26, 19, 13,  9)
  tasmin_vals <- c(1,  3,  6, 10, 13, 17, 20, 19, 15, 10,  6,  2)
  pr_vals     <- c(60, 55, 48, 35, 28, 22, 18, 20, 35, 55, 65, 68)

  # Create one 12-band file per variable.
  make_multiband <- function(vals, path) {
    r <- terra::rast(
      nrows = 3L, ncols = 3L, nlyrs = 12L,
      xmin = 0, xmax = 1, ymin = 0, ymax = 1,
      crs = "EPSG:4326"
    )
    for (m in seq_len(12L))
      terra::values(r[[m]]) <- vals[m]
    terra::writeRaster(r, path, overwrite = TRUE)
    path
  }

  tas_file    <- make_multiband(tas_vals,    file.path(tmpdir, "tas.tif"))
  tasmax_file <- make_multiband(tasmax_vals, file.path(tmpdir, "tasmax.tif"))
  tasmin_file <- make_multiband(tasmin_vals, file.path(tmpdir, "tasmin.tif"))
  pr_file     <- make_multiband(pr_vals,     file.path(tmpdir, "pr.tif"))

  output_path <- file.path(tmpdir, "bioclim_mb.tif")

  ptr <- engine_create()
  engine_open(ptr, tas_file, tasmax_file, tasmin_file, pr_file)
  engine_set_output(ptr, output_path)
  engine_set_tile_size(ptr, 2L)

  result_path <- engine_compute(ptr)

  expect_true(file.exists(result_path))
  out_rast <- terra::rast(result_path)
  expect_equal(terra::nlyr(out_rast), 19L)
})

test_that("engine_compute stops when output path not set", {
  skip_without_gdal()

  ptr <- engine_create()
  # open() with dummy paths — will fail at compute() because output is missing
  # before files are even checked
  expect_error(
    engine_compute(ptr),
    regexp = "output path not set"
  )
})

test_that("engine_compute stops when file list has wrong length (GDAL)", {
  skip_without_gdal()

  ptr <- engine_create()
  engine_set_output(ptr, tempfile(fileext = ".tif"))

  # 3 files — neither 1 nor 12
  bad_files <- rep(tempfile(fileext = ".tif"), 3L)
  engine_open(ptr, bad_files, bad_files, bad_files, bad_files)

  expect_error(
    engine_compute(ptr),
    regexp = "1 multi-band file or 12 single-band files"
  )
})
