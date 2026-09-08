# Tests for bioclim_engine() — R user-facing API
#
# Test structure:
#   1. Input validation (no GDAL needed for most checks).
#   2. GDAL-gated engine round-trip tests.
#   3. Variable selection tests.
#   4. Mask parameter tests.
#   5. Overwrite behavior.

# ── Helpers ──────────────────────────────────────────────────────────────────

skip_without_gdal <- function() {
  if (!has_gdal()) testthat::skip("Package built without GDAL support")
}

skip_without_terra <- function() {
  testthat::skip_if_not_installed("terra")
}

# Create a tiny (3×3, nbands-layer) GeoTIFF filled with a constant value.
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

# Create a list of 12 single-band file paths for one climate variable.
make_monthly_files <- function(vals, dir, prefix) {
  paths <- character(12L)
  for (m in seq_len(12L)) {
    paths[m] <- file.path(dir, sprintf("%s_%02d.tif", prefix, m))
    make_tiny_raster(vals[m], paths[m], nbands = 1L)
  }
  paths
}

# Create one 12-band multi-band GeoTIFF for one climate variable.
make_multiband_file <- function(vals, path) {
  r <- terra::rast(
    nrows = 3L, ncols = 3L, nlyrs = 12L,
    xmin = 0, xmax = 1, ymin = 0, ymax = 1,
    crs = "EPSG:4326"
  )
  for (m in seq_len(12L)) terra::values(r[[m]]) <- vals[m]
  terra::writeRaster(r, path, overwrite = TRUE)
  path
}

# Standard monthly climate values for round-trip tests.
std_tas_vals    <- c(5, 7, 10, 14, 18, 22, 25, 24, 20, 15, 10,  6)
std_tasmax_vals <- c(8, 10, 14, 18, 23, 28, 32, 31, 26, 19, 13,  9)
std_tasmin_vals <- c(1,  3,  6, 10, 13, 17, 20, 19, 15, 10,  6,  2)
std_pr_vals     <- c(60, 55, 48, 35, 28, 22, 18, 20, 35, 55, 65, 68)

# ── 1. Input validation ───────────────────────────────────────────────────────

test_that("bioclim_engine() stops when GDAL is absent", {
  if (has_gdal()) testthat::skip("GDAL is present — absence path N/A")
  skip_without_terra()

  tmpdir <- tempfile("be_val_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  tas_f    <- make_monthly_files(std_tas_vals,    tmpdir, "tas")
  tasmax_f <- make_monthly_files(std_tasmax_vals, tmpdir, "tasmax")
  tasmin_f <- make_monthly_files(std_tasmin_vals, tmpdir, "tasmin")
  pr_f     <- make_monthly_files(std_pr_vals,     tmpdir, "pr")

  expect_error(
    bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f),
    regexp = "GDAL"
  )
})

test_that("bioclim_engine() stops for invalid 'output' argument", {
  skip_without_gdal()
  skip_without_terra()

  tmpdir <- tempfile("be_val_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  tas_f    <- make_monthly_files(std_tas_vals,    tmpdir, "tas")
  tasmax_f <- make_monthly_files(std_tasmax_vals, tmpdir, "tasmax")
  tasmin_f <- make_monthly_files(std_tasmin_vals, tmpdir, "tasmin")
  pr_f     <- make_monthly_files(std_pr_vals,     tmpdir, "pr")

  expect_error(
    bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f, output = 123L),
    regexp = "'output' must be a single non-NA character string"
  )
  expect_error(
    bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f, output = NA_character_),
    regexp = "'output' must be a single non-NA character string"
  )
  expect_error(
    bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f, output = c("a", "b")),
    regexp = "'output' must be a single non-NA character string"
  )
})

test_that("bioclim_engine() stops for wrong number of climate files", {
  skip_without_gdal()
  skip_without_terra()

  tmpdir <- tempfile("be_val_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  # Create only 3 files (invalid — must be 1 or 12)
  bad_files <- make_monthly_files(rep(10, 12L), tmpdir, "bad")[1:3]
  good_f    <- make_monthly_files(std_tas_vals, tmpdir, "tas")

  expect_error(
    bioclim_engine(bad_files, good_f, good_f, good_f),
    regexp = "character vector of length 1 or 12"
  )
})

test_that("bioclim_engine() stops for missing climate files", {
  skip_without_gdal()
  skip_without_terra()

  tmpdir <- tempfile("be_val_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  tas_f    <- make_monthly_files(std_tas_vals,    tmpdir, "tas")
  tasmax_f <- make_monthly_files(std_tasmax_vals, tmpdir, "tasmax")
  tasmin_f <- make_monthly_files(std_tasmin_vals, tmpdir, "tasmin")
  pr_f     <- make_monthly_files(std_pr_vals,     tmpdir, "pr")

  # Replace one real path with a non-existent path
  bad_pr_f <- pr_f
  bad_pr_f[6L] <- "/no/such/file_06.tif"

  expect_error(
    bioclim_engine(tas_f, tasmax_f, tasmin_f, bad_pr_f),
    regexp = "do not exist"
  )
})

test_that("bioclim_engine() stops for invalid 'threads' argument", {
  skip_without_gdal()
  skip_without_terra()

  tmpdir <- tempfile("be_val_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  tas_f    <- make_monthly_files(std_tas_vals,    tmpdir, "tas")
  tasmax_f <- make_monthly_files(std_tasmax_vals, tmpdir, "tasmax")
  tasmin_f <- make_monthly_files(std_tasmin_vals, tmpdir, "tasmin")
  pr_f     <- make_monthly_files(std_pr_vals,     tmpdir, "pr")

  expect_error(
    bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f, threads = 0L),
    regexp = "'threads' must be a finite scalar integer >= 1"
  )
  expect_error(
    bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f, threads = -1L),
    regexp = "'threads' must be a finite scalar integer >= 1"
  )
  expect_error(
    bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f, threads = NA_integer_),
    regexp = "'threads' must be a finite scalar integer >= 1"
  )
})

test_that("bioclim_engine() stops for invalid 'tile_size' argument", {
  skip_without_gdal()
  skip_without_terra()

  tmpdir <- tempfile("be_val_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  tas_f    <- make_monthly_files(std_tas_vals,    tmpdir, "tas")
  tasmax_f <- make_monthly_files(std_tasmax_vals, tmpdir, "tasmax")
  tasmin_f <- make_monthly_files(std_tasmin_vals, tmpdir, "tasmin")
  pr_f     <- make_monthly_files(std_pr_vals,     tmpdir, "pr")

  expect_error(
    bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f, tile_size = 0L),
    regexp = "'tile_size' must be a finite scalar integer >= 1"
  )
})

test_that("bioclim_engine() stops for invalid 'variables' argument", {
  skip_without_gdal()
  skip_without_terra()

  tmpdir <- tempfile("be_val_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  tas_f    <- make_monthly_files(std_tas_vals,    tmpdir, "tas")
  tasmax_f <- make_monthly_files(std_tasmax_vals, tmpdir, "tasmax")
  tasmin_f <- make_monthly_files(std_tasmin_vals, tmpdir, "tasmin")
  pr_f     <- make_monthly_files(std_pr_vals,     tmpdir, "pr")

  expect_error(
    bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f, variables = integer(0)),
    regexp = "'variables' must be a non-empty integer vector"
  )
  expect_error(
    bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f, variables = 0L),
    regexp = "'variables' must contain values between 1 and 19"
  )
  expect_error(
    bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f, variables = 20L),
    regexp = "'variables' must contain values between 1 and 19"
  )
  expect_error(
    bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f, variables = c(1L, 1L)),
    regexp = "'variables' must not contain duplicates"
  )
})

# ── 2. Overwrite behavior ─────────────────────────────────────────────────────

test_that("bioclim_engine() stops when output files exist and overwrite = FALSE", {
  skip_without_gdal()
  skip_without_terra()

  tmpdir <- tempfile("be_ow_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  tas_f    <- make_monthly_files(std_tas_vals,    tmpdir, "tas")
  tasmax_f <- make_monthly_files(std_tasmax_vals, tmpdir, "tasmax")
  tasmin_f <- make_monthly_files(std_tasmin_vals, tmpdir, "tasmin")
  pr_f     <- make_monthly_files(std_pr_vals,     tmpdir, "pr")

  out <- file.path(tmpdir, "out_dir")
  dir.create(out)
  writeLines("dummy", file.path(out, "bio.tif"))

  expect_error(
    bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f, output = out,
                   overwrite = FALSE),
    regexp = "already exist"
  )
})

test_that("bioclim_engine() overwrites existing files when overwrite = TRUE", {
  skip_without_gdal()
  skip_without_terra()

  tmpdir <- tempfile("be_ow2_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  tas_f    <- make_monthly_files(std_tas_vals,    tmpdir, "tas")
  tasmax_f <- make_monthly_files(std_tasmax_vals, tmpdir, "tasmax")
  tasmin_f <- make_monthly_files(std_tasmin_vals, tmpdir, "tasmin")
  pr_f     <- make_monthly_files(std_pr_vals,     tmpdir, "pr")

  out <- file.path(tmpdir, "out_dir")
  dir.create(out)
  writeLines("dummy", file.path(out, "bio.tif"))

  expect_no_error(
    bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f,
                   output = out, overwrite = TRUE, tile_size = 2L)
  )
  expect_true(file.exists(file.path(out, "bio.tif")))
})

# ── 3. Full round-trip (GDAL + terra required) ────────────────────────────────

test_that("bioclim_engine() returns SpatRaster with 19 layers (12 monthly files)", {
  skip_without_gdal()
  skip_without_terra()

  tmpdir <- tempfile("be_rt_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  tas_f    <- make_monthly_files(std_tas_vals,    tmpdir, "tas")
  tasmax_f <- make_monthly_files(std_tasmax_vals, tmpdir, "tasmax")
  tasmin_f <- make_monthly_files(std_tasmin_vals, tmpdir, "tasmin")
  pr_f     <- make_monthly_files(std_pr_vals,     tmpdir, "pr")
  out      <- file.path(tmpdir, "bioclim_out")

  result <- bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f,
                            output = out, tile_size = 2L)

  expect_true(inherits(result, "SpatRaster"))
  expect_equal(terra::nlyr(result), 19L)
  expect_equal(terra::nrow(result), 3L)
  expect_equal(terra::ncol(result), 3L)

  # Output is a single multi-band GeoTIFF
  tif_files <- list.files(out, pattern = "\\.tif$")
  expect_equal(tif_files, "bio.tif")
})

test_that("bioclim_engine() BIO01 is finite and BIO12 approx sum(pr)", {
  skip_without_gdal()
  skip_without_terra()

  tmpdir <- tempfile("be_bio_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  tas_f    <- make_monthly_files(std_tas_vals,    tmpdir, "tas")
  tasmax_f <- make_monthly_files(std_tasmax_vals, tmpdir, "tasmax")
  tasmin_f <- make_monthly_files(std_tasmin_vals, tmpdir, "tasmin")
  pr_f     <- make_monthly_files(std_pr_vals,     tmpdir, "pr")

  result <- bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f, tile_size = 2L)

  bio01_vals <- as.numeric(terra::values(result[[1L]]))
  expect_true(all(is.finite(bio01_vals)),
              info = "BIO01 should be finite for valid inputs")

  expected_pr_sum <- sum(std_pr_vals)
  # BIO12 is the 12th variable in default 1:19
  bio12_layer <- result[[which(names(result) == "bio12")]]
  bio12_vals <- as.numeric(terra::values(bio12_layer))
  expect_true(all(abs(bio12_vals - expected_pr_sum) < 1e-6),
              info = "BIO12 should equal the sum of monthly precipitation")
})

test_that("bioclim_engine() works with 1 multi-band file per variable", {
  skip_without_gdal()
  skip_without_terra()

  tmpdir <- tempfile("be_mb_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  tas_f    <- make_multiband_file(std_tas_vals,    file.path(tmpdir, "tas.tif"))
  tasmax_f <- make_multiband_file(std_tasmax_vals, file.path(tmpdir, "tasmax.tif"))
  tasmin_f <- make_multiband_file(std_tasmin_vals, file.path(tmpdir, "tasmin.tif"))
  pr_f     <- make_multiband_file(std_pr_vals,     file.path(tmpdir, "pr.tif"))

  result <- bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f, tile_size = 2L)

  expect_true(inherits(result, "SpatRaster"))
  expect_equal(terra::nlyr(result), 19L)
})

test_that("bioclim_engine() accepts terra::SpatRaster inputs", {
  skip_without_gdal()
  skip_without_terra()

  tmpdir <- tempfile("be_sr_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  # Write files to disk so terra::sources() returns real paths
  tas_path    <- make_multiband_file(std_tas_vals,    file.path(tmpdir, "tas.tif"))
  tasmax_path <- make_multiband_file(std_tasmax_vals, file.path(tmpdir, "tasmax.tif"))
  tasmin_path <- make_multiband_file(std_tasmin_vals, file.path(tmpdir, "tasmin.tif"))
  pr_path     <- make_multiband_file(std_pr_vals,     file.path(tmpdir, "pr.tif"))

  tas_sr    <- terra::rast(tas_path)
  tasmax_sr <- terra::rast(tasmax_path)
  tasmin_sr <- terra::rast(tasmin_path)
  pr_sr     <- terra::rast(pr_path)

  result <- bioclim_engine(tas_sr, tasmax_sr, tasmin_sr, pr_sr,
                            tile_size = 2L)

  expect_true(inherits(result, "SpatRaster"))
  expect_equal(terra::nlyr(result), 19L)
})

# ── 4. Variable selection ─────────────────────────────────────────────────────

test_that("bioclim_engine() writes only requested variables", {
  skip_without_gdal()
  skip_without_terra()

  tmpdir <- tempfile("be_vars_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  tas_f    <- make_monthly_files(std_tas_vals,    tmpdir, "tas")
  tasmax_f <- make_monthly_files(std_tasmax_vals, tmpdir, "tasmax")
  tasmin_f <- make_monthly_files(std_tasmin_vals, tmpdir, "tasmin")
  pr_f     <- make_monthly_files(std_pr_vals,     tmpdir, "pr")
  out      <- file.path(tmpdir, "subset_out")

  result <- bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f,
                            output = out, variables = c(1L, 12L, 15L),
                            tile_size = 2L)

  expect_true(inherits(result, "SpatRaster"))
  expect_equal(terra::nlyr(result), 3L)
  expect_equal(names(result), c("bio01", "bio12", "bio15"))

  # Output is a single multi-band GeoTIFF
  tif_files <- list.files(out, pattern = "\\.tif$")
  expect_equal(tif_files, "bio.tif")
})

test_that("bioclim_engine() single variable produces 1-layer result", {
  skip_without_gdal()
  skip_without_terra()

  tmpdir <- tempfile("be_single_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  tas_f    <- make_monthly_files(std_tas_vals,    tmpdir, "tas")
  tasmax_f <- make_monthly_files(std_tasmax_vals, tmpdir, "tasmax")
  tasmin_f <- make_monthly_files(std_tasmin_vals, tmpdir, "tasmin")
  pr_f     <- make_monthly_files(std_pr_vals,     tmpdir, "pr")
  out      <- file.path(tmpdir, "single_out")

  result <- bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f,
                            output = out, variables = 12L,
                            tile_size = 2L)

  expect_true(inherits(result, "SpatRaster"))
  expect_equal(terra::nlyr(result), 1L)
  expect_equal(names(result), "bio12")

  # Output is a single multi-band GeoTIFF
  tif_files <- list.files(out, pattern = "\\.tif$")
  expect_equal(tif_files, "bio.tif")

  # Value check: BIO12 = sum of precipitation
  bio12_vals <- as.numeric(terra::values(result))
  expect_true(all(abs(bio12_vals - sum(std_pr_vals)) < 1e-6))
})

# ── 5. Mask parameter ─────────────────────────────────────────────────────────

test_that("bioclim_engine() accepts a raster mask file path", {
  skip_without_gdal()
  skip_without_terra()

  tmpdir <- tempfile("be_msk_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  tas_f    <- make_monthly_files(std_tas_vals,    tmpdir, "tas")
  tasmax_f <- make_monthly_files(std_tasmax_vals, tmpdir, "tasmax")
  tasmin_f <- make_monthly_files(std_tasmin_vals, tmpdir, "tasmin")
  pr_f     <- make_monthly_files(std_pr_vals,     tmpdir, "pr")

  # Create a binary mask raster (all 1s — fully covering)
  mask_path <- file.path(tmpdir, "mask.tif")
  mask_r <- terra::rast(
    nrows = 3L, ncols = 3L, nlyrs = 1L,
    xmin = 0, xmax = 1, ymin = 0, ymax = 1,
    crs = "EPSG:4326"
  )
  terra::values(mask_r) <- 1L
  terra::writeRaster(mask_r, mask_path, overwrite = TRUE)

  result <- bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f,
                            mask = mask_path, tile_size = 2L)

  expect_true(inherits(result, "SpatRaster"))
  expect_equal(terra::nlyr(result), 19L)
})

test_that("bioclim_engine() accepts a SpatRaster mask", {
  skip_without_gdal()
  skip_without_terra()

  tmpdir <- tempfile("be_msksr_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  tas_f    <- make_monthly_files(std_tas_vals,    tmpdir, "tas")
  tasmax_f <- make_monthly_files(std_tasmax_vals, tmpdir, "tasmax")
  tasmin_f <- make_monthly_files(std_tasmin_vals, tmpdir, "tasmin")
  pr_f     <- make_monthly_files(std_pr_vals,     tmpdir, "pr")

  mask_r <- terra::rast(
    nrows = 3L, ncols = 3L, nlyrs = 1L,
    xmin = 0, xmax = 1, ymin = 0, ymax = 1,
    crs = "EPSG:4326"
  )
  terra::values(mask_r) <- 1L

  result <- bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f,
                            mask = mask_r, tile_size = 2L)

  expect_true(inherits(result, "SpatRaster"))
  expect_equal(terra::nlyr(result), 19L)
})

test_that("bioclim_engine() stops when mask file path does not exist", {
  skip_without_gdal()
  skip_without_terra()

  tmpdir <- tempfile("be_mskbad_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  tas_f    <- make_monthly_files(std_tas_vals,    tmpdir, "tas")
  tasmax_f <- make_monthly_files(std_tasmax_vals, tmpdir, "tasmax")
  tasmin_f <- make_monthly_files(std_tasmin_vals, tmpdir, "tasmin")
  pr_f     <- make_monthly_files(std_pr_vals,     tmpdir, "pr")

  expect_error(
    bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f,
                   mask = "/no/such/mask.tif"),
    regexp = "'mask' file does not exist"
  )
})

test_that("bioclim_engine() stops for invalid mask type", {
  skip_without_gdal()
  skip_without_terra()

  tmpdir <- tempfile("be_msktype_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  tas_f    <- make_monthly_files(std_tas_vals,    tmpdir, "tas")
  tasmax_f <- make_monthly_files(std_tasmax_vals, tmpdir, "tasmax")
  tasmin_f <- make_monthly_files(std_tasmin_vals, tmpdir, "tasmin")
  pr_f     <- make_monthly_files(std_pr_vals,     tmpdir, "pr")

  expect_error(
    bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f, mask = 42L),
    regexp = "'mask' must be"
  )
})
