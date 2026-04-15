# tests/testthat/test-cuda.R
# Tests for CUDA device detection and device dispatch in bioclim_engine().
#
# All GPU-dependent tests skip gracefully when no CUDA device is present.
# These tests focus on the R-level API and fallback behaviour; they do NOT
# require actual GPU hardware to pass on CI.

# ── Helpers ──────────────────────────────────────────────────────────────────

skip_without_gdal <- function() {
  if (!has_gdal()) testthat::skip("Package built without GDAL support")
}

skip_without_terra <- function() {
  testthat::skip_if_not_installed("terra")
}

# Build a tiny (3×3, 12-band) GeoTIFF filled with a constant value.
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

# ── cuda_device_count ────────────────────────────────────────────────────────

test_that("cuda_device_count returns a non-negative integer", {
  result <- cuda_device_count()
  expect_type(result, "integer")
  expect_true(result >= 0L)
})

# ── has_cuda ────────────────────────────────────────────────────────────────

test_that("has_cuda returns a logical scalar", {
  result <- has_cuda()
  expect_type(result, "logical")
  expect_length(result, 1L)
  expect_false(is.na(result))
})

test_that("has_cuda agrees with cuda_device_count", {
  expect_equal(has_cuda(), cuda_device_count() > 0L)
})

# ── cuda_device_info / cuda_info ─────────────────────────────────────────────

test_that("cuda_device_info returns a list", {
  info <- cuda_device_info()
  expect_type(info, "list")
})

test_that("cuda_info returns a list", {
  info <- cuda_info()
  expect_type(info, "list")
})

test_that("cuda_info fields are correct when a GPU is present", {
  skip_if(!has_cuda(), "No CUDA device present")
  info <- cuda_info()
  expect_named(info, c("name", "memory_gb", "compute_capability"))
  expect_type(info$name, "character")
  expect_type(info$memory_gb, "double")
  expect_true(info$memory_gb > 0)
  expect_type(info$compute_capability, "character")
  expect_match(info$compute_capability, "^[0-9]+\\.[0-9]+$")
})

test_that("cuda_info returns empty list when no GPU", {
  skip_if(has_cuda(), "CUDA device present — skipping no-GPU test")
  expect_length(cuda_info(), 0L)
})

# ── engine_set_device ─────────────────────────────────────────────────────────

test_that("engine_set_device accepts valid device strings", {
  eng <- engine_create()
  expect_no_error(engine_set_device(eng, "auto"))
  expect_no_error(engine_set_device(eng, "cpu"))
  expect_no_error(engine_set_device(eng, "gpu"))
})

test_that("engine_set_device rejects unknown device strings", {
  eng <- engine_create()
  expect_error(engine_set_device(eng, "tpu"), "unknown device")
})

# ── bioclim_engine device parameter ──────────────────────────────────────────

test_that("bioclim_engine device='cpu' always works", {
  skip_without_gdal()
  skip_without_terra()

  tmp <- tempdir()
  tas    <- make_tiny_raster(15.0, file.path(tmp, "cuda_tas.tif"))
  tasmax <- make_tiny_raster(20.0, file.path(tmp, "cuda_tasmax.tif"))
  tasmin <- make_tiny_raster(10.0, file.path(tmp, "cuda_tasmin.tif"))
  pr     <- make_tiny_raster(50.0, file.path(tmp, "cuda_pr.tif"))
  out    <- tempfile(fileext = ".tif")

  result <- bioclim_engine(
    tas, tasmax, tasmin, pr,
    output   = out,
    device   = "cpu",
    overwrite = TRUE
  )
  expect_true(file.exists(out))
  if (inherits(result, "SpatRaster")) {
    expect_equal(terra::nlyr(result), 19L)
  } else {
    expect_type(result, "character")
  }
})

test_that("bioclim_engine device='gpu' falls back to CPU with a warning when no GPU", {
  skip_without_gdal()
  skip_without_terra()
  skip_if(has_cuda(), "CUDA device present — fallback test only applies when no GPU")

  tmp <- tempdir()
  tas    <- make_tiny_raster(15.0, file.path(tmp, "cuda_fb_tas.tif"))
  tasmax <- make_tiny_raster(20.0, file.path(tmp, "cuda_fb_tasmax.tif"))
  tasmin <- make_tiny_raster(10.0, file.path(tmp, "cuda_fb_tasmin.tif"))
  pr     <- make_tiny_raster(50.0, file.path(tmp, "cuda_fb_pr.tif"))
  out    <- tempfile(fileext = ".tif")

  expect_warning(
    bioclim_engine(
      tas, tasmax, tasmin, pr,
      output    = out,
      device    = "gpu",
      overwrite = TRUE
    ),
    regexp = "CUDA GPU requested but not available",
    fixed  = TRUE
  )
  expect_true(file.exists(out))
})

test_that("bioclim_engine device='auto' runs without error", {
  skip_without_gdal()
  skip_without_terra()

  tmp <- tempdir()
  tas    <- make_tiny_raster(15.0, file.path(tmp, "cuda_auto_tas.tif"))
  tasmax <- make_tiny_raster(20.0, file.path(tmp, "cuda_auto_tasmax.tif"))
  tasmin <- make_tiny_raster(10.0, file.path(tmp, "cuda_auto_tasmin.tif"))
  pr     <- make_tiny_raster(50.0, file.path(tmp, "cuda_auto_pr.tif"))
  out    <- tempfile(fileext = ".tif")

  expect_no_error(
    bioclim_engine(
      tas, tasmax, tasmin, pr,
      output    = out,
      device    = "auto",
      overwrite = TRUE
    )
  )
  expect_true(file.exists(out))
})

# ── BIO15 NaN/formula alignment: GPU == CPU ───────────────────────────────────

test_that("BIO15 is NaN for zero precipitation (CPU path)", {
  skip_without_gdal()
  skip_without_terra()

  tmp <- tempdir()
  tas    <- make_tiny_raster(15.0, file.path(tmp, "bio15_nan_tas.tif"))
  tasmax <- make_tiny_raster(20.0, file.path(tmp, "bio15_nan_tasmax.tif"))
  tasmin <- make_tiny_raster(10.0, file.path(tmp, "bio15_nan_tasmin.tif"))
  pr     <- make_tiny_raster(0.0,  file.path(tmp, "bio15_nan_pr.tif"))
  out    <- tempfile(fileext = ".tif")

  bioclim_engine(
    tas, tasmax, tasmin, pr,
    output   = out,
    device   = "cpu",
    overwrite = TRUE
  )
  result <- terra::rast(out)
  bio15_vals <- as.numeric(terra::values(result[[15]]))
  expect_true(all(is.nan(bio15_vals) | is.na(bio15_vals)))
})

test_that("BIO15 GPU and CPU produce identical values (nonzero precipitation)", {
  skip_without_gdal()
  skip_without_terra()
  skip_if(!has_cuda(), "No CUDA device present")

  tmp <- tempdir()
  tas    <- make_tiny_raster(15.0, file.path(tmp, "bio15_cmp_tas.tif"))
  tasmax <- make_tiny_raster(20.0, file.path(tmp, "bio15_cmp_tasmax.tif"))
  tasmin <- make_tiny_raster(10.0, file.path(tmp, "bio15_cmp_tasmin.tif"))
  pr     <- make_tiny_raster(50.0, file.path(tmp, "bio15_cmp_pr.tif"))
  out_cpu <- tempfile(fileext = ".tif")
  out_gpu <- tempfile(fileext = ".tif")

  bioclim_engine(tas, tasmax, tasmin, pr, output = out_cpu,
                 device = "cpu", overwrite = TRUE)
  bioclim_engine(tas, tasmax, tasmin, pr, output = out_gpu,
                 device = "gpu", overwrite = TRUE)

  bio15_cpu <- as.numeric(terra::values(terra::rast(out_cpu)[[15]]))
  bio15_gpu <- as.numeric(terra::values(terra::rast(out_gpu)[[15]]))
  expect_equal(bio15_cpu, bio15_gpu, tolerance = 1e-9)
})

test_that("BIO15 GPU returns NaN for zero precipitation", {
  skip_without_gdal()
  skip_without_terra()
  skip_if(!has_cuda(), "No CUDA device present")

  tmp <- tempdir()
  tas    <- make_tiny_raster(15.0, file.path(tmp, "bio15_gpu_nan_tas.tif"))
  tasmax <- make_tiny_raster(20.0, file.path(tmp, "bio15_gpu_nan_tasmax.tif"))
  tasmin <- make_tiny_raster(10.0, file.path(tmp, "bio15_gpu_nan_tasmin.tif"))
  pr     <- make_tiny_raster(0.0,  file.path(tmp, "bio15_gpu_nan_pr.tif"))
  out    <- tempfile(fileext = ".tif")

  bioclim_engine(
    tas, tasmax, tasmin, pr,
    output   = out,
    device   = "gpu",
    overwrite = TRUE
  )
  result <- terra::rast(out)
  bio15_vals <- as.numeric(terra::values(result[[15]]))
  expect_true(all(is.nan(bio15_vals) | is.na(bio15_vals)))
})
