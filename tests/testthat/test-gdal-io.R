# Tests for GDAL I/O integration (Issue #22).
#
# All tests skip when the package was built without GDAL support so the test
# suite remains green on platforms where GDAL is absent (Windows CI, basic
# Linux without libgdal-dev, etc.).
#
# A tiny 3×3 Float32 GeoTIFF is shipped in inst/extdata/tiny.tif for smoke
# testing.  It contains pixel values 1..9 (row-major) with no geotransform
# or CRS set.

# ── Helpers ──────────────────────────────────────────────────────────────────

# Return TRUE if rxbioclim was compiled with GDAL support.
.has_gdal <- function() {
  tryCatch({
    # gdal_can_open() stops() when GDAL is absent.
    gdal_can_open(tiny_tif())
    TRUE
  }, error = function(e) {
    msg <- conditionMessage(e)
    if (grepl("built without GDAL", msg, fixed = TRUE)) {
      FALSE
    } else {
      # GDAL is present but the file is invalid — that's fine.
      TRUE
    }
  })
}

skip_without_gdal <- function() {
  if (!.has_gdal()) {
    testthat::skip("Package built without GDAL support")
  }
}

tiny_tif <- function() {
  system.file("extdata", "tiny.tif", package = "rxbioclim")
}

# ── gdal_can_open ────────────────────────────────────────────────────────────

test_that("gdal_can_open stops with clear message when GDAL is absent", {
  if (.has_gdal()) {
    testthat::skip("GDAL is present — testing absence path not applicable")
  }
  expect_error(gdal_can_open("anything.tif"),
               regexp = "built without GDAL")
})

test_that("gdal_can_open returns TRUE for valid TIFF when GDAL present", {
  skip_without_gdal()
  path <- tiny_tif()
  expect_true(file.exists(path))
  expect_true(gdal_can_open(path))
})

test_that("gdal_can_open returns FALSE for non-existent file when GDAL present", {
  skip_without_gdal()
  expect_false(gdal_can_open(tempfile(fileext = ".tif")))
})

# ── gdal_info ────────────────────────────────────────────────────────────────

test_that("gdal_info stops with clear message when GDAL is absent", {
  if (.has_gdal()) {
    testthat::skip("GDAL is present — testing absence path not applicable")
  }
  expect_error(gdal_info("anything.tif"),
               regexp = "built without GDAL")
})

test_that("gdal_info returns list with expected names for tiny.tif", {
  skip_without_gdal()
  path <- tiny_tif()
  info  <- gdal_info(path)

  # Must be a named list.
  expect_type(info, "list")
  expect_named(info,
    c("path", "nrows", "ncols", "nbands",
      "geotransform", "crs", "scale", "offset"),
    ignore.order = FALSE)
})

test_that("gdal_info returns correct dimensions for 3x3 single-band tiny.tif", {
  skip_without_gdal()
  info <- gdal_info(tiny_tif())

  expect_equal(info$nrows,  3L)
  expect_equal(info$ncols,  3L)
  expect_equal(info$nbands, 1L)
})

test_that("gdal_info geotransform is numeric vector of length 6", {
  skip_without_gdal()
  info <- gdal_info(tiny_tif())
  expect_type(info$geotransform, "double")
  expect_length(info$geotransform, 6L)
})

test_that("gdal_info scale and offset vectors have one element per band", {
  skip_without_gdal()
  info <- gdal_info(tiny_tif())
  expect_length(info$scale,  info$nbands)
  expect_length(info$offset, info$nbands)
})

test_that("gdal_info default scale is 1 and offset is 0 for unscaled raster", {
  skip_without_gdal()
  info <- gdal_info(tiny_tif())
  expect_equal(info$scale,  rep(1.0, info$nbands), tolerance = 1e-12)
  expect_equal(info$offset, rep(0.0, info$nbands), tolerance = 1e-12)
})

test_that("gdal_info path element matches input", {
  skip_without_gdal()
  path <- tiny_tif()
  info <- gdal_info(path)
  expect_equal(info$path, path)
})
