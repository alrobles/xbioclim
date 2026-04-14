# Tests for polygon masking support (Issue #25).
#
# All GDAL-dependent tests skip when the package was built without GDAL
# support, keeping the test suite green on platforms without libgdal-dev.
#
# The tiny.tif shipped in inst/extdata (3×3 pixels, no CRS) is used as the
# reference raster.  A minimal GeoJSON polygon covering the entire raster
# extent is generated on the fly — no external test files are needed.

# ── Helpers ──────────────────────────────────────────────────────────────────

.has_gdal_mask <- function() {
  tryCatch({
    gdal_can_open("/dev/null")
    TRUE
  }, error = function(e) {
    msg <- conditionMessage(e)
    if (grepl("built without GDAL", msg, fixed = TRUE)) FALSE else TRUE
  })
}

skip_without_gdal <- function() {
  if (!.has_gdal_mask()) {
    testthat::skip("Package built without GDAL support")
  }
}

tiny_tif <- function() {
  system.file("extdata", "tiny.tif", package = "rxbioclim")
}

# A GeoJSON polygon that covers the full extent of tiny.tif.
# tiny.tif has the default identity geotransform, so pixel coordinates are
# the spatial coordinates.  We cover [0,3] x [0,-3] in GDAL convention.
write_covering_geojson <- function(path) {
  geojson <- paste0(
    '{"type":"FeatureCollection","features":[{"type":"Feature",',
    '"geometry":{"type":"Polygon",',
    '"coordinates":[[[0,0],[3,0],[3,-3],[0,-3],[0,0]]]},',
    '"properties":{}}]}'
  )
  writeLines(geojson, path)
}

write_empty_geojson <- function(path) {
  geojson <- paste0(
    '{"type":"FeatureCollection","features":[]}'
  )
  writeLines(geojson, path)
}

# ── rasterize_mask_cpp ───────────────────────────────────────────────────────

test_that("rasterize_mask_cpp stops with clear message when GDAL is absent", {
  if (.has_gdal_mask()) {
    testthat::skip("GDAL is present — absence path not applicable")
  }
  expect_error(
    rasterize_mask_cpp("poly.geojson", "ref.tif", tempfile(fileext = ".tif")),
    regexp = "built without GDAL"
  )
})

test_that("rasterize_mask_cpp creates a mask raster for a covering polygon", {
  skip_without_gdal()

  ref  <- tiny_tif()
  poly <- tempfile(fileext = ".geojson")
  mask <- tempfile(fileext = ".tif")
  on.exit(unlink(c(poly, mask)), add = TRUE)

  write_covering_geojson(poly)
  rasterize_mask_cpp(poly, ref, mask)

  expect_true(file.exists(mask))
  expect_true(gdal_can_open(mask))
})

test_that("rasterize_mask_cpp mask has same dimensions as reference raster", {
  skip_without_gdal()

  ref  <- tiny_tif()
  poly <- tempfile(fileext = ".geojson")
  mask <- tempfile(fileext = ".tif")
  on.exit(unlink(c(poly, mask)), add = TRUE)

  write_covering_geojson(poly)
  rasterize_mask_cpp(poly, ref, mask)

  ref_info  <- gdal_info(ref)
  mask_info <- gdal_info(mask)

  expect_equal(mask_info$nrows,  ref_info$nrows)
  expect_equal(mask_info$ncols,  ref_info$ncols)
  expect_equal(mask_info$nbands, 1L)
})

test_that("rasterize_mask_cpp all pixels are 1 for fully covering polygon", {
  skip_without_gdal()
  testthat::skip_if_not_installed("terra")

  ref  <- tiny_tif()
  poly <- tempfile(fileext = ".geojson")
  mask <- tempfile(fileext = ".tif")
  on.exit(unlink(c(poly, mask)), add = TRUE)

  write_covering_geojson(poly)
  rasterize_mask_cpp(poly, ref, mask)

  r_mask <- terra::rast(mask)
  vals   <- as.numeric(terra::values(r_mask))
  expect_true(all(vals == 1L))
})

test_that("rasterize_mask_cpp stops for a non-existent vector source", {
  skip_without_gdal()

  ref  <- tiny_tif()
  mask <- tempfile(fileext = ".tif")
  on.exit(unlink(mask), add = TRUE)

  expect_error(
    rasterize_mask_cpp(tempfile(fileext = ".geojson"), ref, mask)
  )
})

# ── apply_mask_cpp ────────────────────────────────────────────────────────────

test_that("apply_mask_cpp stops with clear message when GDAL is absent", {
  if (.has_gdal_mask()) {
    testthat::skip("GDAL is present — absence path not applicable")
  }
  expect_error(
    apply_mask_cpp("input.tif", "mask.tif", tempfile(fileext = ".tif")),
    regexp = "built without GDAL"
  )
})

test_that("apply_mask_cpp creates output raster with same dimensions", {
  skip_without_gdal()

  ref    <- tiny_tif()
  poly   <- tempfile(fileext = ".geojson")
  mask   <- tempfile(fileext = ".tif")
  output <- tempfile(fileext = ".tif")
  on.exit(unlink(c(poly, mask, output)), add = TRUE)

  write_covering_geojson(poly)
  rasterize_mask_cpp(poly, ref, mask)
  apply_mask_cpp(ref, mask, output)

  expect_true(file.exists(output))

  ref_info <- gdal_info(ref)
  out_info <- gdal_info(output)
  expect_equal(out_info$nrows,  ref_info$nrows)
  expect_equal(out_info$ncols,  ref_info$ncols)
  expect_equal(out_info$nbands, ref_info$nbands)
})

test_that("apply_mask_cpp preserves all pixel values for a fully covering mask", {
  skip_without_gdal()
  testthat::skip_if_not_installed("terra")

  ref    <- tiny_tif()
  poly   <- tempfile(fileext = ".geojson")
  mask   <- tempfile(fileext = ".tif")
  output <- tempfile(fileext = ".tif")
  on.exit(unlink(c(poly, mask, output)), add = TRUE)

  write_covering_geojson(poly)
  rasterize_mask_cpp(poly, ref, mask)
  apply_mask_cpp(ref, mask, output)

  ref_vals <- as.numeric(terra::values(terra::rast(ref)))
  out_vals <- as.numeric(terra::values(terra::rast(output)))

  # No NAs expected — all pixels covered.
  expect_false(anyNA(out_vals))
  expect_equal(out_vals, ref_vals, tolerance = 1e-6)
})

test_that("apply_mask_cpp sets pixels to NA where mask == 0", {
  skip_without_gdal()
  testthat::skip_if_not_installed("terra")

  ref    <- tiny_tif()
  poly   <- tempfile(fileext = ".geojson")
  mask   <- tempfile(fileext = ".tif")
  output <- tempfile(fileext = ".tif")
  on.exit(unlink(c(poly, mask, output)), add = TRUE)

  # A polygon that covers only the top-left pixel (roughly).
  geojson <- paste0(
    '{"type":"FeatureCollection","features":[{"type":"Feature",',
    '"geometry":{"type":"Polygon",',
    '"coordinates":[[[0,0],[1,0],[1,-1],[0,-1],[0,0]]]},',
    '"properties":{}}]}'
  )
  writeLines(geojson, poly)

  rasterize_mask_cpp(poly, ref, mask)
  apply_mask_cpp(ref, mask, output)

  r_out  <- terra::rast(output)
  vals   <- as.numeric(terra::values(r_out))

  # At least some pixels should be NA (those outside the tiny polygon).
  expect_true(anyNA(vals))
})

test_that("apply_mask_cpp stops when mask and input dimensions differ", {
  skip_without_gdal()

  ref  <- tiny_tif()

  # Create a 1×1 mask — does not match 3×3 tiny.tif.
  small_mask <- tempfile(fileext = ".tif")
  output     <- tempfile(fileext = ".tif")
  on.exit(unlink(c(small_mask, output)), add = TRUE)

  # Build a 1×1 GeoTIFF mask manually via GdalWriter through gdal_info proxy.
  # We use terra (if available) to create it.
  testthat::skip_if_not_installed("terra")
  terra::writeRaster(
    terra::rast(matrix(1L, 1, 1)),
    small_mask, overwrite = TRUE,
    datatype = "INT1U"
  )

  expect_error(apply_mask_cpp(ref, small_mask, output))
})

# ── create_mask (R wrapper) ───────────────────────────────────────────────────

test_that("create_mask returns mask path invisibly for a file-path polygon", {
  skip_without_gdal()

  ref  <- tiny_tif()
  poly <- tempfile(fileext = ".geojson")
  mask <- tempfile(fileext = ".tif")
  on.exit(unlink(c(poly, mask)), add = TRUE)

  write_covering_geojson(poly)

  result <- create_mask(poly, ref, output = mask)
  expect_equal(result, mask)
  expect_true(file.exists(mask))
})

test_that("create_mask with apply_to returns applied-output path", {
  skip_without_gdal()

  ref    <- tiny_tif()
  poly   <- tempfile(fileext = ".geojson")
  mask   <- tempfile(fileext = ".tif")
  output <- tempfile(fileext = ".tif")
  on.exit(unlink(c(poly, mask, output)), add = TRUE)

  write_covering_geojson(poly)

  result <- create_mask(poly, ref, output = mask,
                        apply_to = ref, apply_output = output)
  expect_equal(result, output)
  expect_true(file.exists(output))
})

test_that("create_mask accepts an sf polygon object", {
  skip_without_gdal()
  testthat::skip_if_not_installed("sf")

  ref  <- tiny_tif()
  mask <- tempfile(fileext = ".tif")
  on.exit(unlink(mask), add = TRUE)

  poly_sf <- sf::st_sfc(
    sf::st_polygon(list(matrix(
      c(0, 0, 3, 0, 3, -3, 0, -3, 0, 0),
      ncol = 2, byrow = TRUE
    )))
  )
  poly_sf <- sf::st_sf(geometry = poly_sf)

  result <- create_mask(poly_sf, ref, output = mask)
  expect_equal(result, mask)
  expect_true(file.exists(mask))
})

test_that("create_mask accepts a SpatRaster as a pre-made mask", {
  testthat::skip_if_not_installed("terra")

  mask_rast <- terra::rast(matrix(c(1L, 0L, 1L, 0L), nrow = 2))
  output    <- tempfile(fileext = ".tif")
  on.exit(unlink(output), add = TRUE)

  result <- create_mask(mask_rast, output = output)
  expect_equal(result, output)
  expect_true(file.exists(output))
})
