# Tests for block-based raster processing (bioclim_raster)
# Requires the 'terra' package; tests are skipped if terra is not installed.

skip_if_no_terra <- function() {
  skip_if_not_installed("terra")
}

# ── Helper: build minimal SpatRasters ───────────────────────────────────────

make_test_rasters <- function(nrows = 4L, ncols = 3L) {
  # Each cell gets the same linear monthly pattern: tas[m] = m, etc.
  n_cells <- nrows * ncols
  # Replicate the mock monthly vectors across all cells
  make_rast <- function(monthly_vals) {
    r <- terra::rast(
      nrows = nrows, ncols = ncols, nlyr = 12L,
      xmin = 0, xmax = ncols, ymin = 0, ymax = nrows,
      crs = "EPSG:4326"
    )
    # values() layout: row-major, one row per cell, one column per layer
    terra::values(r) <- matrix(
      rep(monthly_vals, times = n_cells),
      nrow = n_cells, ncol = 12L, byrow = TRUE
    )
    r
  }
  list(
    tas    = make_rast(1:12),
    tasmax = make_rast(2:13),
    tasmin = make_rast(0:11),
    pr     = make_rast(1:12)
  )
}

# ── bioclim_block() ──────────────────────────────────────────────────────────

test_that("bioclim_block returns a 19-column matrix", {
  m_tas    <- matrix(rep(1:12, 3), nrow = 3, byrow = TRUE)
  m_tasmax <- matrix(rep(2:13, 3), nrow = 3, byrow = TRUE)
  m_tasmin <- matrix(rep(0:11, 3), nrow = 3, byrow = TRUE)
  m_pr     <- matrix(rep(1:12, 3), nrow = 3, byrow = TRUE)

  result <- bioclim_block(m_tas, m_tasmax, m_tasmin, m_pr)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3L)
  expect_equal(ncol(result), 19L)
})

test_that("bioclim_block matches bioclim() for each row", {
  ref <- bioclim(1:12, 2:13, 0:11, 1:12)

  m_tas    <- matrix(rep(1:12, 5), nrow = 5, byrow = TRUE)
  m_tasmax <- matrix(rep(2:13, 5), nrow = 5, byrow = TRUE)
  m_tasmin <- matrix(rep(0:11, 5), nrow = 5, byrow = TRUE)
  m_pr     <- matrix(rep(1:12, 5), nrow = 5, byrow = TRUE)

  result <- bioclim_block(m_tas, m_tasmax, m_tasmin, m_pr)
  for (i in seq_len(5)) {
    expect_equal(result[i, ], unname(ref), tolerance = 1e-6)
  }
})

test_that("bioclim_block returns NA row when any input has NA", {
  m_tas    <- matrix(rep(1:12, 3), nrow = 3, byrow = TRUE)
  m_tasmax <- matrix(rep(2:13, 3), nrow = 3, byrow = TRUE)
  m_tasmin <- matrix(rep(0:11, 3), nrow = 3, byrow = TRUE)
  m_pr     <- matrix(rep(1:12, 3), nrow = 3, byrow = TRUE)

  # Insert NA into the second cell's tas values
  m_tas[2, 1] <- NA

  result <- bioclim_block(m_tas, m_tasmax, m_tasmin, m_pr)
  expect_false(anyNA(result[1, ]))
  expect_true(all(is.na(result[2, ])))
  expect_false(anyNA(result[3, ]))
})

# ── validate_spatraster() ────────────────────────────────────────────────────

test_that("validate_spatraster errors on non-SpatRaster input", {
  skip_if_no_terra()
  expect_error(validate_spatraster(matrix(1:12, 1, 12), "x"), "must be a SpatRaster")
})

test_that("validate_spatraster errors on wrong number of layers", {
  skip_if_no_terra()
  r_wrong <- terra::rast(nrows = 2, ncols = 2, nlyr = 6L)
  expect_error(validate_spatraster(r_wrong, "tas"), "must have 12 layers")
})

test_that("validate_spatraster passes for valid 12-layer SpatRaster", {
  skip_if_no_terra()
  r_ok <- terra::rast(nrows = 2, ncols = 2, nlyr = 12L)
  expect_invisible(validate_spatraster(r_ok, "tas"))
})

# ── bioclim_raster() ─────────────────────────────────────────────────────────

test_that("bioclim_raster is an exported function in the package namespace", {
  expect_true(is.function(bioclim_raster))
  expect_true("bioclim_raster" %in% getNamespaceExports("rxbioclim"))
})

test_that("bioclim_raster returns SpatRaster with 19 layers", {
  skip_if_no_terra()
  rasts  <- make_test_rasters()
  result <- bioclim_raster(rasts$tas, rasts$tasmax, rasts$tasmin, rasts$pr)
  expect_true(inherits(result, "SpatRaster"))
  expect_equal(terra::nlyr(result), 19L)
})

test_that("bioclim_raster output layer names are bio01 through bio19", {
  skip_if_no_terra()
  rasts  <- make_test_rasters()
  result <- bioclim_raster(rasts$tas, rasts$tasmax, rasts$tasmin, rasts$pr)
  expect_equal(names(result), paste0("bio", sprintf("%02d", 1:19)))
})

test_that("bioclim_raster result matches bioclim() per-pixel", {
  skip_if_no_terra()
  rasts  <- make_test_rasters()
  result <- bioclim_raster(rasts$tas, rasts$tasmax, rasts$tasmin, rasts$pr)

  ref <- bioclim(1:12, 2:13, 0:11, 1:12)
  result_vals <- terra::values(result)

  # Every cell should match the reference (all cells have identical monthly data)
  for (i in seq_len(nrow(result_vals))) {
    expect_equal(result_vals[i, ], unname(ref), tolerance = 1e-4)
  }
})

test_that("bioclim_raster preserves spatial extent and CRS", {
  skip_if_no_terra()
  rasts  <- make_test_rasters()
  result <- bioclim_raster(rasts$tas, rasts$tasmax, rasts$tasmin, rasts$pr)

  expect_equal(terra::ext(result), terra::ext(rasts$tas))
  expect_equal(terra::crs(result), terra::crs(rasts$tas))
  expect_equal(terra::nrow(result), terra::nrow(rasts$tas))
  expect_equal(terra::ncol(result), terra::ncol(rasts$tas))
})

test_that("bioclim_raster handles NA cells correctly", {
  skip_if_no_terra()
  rasts <- make_test_rasters()
  # Introduce NA into the first layer of tas for cell 1
  v <- terra::values(rasts$tas)
  v[1, 1] <- NA
  terra::values(rasts$tas) <- v

  result <- bioclim_raster(rasts$tas, rasts$tasmax, rasts$tasmin, rasts$pr)
  result_vals <- terra::values(result)

  # Cell 1 should be all NA
  expect_true(all(is.na(result_vals[1, ])))
  # Other cells should be non-NA
  expect_false(anyNA(result_vals[2, ]))
})

test_that("bioclim_raster works with explicit n_blocks parameter", {
  skip_if_no_terra()
  rasts  <- make_test_rasters()
  result <- bioclim_raster(rasts$tas, rasts$tasmax, rasts$tasmin, rasts$pr,
                            n_blocks = 2L)
  expect_equal(terra::nlyr(result), 19L)
  ref <- bioclim(1:12, 2:13, 0:11, 1:12)
  result_vals <- terra::values(result)
  expect_equal(result_vals[1, ], unname(ref), tolerance = 1e-4)
})

test_that("bioclim_raster rejects non-SpatRaster input", {
  skip_if_no_terra()
  rasts <- make_test_rasters()
  expect_error(
    bioclim_raster(matrix(1, 1, 12), rasts$tasmax, rasts$tasmin, rasts$pr),
    "must be a SpatRaster"
  )
})

test_that("bioclim_raster rejects SpatRaster with wrong layer count", {
  skip_if_no_terra()
  rasts   <- make_test_rasters()
  r_wrong <- terra::rast(nrows = 4, ncols = 3, nlyr = 6L)
  expect_error(
    bioclim_raster(r_wrong, rasts$tasmax, rasts$tasmin, rasts$pr),
    "must have 12 layers"
  )
})
