# Tests for block-based raster processing (bioclim_raster)
# Requires the 'terra' package; tests are skipped if terra is not installed.

skip_if_no_terra <- function() {
  skip_if_not_installed("terra")
}

# Reference bioclim result for the standard linear monthly test pattern
# (tas = 1:12, tasmax = 2:13, tasmin = 0:11, pr = 1:12)
ref_bioclim <- bioclim(1:12, 2:13, 0:11, 1:12)

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
  m_tas    <- matrix(rep(1:12, 5), nrow = 5, byrow = TRUE)
  m_tasmax <- matrix(rep(2:13, 5), nrow = 5, byrow = TRUE)
  m_tasmin <- matrix(rep(0:11, 5), nrow = 5, byrow = TRUE)
  m_pr     <- matrix(rep(1:12, 5), nrow = 5, byrow = TRUE)

  result <- bioclim_block(m_tas, m_tasmax, m_tasmin, m_pr)
  for (i in seq_len(5)) {
    expect_equal(result[i, ], unname(ref_bioclim), tolerance = 1e-6)
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
  expect_error(rxbioclim:::validate_spatraster(matrix(1:12, 1, 12), "x"), "must be a SpatRaster")
})

test_that("validate_spatraster errors on wrong number of layers", {
  skip_if_no_terra()
  r_wrong <- terra::rast(nrows = 2, ncols = 2, nlyr = 6L)
  expect_error(rxbioclim:::validate_spatraster(r_wrong, "tas"), "must have 12 layers")
})

test_that("validate_spatraster passes for valid 12-layer SpatRaster", {
  skip_if_no_terra()
  r_ok <- terra::rast(nrows = 2, ncols = 2, nlyr = 12L)
  expect_invisible(rxbioclim:::validate_spatraster(r_ok, "tas"))
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

  result_vals <- terra::values(result)

  # Every cell should match the reference (all cells have identical monthly data)
  for (i in seq_len(nrow(result_vals))) {
    expect_equal(result_vals[i, ], unname(ref_bioclim), tolerance = 1e-4)
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
  result_vals <- terra::values(result)
  expect_equal(result_vals[1, ], unname(ref_bioclim), tolerance = 1e-4)
})

test_that("bioclim_raster writes output to disk when filename is supplied", {
  skip_if_no_terra()
  rasts <- make_test_rasters()
  outfile <- tempfile(fileext = ".tif")
  on.exit(unlink(outfile), add = TRUE)

  result <- bioclim_raster(
    rasts$tas, rasts$tasmax, rasts$tasmin, rasts$pr,
    filename = outfile,
    overwrite = TRUE
  )

  expect_s4_class(result, "SpatRaster")
  expect_true(file.exists(outfile))
  expect_equal(terra::nlyr(result), 19L)

  result_vals <- terra::values(result)
  expect_equal(result_vals[1, ], unname(ref_bioclim), tolerance = 1e-4)
})

test_that("bioclim_raster works with ncores > 1", {
  skip_if_no_terra()
  if (is.na(parallel::detectCores()) || parallel::detectCores() < 2L) {
    skip("Parallel test requires at least 2 cores")
  }

  rasts <- make_test_rasters()
  result <- bioclim_raster(
    rasts$tas, rasts$tasmax, rasts$tasmin, rasts$pr,
    ncores = 2L
  )

  expect_equal(terra::nlyr(result), 19L)
  result_vals <- terra::values(result)
  expect_equal(result_vals[1, ], unname(ref_bioclim), tolerance = 1e-4)
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

test_that("bioclim_raster rejects geometrically mismatched rasters", {
  skip_if_no_terra()
  rasts  <- make_test_rasters()
  r_diff <- make_test_rasters(nrows = 2L, ncols = 3L)$tas
  expect_error(
    bioclim_raster(rasts$tas, r_diff, rasts$tasmin, rasts$pr),
    "same extent"
  )
})

test_that("bioclim_raster rejects invalid ncores", {
  skip_if_no_terra()
  rasts <- make_test_rasters()
  expect_error(
    bioclim_raster(rasts$tas, rasts$tasmax, rasts$tasmin, rasts$pr, ncores = NA),
    "'ncores' must be a finite scalar integer >= 1"
  )
  expect_error(
    bioclim_raster(rasts$tas, rasts$tasmax, rasts$tasmin, rasts$pr, ncores = 0L),
    "'ncores' must be a finite scalar integer >= 1"
  )
})

test_that("bioclim_raster rejects invalid n_blocks", {
  skip_if_no_terra()
  rasts <- make_test_rasters()
  expect_error(
    bioclim_raster(rasts$tas, rasts$tasmax, rasts$tasmin, rasts$pr, n_blocks = 0L),
    "'n_blocks' must be a finite scalar integer >= 1"
  )
  expect_error(
    bioclim_raster(rasts$tas, rasts$tasmax, rasts$tasmin, rasts$pr, n_blocks = NA),
    "'n_blocks' must be a finite scalar integer >= 1"
  )
})

test_that("bioclim_raster writes output to disk when filename is supplied", {
  skip_if_no_terra()
  rasts   <- make_test_rasters()
  outfile <- tempfile(fileext = ".tif")
  on.exit(unlink(outfile), add = TRUE)

  result <- bioclim_raster(
    rasts$tas, rasts$tasmax, rasts$tasmin, rasts$pr,
    filename  = outfile,
    overwrite = TRUE
  )

  expect_true(inherits(result, "SpatRaster"))
  expect_true(file.exists(outfile))
  expect_equal(terra::nlyr(result), 19L)

  result_vals <- terra::values(result)
  expect_equal(result_vals[1, ], unname(ref_bioclim), tolerance = 1e-4)
})

test_that("bioclim_raster works with ncores > 1", {
  skip_if_no_terra()
  skip_if_not_installed("parallel")
  n_avail <- parallel::detectCores()
  if (is.na(n_avail) || n_avail < 2L) skip("Parallel test requires at least 2 cores")

  rasts  <- make_test_rasters()
  result <- bioclim_raster(
    rasts$tas, rasts$tasmax, rasts$tasmin, rasts$pr,
    ncores = 2L
  )

  expect_equal(terra::nlyr(result), 19L)
  result_vals <- terra::values(result)
  expect_equal(result_vals[1, ], unname(ref_bioclim), tolerance = 1e-4)
})
