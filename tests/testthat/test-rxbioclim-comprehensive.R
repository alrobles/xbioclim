# ============================================================================
# Comprehensive Integration Test Suite for rxbioclim
# ============================================================================
#
# This script exercises every layer of the package:
#   Layer 1: Pure R primitives (primitives.R)
#   Layer 2: Individual BIO functions (bioclim.R) — R path
#   Layer 3: S4 class BioclimData — C++ backed via Rcpp
#   Layer 4: BioclimModel — external pointer (XPtr) path
#   Layer 5: ClimateBlock Rcpp module
#   Layer 6: bioclim_xt() — zero-copy C++ vectorized path
#   Layer 7: bioclim_raster() — terra block-loop path
#   Layer 8: BioclimEngine — GDAL tiled XPtr engine
#   Layer 9: bioclim_engine() — R user-facing API
#   Layer 10: Masking support (GDAL polygon masking)
#   Layer 11: Cross-layer numerical consistency
#   Layer 12: Edge cases, NA handling, stress tests
#
# All GDAL/terra-dependent tests skip gracefully on platforms without them.
# Designed for Issue #28 (backward compatibility) validation.
# ============================================================================

library(testthat)

# ── Shared test data ─────────────────────────────────────────────────────────

# Standard linear pattern (xbioclim convention)
tas_lin    <- as.numeric(1:12)
tasmax_lin <- as.numeric(2:13)
tasmin_lin <- as.numeric(0:11)
pr_lin     <- as.numeric(1:12)

# Reversed precipitation (anti-correlated with temperature)
pr_rev <- as.numeric(12:1)

# Constant climate
tas_const    <- rep(20, 12)
tasmax_const <- rep(25, 12)
tasmin_const <- rep(15, 12)
pr_const     <- rep(50, 12)

# Realistic temperate climate (Issue #24 test pattern)
tas_real    <- c(5, 7, 10, 14, 18, 22, 25, 24, 20, 15, 10,  6)
tasmax_real <- c(8, 10, 14, 18, 23, 28, 32, 31, 26, 19, 13,  9)
tasmin_real <- c(1,  3,  6, 10, 13, 17, 20, 19, 15, 10,  6,  2)
pr_real     <- c(60, 55, 48, 35, 28, 22, 18, 20, 35, 55, 65, 68)

# Reference values for linear pattern (verified against xbioclim C++)
ref_linear <- c(
  bio01 = 6.5,      bio02 = 2.0,      bio03 = 15.3846,
  bio04 = 345.2053, bio05 = 13.0,     bio06 = 0.0,
  bio07 = 13.0,     bio08 = 11.0,     bio09 = 2.0,
  bio10 = 11.0,     bio11 = 2.0,      bio12 = 78.0,
  bio13 = 12.0,     bio14 = 1.0,      bio15 = 53.1085,
  bio16 = 33.0,     bio17 = 6.0,      bio18 = 33.0,
  bio19 = 6.0
)

tol       <- 1e-4
tol_loose <- 1e-3

to_mat <- function(x) matrix(as.double(x), nrow = 1L, ncol = 12L)

# ── Skip helpers ─────────────────────────────────────────────────────────────

skip_without_terra <- function() testthat::skip_if_not_installed("terra")
skip_without_gdal  <- function() if (!has_gdal()) testthat::skip("No GDAL")
skip_without_sf    <- function() testthat::skip_if_not_installed("sf")

# ============================================================================
# LAYER 1: Pure R primitives
# ============================================================================
test_that("L1: sd_pop matches population SD formula", {
  x <- 1:12
  expected <- sqrt(sum((x - mean(x))^2) / length(x))
  expect_equal(sd_pop(x), expected, tolerance = 1e-10)
})
test_that("L1: sd_pop is zero for constant vector", {
  expect_equal(sd_pop(rep(42, 12)), 0)
})
test_that("L1: rolling_quarter_sum circular wrapping", {
  qs <- rolling_quarter_sum(1:12)
  expect_length(qs, 12)
  expect_equal(qs[1],  6)    # 1+2+3
  expect_equal(qs[10], 33)   # 10+11+12
  expect_equal(qs[11], 24)   # 11+12+1
  expect_equal(qs[12], 15)   # 12+1+2
})
test_that("L1: rolling_quarter_mean = rolling_quarter_sum / 3", {
  x <- runif(12) * 100
  qs <- rolling_quarter_sum(x)
  qm <- rolling_quarter_mean(x)
  expect_equal(qm, qs / 3, tolerance = 1e-12)
})
test_that("L1: quarter_argmax and quarter_argmin correct for linear", {
  expect_equal(quarter_argmax(1:12), 10)
  expect_equal(quarter_argmin(1:12), 1)
})
test_that("L1: quarter_values extracts correct months", {
  expect_equal(quarter_values(1:12, 1), c(1, 2, 3))
  expect_equal(quarter_values(1:12, 11), c(11, 12, 1))
})
test_that("L1: validate_monthly rejects bad input", {
  expect_error(validate_monthly(letters[1:12]), "must be numeric")
  expect_error(validate_monthly(1:6), "must have length 12")
  expect_silent(validate_monthly(1:12))
})

# ============================================================================
# LAYER 2: Individual BIO functions (R path)
# ============================================================================
test_that("L2: All 19 BIO functions produce correct linear-pattern values", {
  expect_equal(bio01(tas_lin),                ref_linear[["bio01"]], tolerance = tol)
  expect_equal(bio02(tasmax_lin, tasmin_lin),  ref_linear[["bio02"]], tolerance = tol)
  expect_equal(bio03(tasmax_lin, tasmin_lin),  ref_linear[["bio03"]], tolerance = tol_loose)
  expect_equal(bio04(tas_lin),                ref_linear[["bio04"]], tolerance = tol_loose)
  expect_equal(bio05(tasmax_lin),             ref_linear[["bio05"]], tolerance = tol)
  expect_equal(bio06(tasmin_lin),             ref_linear[["bio06"]], tolerance = tol)
  expect_equal(bio07(tasmax_lin, tasmin_lin),  ref_linear[["bio07"]], tolerance = tol)
  expect_equal(bio08(tas_lin, pr_lin),        ref_linear[["bio08"]], tolerance = tol)
  expect_equal(bio09(tas_lin, pr_lin),        ref_linear[["bio09"]], tolerance = tol)
  expect_equal(bio10(tas_lin),                ref_linear[["bio10"]], tolerance = tol)
  expect_equal(bio11(tas_lin),                ref_linear[["bio11"]], tolerance = tol)
  expect_equal(bio12(pr_lin),                 ref_linear[["bio12"]], tolerance = tol)
  expect_equal(bio13(pr_lin),                 ref_linear[["bio13"]], tolerance = tol)
  expect_equal(bio14(pr_lin),                 ref_linear[["bio14"]], tolerance = tol)
  expect_equal(bio15(pr_lin),                 ref_linear[["bio15"]], tolerance = tol_loose)
  expect_equal(bio16(pr_lin),                 ref_linear[["bio16"]], tolerance = tol)
  expect_equal(bio17(pr_lin),                 ref_linear[["bio17"]], tolerance = tol)
  expect_equal(bio18(tas_lin, pr_lin),        ref_linear[["bio18"]], tolerance = tol)
  expect_equal(bio19(tas_lin, pr_lin),        ref_linear[["bio19"]], tolerance = tol)
})
test_that("L2: bioclim() unified function matches all individual functions", {
  result <- bioclim(tas_lin, tasmax_lin, tasmin_lin, pr_lin)
  expect_length(result, 19)
  expect_named(result, paste0("bio", sprintf("%02d", 1:19)))
  for (nm in names(ref_linear)) {
    expect_equal(result[[nm]], ref_linear[[nm]],
                 tolerance = tol_loose,
                 label = paste("bioclim()", nm))
  }
})
test_that("L2: Constant climate produces expected edge values", {
  expect_equal(bio01(tas_const), 20.0, tolerance = tol)
  expect_equal(bio04(tas_const), 0.0, tolerance = tol)            # zero seasonality
  expect_equal(bio07(tasmax_const, tasmin_const), 10.0, tolerance = tol)
  expect_equal(bio03(tasmax_const, tasmin_const), 100.0, tolerance = tol) # isothermal
  expect_equal(bio15(pr_const), 0.0, tolerance = tol)             # zero CV
})
test_that("L2: Anti-correlated precip flips BIO08/BIO09", {
  expect_equal(bio08(tas_lin, pr_rev), 2.0, tolerance = tol)
  expect_equal(bio09(tas_lin, pr_rev), 11.0, tolerance = tol)
})

# ============================================================================
# LAYER 3: BioclimData S4 class
# ============================================================================
test_that("L3: BioclimData constructor creates valid S4 object", {
  bd <- BioclimData(tas_lin, tasmax_lin, tasmin_lin, pr_lin)
  expect_s4_class(bd, "BioclimData")
  expect_true(is.matrix(bd@tas))
  expect_equal(nrow(bd@tas), 1L)
  expect_equal(ncol(bd@tas), 12L)
})
test_that("L3: BioclimData dispatches all 19 bio functions via C++", {
  bd <- BioclimData(tas_lin, tasmax_lin, tasmin_lin, pr_lin)
  expect_equal(bio01(bd), ref_linear[["bio01"]], tolerance = tol)
  expect_equal(bio05(bd), ref_linear[["bio05"]], tolerance = tol)
  expect_equal(bio12(bd), ref_linear[["bio12"]], tolerance = tol)
  expect_equal(bio16(bd), ref_linear[["bio16"]], tolerance = tol)
})
test_that("L3: BioclimData bioclim() returns matrix with all 19 cols", {
  bd <- BioclimData(tas_lin, tasmax_lin, tasmin_lin, pr_lin)
  result <- bioclim(bd)
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 19L)
  expect_equal(nrow(result), 1L)
})

# ============================================================================
# LAYER 4: BioclimModel (XPtr-based)
# ============================================================================
test_that("L4: BioclimModel creates valid object with external pointer", {
  m <- BioclimModel(tas_lin, tasmax_lin, tasmin_lin, pr_lin)
  expect_s4_class(m, "BioclimModel")
  expect_true(is(m@pntr, "externalptr"))
  expect_false(bioclim_model_is_null(m@pntr))
})
test_that("L4: BioclimModel matches R bioclim() numerically", {
  m <- BioclimModel(tas_lin, tasmax_lin, tasmin_lin, pr_lin)
  r_result <- bioclim(tas_lin, tasmax_lin, tasmin_lin, pr_lin)
  expect_equal(bio01(m), r_result[["bio01"]], tolerance = tol)
  expect_equal(bio12(m), r_result[["bio12"]], tolerance = tol)
  expect_equal(bio19(m), r_result[["bio19"]], tolerance = tol)
})
test_that("L4: BioclimModel rejects invalid input lengths", {
  expect_error(BioclimModel(1:6, tasmax_lin, tasmin_lin, pr_lin), "must have length 12")
  expect_error(BioclimModel(tas_lin, 1:6, tasmin_lin, pr_lin), "must have length 12")
})

# ============================================================================
# LAYER 5: ClimateBlock Rcpp module
# ============================================================================
test_that("L5: ClimateBlock creates and computes correctly", {
  block <- new(ClimateBlock,
               to_mat(tas_lin), to_mat(tasmax_lin),
               to_mat(tasmin_lin), to_mat(pr_lin))
  expect_true(is(block, "Rcpp_ClimateBlock"))
  expect_equal(block$n_pixels(), 1L)
})
test_that("L5: ClimateBlock$compute() returns 19-col named matrix", {
  block <- new(ClimateBlock,
               to_mat(tas_lin), to_mat(tasmax_lin),
               to_mat(tasmin_lin), to_mat(pr_lin))
  result <- block$compute()
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 19L)
  expect_equal(colnames(result), paste0("bio", sprintf("%02d", 1:19)))
})
test_that("L5: ClimateBlock matches bioclim() for all 19 vars", {
  block <- new(ClimateBlock,
               to_mat(tas_lin), to_mat(tasmax_lin),
               to_mat(tasmin_lin), to_mat(pr_lin))
  cpp_res <- block$compute()[1L, ]
  r_res   <- bioclim(tas_lin, tasmax_lin, tasmin_lin, pr_lin)
  for (i in seq_along(r_res)) {
    if (is.nan(r_res[i])) {
      expect_true(is.nan(cpp_res[i]))
    } else {
      expect_equal(cpp_res[i], r_res[i], tolerance = tol,
                   label = paste("ClimateBlock var", i))
    }
  }
})
test_that("L5: ClimateBlock multi-pixel consistency", {
  n <- 50
  mat_tas    <- matrix(rep(tas_lin, n), nrow = n, byrow = TRUE)
  mat_tasmax <- matrix(rep(tasmax_lin, n), nrow = n, byrow = TRUE)
  mat_tasmin <- matrix(rep(tasmin_lin, n), nrow = n, byrow = TRUE)
  mat_pr     <- matrix(rep(pr_lin, n), nrow = n, byrow = TRUE)
  
  block  <- new(ClimateBlock, mat_tas, mat_tasmax, mat_tasmin, mat_pr)
  result <- block$compute()
  expect_equal(nrow(result), n)
  # All rows should be identical since input is replicated
  for (i in 2:n) {
    expect_equal(result[i, ], result[1, ], tolerance = 1e-12,
                 label = paste("row", i, "matches row 1"))
  }
})

# ============================================================================
# LAYER 6: bioclim_xt() — zero-copy vectorized C++ path
# ============================================================================
test_that("L6: bioclim_xt returns n_pixels x 19 matrix", {
  result <- bioclim_xt(to_mat(tas_lin), to_mat(tasmax_lin),
                       to_mat(tasmin_lin), to_mat(pr_lin))
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1L)
  expect_equal(ncol(result), 19L)
})
test_that("L6: bioclim_xt matches reference values for linear pattern", {
  result <- bioclim_xt(to_mat(tas_lin), to_mat(tasmax_lin),
                       to_mat(tasmin_lin), to_mat(pr_lin))
  for (j in 1:19) {
    nm <- paste0("bio", sprintf("%02d", j))
    expect_equal(unname(result[1, j]), ref_linear[[nm]], tolerance = tol_loose,
                 label = paste("bioclim_xt", nm))
  }
})
test_that("L6: bioclim_xt handles multi-pixel batch correctly", {
  n <- 100
  mat_tas    <- matrix(rep(tas_lin, n), nrow = n, byrow = TRUE)
  mat_tasmax <- matrix(rep(tasmax_lin, n), nrow = n, byrow = TRUE)
  mat_tasmin <- matrix(rep(tasmin_lin, n), nrow = n, byrow = TRUE)
  mat_pr     <- matrix(rep(pr_lin, n), nrow = n, byrow = TRUE)
  
  result <- bioclim_xt(mat_tas, mat_tasmax, mat_tasmin, mat_pr)
  expect_equal(nrow(result), n)
  expect_equal(ncol(result), 19L)
  # All rows identical
  for (i in 2:n) {
    expect_equal(result[i, ], result[1, ], tolerance = 1e-12)
  }
})
test_that("L6: bioclim_xt with OpenMP (ncores=2) matches single-threaded", {
  n <- 200
  mat_tas    <- matrix(rep(tas_lin, n), nrow = n, byrow = TRUE)
  mat_tasmax <- matrix(rep(tasmax_lin, n), nrow = n, byrow = TRUE)
  mat_tasmin <- matrix(rep(tasmin_lin, n), nrow = n, byrow = TRUE)
  mat_pr     <- matrix(rep(pr_lin, n), nrow = n, byrow = TRUE)
  
  res1 <- bioclim_xt(mat_tas, mat_tasmax, mat_tasmin, mat_pr, ncores = 1L)
  res2 <- bioclim_xt(mat_tas, mat_tasmax, mat_tasmin, mat_pr, ncores = 2L)
  expect_equal(res1, res2, tolerance = 1e-12)
})
test_that("L6: bioclim_xt with realistic climate data", {
  result <- bioclim_xt(to_mat(tas_real), to_mat(tasmax_real),
                       to_mat(tasmin_real), to_mat(pr_real))
  # BIO01 = mean of monthly mean temp
  expect_equal(unname(result[1, 1]), mean(tas_real), tolerance = tol)
  # BIO12 = sum of monthly precip
  expect_equal(unname(result[1, 12]), sum(pr_real), tolerance = tol)
  # BIO05 = max of monthly max temp
  expect_equal(unname(result[1, 5]), max(tasmax_real), tolerance = tol)
  # BIO06 = min of monthly min temp
  expect_equal(unname(result[1, 6]), min(tasmin_real), tolerance = tol)
})

# ============================================================================
# LAYER 7: bioclim_raster() — terra block-loop path
# ============================================================================
test_that("L7: bioclim_raster produces 19-layer SpatRaster", {
  skip_without_terra()
  n_cells <- 4 * 3
  make_rast <- function(vals) {
    r <- terra::rast(nrows = 4, ncols = 3, nlyr = 12,
                     xmin = 0, xmax = 3, ymin = 0, ymax = 4, crs = "EPSG:4326")
    terra::values(r) <- matrix(rep(vals, n_cells), nrow = n_cells, ncol = 12, byrow = TRUE)
    r
  }
  result <- bioclim_raster(make_rast(tas_lin), make_rast(tasmax_lin),
                           make_rast(tasmin_lin), make_rast(pr_lin))
  expect_equal(terra::nlyr(result), 19L)
  expect_equal(names(result), paste0("bio", sprintf("%02d", 1:19)))
})
test_that("L7: bioclim_raster pixel values match bioclim()", {
  skip_without_terra()
  n_cells <- 4 * 3
  make_rast <- function(vals) {
    r <- terra::rast(nrows = 4, ncols = 3, nlyr = 12,
                     xmin = 0, xmax = 3, ymin = 0, ymax = 4, crs = "EPSG:4326")
    terra::values(r) <- matrix(rep(vals, n_cells), nrow = n_cells, ncol = 12, byrow = TRUE)
    r
  }
  result <- bioclim_raster(make_rast(tas_lin), make_rast(tasmax_lin),
                           make_rast(tasmin_lin), make_rast(pr_lin))
  vals <- terra::values(result)
  # Every cell should match the reference
  for (cell in 1:n_cells) {
    for (j in 1:19) {
      nm <- paste0("bio", sprintf("%02d", j))
      expect_equal(as.numeric(vals[cell, j]), ref_linear[[nm]], tolerance = tol_loose,
                   label = paste("raster cell", cell, nm))
    }
  }
})
test_that("L7: bioclim_raster handles NA cells", {
  skip_without_terra()
  n_cells <- 4 * 3
  make_rast <- function(vals) {
    r <- terra::rast(nrows = 4, ncols = 3, nlyr = 12,
                     xmin = 0, xmax = 3, ymin = 0, ymax = 4, crs = "EPSG:4326")
    terra::values(r) <- matrix(rep(vals, n_cells), nrow = n_cells, ncol = 12, byrow = TRUE)
    r
  }
  r_tas <- make_rast(tas_lin)
  # Set cell 1, layer 1 to NA
  terra::values(r_tas)[1, 1] <- NA
  result <- bioclim_raster(r_tas, make_rast(tasmax_lin),
                           make_rast(tasmin_lin), make_rast(pr_lin))
  vals <- terra::values(result)
  # Cell 1 should be all NA
  expect_true(all(is.na(vals[1, ])))
  # Cell 2 should be valid
  expect_false(any(is.na(vals[2, ])))
})
test_that("L7: bioclim_raster with multiple blocks matches single block", {
  skip_without_terra()
  n_cells <- 4 * 3
  make_rast <- function(vals) {
    r <- terra::rast(nrows = 4, ncols = 3, nlyr = 12,
                     xmin = 0, xmax = 3, ymin = 0, ymax = 4, crs = "EPSG:4326")
    terra::values(r) <- matrix(rep(vals, n_cells), nrow = n_cells, ncol = 12, byrow = TRUE)
    r
  }
  r_t <- make_rast(tas_lin); r_tx <- make_rast(tasmax_lin)
  r_tn <- make_rast(tasmin_lin); r_p <- make_rast(pr_lin)

  res1 <- bioclim_raster(r_t, r_tx, r_tn, r_p, n_blocks = 1L)
  res4 <- bioclim_raster(r_t, r_tx, r_tn, r_p, n_blocks = 4L)
  expect_equal(terra::values(res1), terra::values(res4), tolerance = 1e-12)
})
test_that("L7: bioclim_raster with ncores=2 matches ncores=1", {
  skip_without_terra()
  n_cells <- 4 * 3
  make_rast <- function(vals) {
    r <- terra::rast(nrows = 4, ncols = 3, nlyr = 12,
                     xmin = 0, xmax = 3, ymin = 0, ymax = 4, crs = "EPSG:4326")
    terra::values(r) <- matrix(rep(vals, n_cells), nrow = n_cells, ncol = 12, byrow = TRUE)
    r
  }
  r_t <- make_rast(tas_lin); r_tx <- make_rast(tasmax_lin)
  r_tn <- make_rast(tasmin_lin); r_p <- make_rast(pr_lin)

  res1 <- bioclim_raster(r_t, r_tx, r_tn, r_p, ncores = 1L)
  res2 <- bioclim_raster(r_t, r_tx, r_tn, r_p, ncores = 2L)
  expect_equal(terra::values(res1), terra::values(res2), tolerance = 1e-12)
})
test_that("L7: bioclim_raster rejects non-SpatRaster", {
  skip_without_terra()
  expect_error(bioclim_raster("not_a_raster", "a", "b", "c"), "must be a SpatRaster")
})

# ============================================================================
# LAYER 8: BioclimEngine — GDAL tiled XPtr engine
# ============================================================================
test_that("L8: engine_create returns externalptr", {
  ptr <- engine_create()
  expect_true(is(ptr, "externalptr"))
})
test_that("L8: all engine setter methods work without error", {
  ptr <- engine_create()
  expect_no_error(engine_set_output(ptr, tempfile(fileext = ".tif")))
  expect_no_error(engine_set_mask(ptr, ""))
  expect_no_error(engine_set_threads(ptr, 4L))
  expect_no_error(engine_set_tile_size(ptr, 64L))
})
test_that("L8: engine_compute stops clearly without GDAL", {
  if (has_gdal()) skip("GDAL present")
  ptr <- engine_create()
  expect_error(engine_compute(ptr), "GDAL")
})
test_that("L8: engine round-trip with tiny rasters", {
  skip_without_gdal()
  skip_without_terra()

  tmpdir <- tempfile("engine_rt_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  make_tif <- function(vals, path) {
    r <- terra::rast(nrows = 3, ncols = 3, nlyrs = 12,
                     xmin = 0, xmax = 1, ymin = 0, ymax = 1, crs = "EPSG:4326")
    for (m in 1:12) terra::values(r[[m]]) <- vals[m]
    terra::writeRaster(r, path, overwrite = TRUE)
    path
  }

  tas_f    <- make_tif(tas_real,    file.path(tmpdir, "tas.tif"))
  tasmax_f <- make_tif(tasmax_real, file.path(tmpdir, "tasmax.tif"))
  tasmin_f <- make_tif(tasmin_real, file.path(tmpdir, "tasmin.tif"))
  pr_f     <- make_tif(pr_real,     file.path(tmpdir, "pr.tif"))
  out_f    <- file.path(tmpdir, "bio_out.tif")

  ptr <- engine_create()
  engine_open(ptr, tas_f, tasmax_f, tasmin_f, pr_f)
  engine_set_output(ptr, out_f)
  engine_set_threads(ptr, 1L)
  engine_set_tile_size(ptr, 2L)  # force edge-tile handling
  result_path <- engine_compute(ptr)

  expect_true(file.exists(result_path))
  out_r <- terra::rast(result_path)
  expect_equal(terra::nlyr(out_r), 19L)

  vals <- terra::values(out_r)
  # BIO01 = mean(tas_real)
  expect_equal(unname(vals[1, 1]), mean(tas_real), tolerance = tol)
  # BIO12 = sum(pr_real)
  expect_equal(unname(vals[1, 12]), sum(pr_real), tolerance = tol)
  # BIO05 = max(tasmax_real)
  expect_equal(unname(vals[1, 5]), max(tasmax_real), tolerance = tol)
})

# ============================================================================
# LAYER 9: bioclim_engine() — R user-facing API
# ============================================================================
test_that("L9: bioclim_engine() stops without GDAL", {
  if (has_gdal()) skip("GDAL present")
  skip_without_terra()
  tmpdir <- tempfile("be_nogdal_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
  # Create dummy files so validation doesn't fail first
  f <- file.path(tmpdir, "dummy.tif")
  r <- terra::rast(nrows = 3, ncols = 3, nlyrs = 12)
  terra::values(r) <- 1
  terra::writeRaster(r, f, overwrite = TRUE)
  expect_error(bioclim_engine(f, f, f, f), "GDAL")
})
test_that("L9: bioclim_engine() input validation", {
  skip_without_gdal()
  # Bad output
  expect_error(bioclim_engine("a", "b", "c", "d", output = 42), "character string")
  # Bad threads
  expect_error(bioclim_engine("a", "b", "c", "d", threads = -1), "integer >= 1")
  # Bad tile_size
  expect_error(bioclim_engine("a", "b", "c", "d", tile_size = 0), "integer >= 1")
})
test_that("L9: bioclim_engine() full round-trip with multi-band files", {
  skip_without_gdal()
  skip_without_terra()

  tmpdir <- tempfile("be_round_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  make_tif <- function(vals, path) {
    r <- terra::rast(nrows = 5, ncols = 5, nlyrs = 12,
                     xmin = 0, xmax = 1, ymin = 0, ymax = 1, crs = "EPSG:4326")
    for (m in 1:12) terra::values(r[[m]]) <- vals[m]
    terra::writeRaster(r, path, overwrite = TRUE)
    path
  }

  tas_f    <- make_tif(tas_real,    file.path(tmpdir, "tas.tif"))
  tasmax_f <- make_tif(tasmax_real, file.path(tmpdir, "tasmax.tif"))
  tasmin_f <- make_tif(tasmin_real, file.path(tmpdir, "tasmin.tif"))
  pr_f     <- make_tif(pr_real,     file.path(tmpdir, "pr.tif"))

  result <- bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f,
                           threads = 2L, tile_size = 3L)

  if (inherits(result, "SpatRaster")) {
    expect_equal(terra::nlyr(result), 19L)
    vals <- terra::values(result)
    # All 25 pixels should have the same values (uniform input)
    expect_equal(unname(vals[1, 1]), mean(tas_real), tolerance = tol)
    expect_equal(unname(vals[1, 12]), sum(pr_real), tolerance = tol)
    # Check consistency: all cells identical
    for (i in 2:25) {
      expect_equal(vals[i, ], vals[1, ], tolerance = 1e-10,
                   label = paste("cell", i))
    }
  } else {
    # Without terra, we get a character vector of file paths
    expect_true(all(file.exists(result)))
  }
})
test_that("L9: bioclim_engine() accepts SpatRaster input", {
  skip_without_gdal()
  skip_without_terra()

  tmpdir <- tempfile("be_spatr_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  make_rast <- function(vals) {
    r <- terra::rast(nrows = 3, ncols = 3, nlyrs = 12,
                     xmin = 0, xmax = 1, ymin = 0, ymax = 1, crs = "EPSG:4326")
    for (m in 1:12) terra::values(r[[m]]) <- vals[m]
    f <- file.path(tmpdir, paste0(c(sample(letters, 8, replace = TRUE), ".tif"), collapse = ""))
    terra::writeRaster(r, f, overwrite = TRUE)
    terra::rast(f)
  }

  result <- bioclim_engine(make_rast(tas_real), make_rast(tasmax_real),
                           make_rast(tasmin_real), make_rast(pr_real))
  if (inherits(result, "SpatRaster")) {
    expect_equal(terra::nlyr(result), 19L)
  }
})
test_that("L9: bioclim_engine() overwrite behavior", {
  skip_without_gdal()
  skip_without_terra()

  tmpdir <- tempfile("be_overwrite_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  make_tif <- function(vals, path) {
    r <- terra::rast(nrows = 3, ncols = 3, nlyrs = 12,
                     xmin = 0, xmax = 1, ymin = 0, ymax = 1, crs = "EPSG:4326")
    for (m in 1:12) terra::values(r[[m]]) <- vals[m]
    terra::writeRaster(r, path, overwrite = TRUE)
    path
  }
  tas_f    <- make_tif(tas_real,    file.path(tmpdir, "tas.tif"))
  tasmax_f <- make_tif(tasmax_real, file.path(tmpdir, "tasmax.tif"))
  tasmin_f <- make_tif(tasmin_real,    file.path(tmpdir, "tasmin.tif"))
  pr_f     <- make_tif(pr_real,     file.path(tmpdir, "pr.tif"))
  out_d    <- file.path(tmpdir, "out_dir")

  bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f, output = out_d)
  expect_true(dir.exists(out_d))
  # Second call without overwrite should fail (files already exist)

  expect_error(bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f, output = out_d),
               "already exist")
  # With overwrite = TRUE should succeed
  expect_no_error(
    bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f, output = out_d, overwrite = TRUE)
  )
})

# ============================================================================
# LAYER 10: Masking support
# ============================================================================
test_that("L10: create_mask with character path", {
  skip_without_gdal()

  ref <- system.file("extdata", "tiny.tif", package = "rxbioclim")
  if (!nzchar(ref)) skip("tiny.tif not found in inst/extdata")

  poly <- tempfile(fileext = ".geojson")
  mask <- tempfile(fileext = ".tif")
  on.exit(unlink(c(poly, mask)), add = TRUE)

  # Write covering polygon
  geojson <- paste0(
    '{"type":"FeatureCollection","features":[{"type":"Feature",
    "geometry":{"type":"Polygon",
    "coordinates":[[[0,0],[3,0],[3,-3],[0,-3],[0,0]]]},
    "properties":{}}]}'
  )
  writeLines(geojson, poly)
  rasterize_mask_cpp(poly, ref, mask)
  expect_true(file.exists(mask))
})
test_that("L10: bioclim_engine with mask parameter", {
  skip_without_gdal()
  skip_without_terra()

  tmpdir <- tempfile("mask_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  make_tif <- function(vals, path) {
    r <- terra::rast(nrows = 4, ncols = 4, nlyrs = 12,
                     xmin = 0, xmax = 1, ymin = 0, ymax = 1, crs = "EPSG:4326")
    for (m in 1:12) terra::values(r[[m]]) <- vals[m]
    terra::writeRaster(r, path, overwrite = TRUE)
    path
  }
  # Create mask: top half = 1, bottom half = 0
  mask_r <- terra::rast(nrows = 4, ncols = 4, nlyrs = 1,
                        xmin = 0, xmax = 1, ymin = 0, ymax = 1, crs = "EPSG:4326")
  mask_vals <- c(rep(1, 8), rep(0, 8))  # 4x4, top 2 rows = 1, bottom 2 rows = 0
  terra::values(mask_r) <- mask_vals
  mask_f <- file.path(tmpdir, "mask.tif")
  terra::writeRaster(mask_r, mask_f, overwrite = TRUE)

  tas_f    <- make_tif(tas_real,    file.path(tmpdir, "tas.tif"))
  tasmax_f <- make_tif(tasmax_real, file.path(tmpdir, "tasmax.tif"))
  tasmin_f <- make_tif(tasmin_real, file.path(tmpdir, "tasmin.tif"))
  pr_f     <- make_tif(pr_real,     file.path(tmpdir, "pr.tif"))

  result <- bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f,
                           mask = mask_f, tile_size = 2L)
  if (inherits(result, "SpatRaster")) {
    vals <- terra::values(result)
    # Top half (rows 1-8) should have real values
    expect_false(any(is.na(vals[1, ])))
    # Bottom half (rows 9-16) should be NA/NaN (masked)
    expect_true(all(is.na(vals[9, ]) | is.nan(vals[9, ])))
  }
})

# ============================================================================
# LAYER 11: Cross-layer numerical consistency
# ============================================================================
test_that("L11: R bioclim() == bioclim_xt() == ClimateBlock == BioclimData", {
  r_result  <- bioclim(tas_lin, tasmax_lin, tasmin_lin, pr_lin)
  xt_result <- bioclim_xt(to_mat(tas_lin), to_mat(tasmax_lin),
                          to_mat(tasmin_lin), to_mat(pr_lin))[1, ]
  cb_result <- new(ClimateBlock, to_mat(tas_lin), to_mat(tasmax_lin),
                   to_mat(tasmin_lin), to_mat(pr_lin))$compute()[1, ]
  bd_result <- bioclim(BioclimData(tas_lin, tasmax_lin, tasmin_lin, pr_lin))[1, ]

  for (j in 1:19) {
    nm <- paste0("bio", sprintf("%02d", j))
    vals <- c(R = r_result[[nm]], xt = unname(xt_result[j]),
              CB = unname(cb_result[j]), BD = unname(bd_result[j]))

    if (is.nan(vals["R"])) {
      for (v in vals) expect_true(is.nan(v), label = paste(nm, "NaN consistency"))
    } else {
      expect_equal(unname(vals["xt"]), unname(vals["R"]), tolerance = 1e-10,
                   label = paste(nm, "xt vs R"))
      expect_equal(unname(vals["CB"]), unname(vals["R"]), tolerance = 1e-10,
                   label = paste(nm, "CB vs R"))
      expect_equal(unname(vals["BD"]), unname(vals["R"]), tolerance = 1e-10,
                   label = paste(nm, "BD vs R"))
    }
  }
})
test_that("L11: bioclim_raster matches bioclim_engine for same input", {
  skip_without_gdal()
  skip_without_terra()

  tmpdir <- tempfile("cross_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  n_cells <- 4 * 3
  make_rast <- function(vals) {
    r <- terra::rast(nrows = 4, ncols = 3, nlyr = 12,
                     xmin = 0, xmax = 3, ymin = 0, ymax = 4, crs = "EPSG:4326")
    terra::values(r) <- matrix(rep(vals, n_cells), nrow = n_cells, ncol = 12, byrow = TRUE)
    r
  }

  make_tif <- function(vals, path) {
    r <- make_rast(vals)
    terra::writeRaster(r, path, overwrite = TRUE)
    path
  }

  r_t  <- make_rast(tas_real);  r_tx <- make_rast(tasmax_real)
  r_tn <- make_rast(tasmin_real); r_p <- make_rast(pr_real)

  raster_result <- bioclim_raster(r_t, r_tx, r_tn, r_p)
  raster_vals   <- terra::values(raster_result)

  tas_f    <- make_tif(tas_real,    file.path(tmpdir, "tas.tif"))
  tasmax_f <- make_tif(tasmax_real, file.path(tmpdir, "tasmax.tif"))
  tasmin_f <- make_tif(tasmin_real, file.path(tmpdir, "tasmin.tif"))
  pr_f     <- make_tif(pr_real,     file.path(tmpdir, "pr.tif"))

  engine_result <- bioclim_engine(tas_f, tasmax_f, tasmin_f, pr_f, tile_size = 2L)
  engine_vals   <- terra::values(engine_result)

  # Both should produce 19 columns, same number of rows
  expect_equal(ncol(raster_vals), 19L)
  expect_equal(ncol(engine_vals), 19L)
  expect_equal(nrow(raster_vals), nrow(engine_vals))

  # Values should match cell-by-cell
  for (i in 1:nrow(raster_vals)) {
    for (j in 1:19) {
      if (is.na(raster_vals[i, j]) || is.nan(raster_vals[i, j])) {
        expect_true(is.na(engine_vals[i, j]) || is.nan(engine_vals[i, j]),
                    label = paste("cell", i, "var", j, "NA consistency"))
      } else {
        expect_equal(unname(engine_vals[i, j]), unname(raster_vals[i, j]), tolerance = 1e-6,
                     label = paste("cell", i, "var", j))
      }
    }
  }
})

# ============================================================================
# LAYER 12: Edge cases, NA handling, stress tests
# ============================================================================
test_that("L12: bioclim_xt all-NA row produces all-NA output", {
  mat <- to_mat(tas_lin)
  mat[1, 1] <- NA
  result <- bioclim_xt(mat, to_mat(tasmax_lin), to_mat(tasmin_lin), to_mat(pr_lin))
  expect_true(all(is.na(result[1, ])))
})
test_that("L12: bioclim_xt mixed NA and valid rows", {
  n <- 10
  mat_tas    <- matrix(rep(tas_lin, n), nrow = n, byrow = TRUE)
  mat_tasmax <- matrix(rep(tasmax_lin, n), nrow = n, byrow = TRUE)
  mat_tasmin <- matrix(rep(tasmin_lin, n), nrow = n, byrow = TRUE)
  mat_pr     <- matrix(rep(pr_lin, n), nrow = n, byrow = TRUE)
  # Inject NA in row 3 and row 7
  mat_tas[3, 5] <- NA
  mat_pr[7, 1]  <- NA

  result <- bioclim_xt(mat_tas, mat_tasmax, mat_tasmin, mat_pr)
  expect_true(all(is.na(result[3, ])))
  expect_true(all(is.na(result[7, ])))
  # Other rows should be valid
  expect_false(any(is.na(result[1, ])))
  expect_false(any(is.na(result[5, ])))
})
test_that("L12: bioclim() with zero precipitation is handled", {
  result <- bioclim(tas_lin, tasmax_lin, tasmin_lin, rep(0, 12))
  expect_equal(result[["bio12"]], 0)
  expect_equal(result[["bio13"]], 0)
  expect_equal(result[["bio14"]], 0)
  expect_equal(result[["bio16"]], 0)
  expect_equal(result[["bio17"]], 0)
})
test_that("L12: bioclim() with identical temperatures across months", {
  result <- bioclim(tas_const, tasmax_const, tasmin_const, pr_lin)
  expect_equal(result[["bio01"]], 20.0, tolerance = tol)
  expect_equal(result[["bio04"]], 0.0, tolerance = tol)   # zero seasonality
  expect_equal(result[["bio05"]], 25.0, tolerance = tol)
  expect_equal(result[["bio06"]], 15.0, tolerance = tol)
  expect_equal(result[["bio07"]], 10.0, tolerance = tol)
  # When all quarters have the same temp, warmest/coldest quarter = mean
  expect_equal(result[["bio10"]], 20.0, tolerance = tol)
  expect_equal(result[["bio11"]], 20.0, tolerance = tol)
})
test_that("L12: bioclim_xt handles large batch (1000 pixels)", {
  n <- 1000
  mat_tas    <- matrix(rep(tas_real, n), nrow = n, byrow = TRUE)
  mat_tasmax <- matrix(rep(tasmax_real, n), nrow = n, byrow = TRUE)
  mat_tasmin <- matrix(rep(tasmin_real, n), nrow = n, byrow = TRUE)
  mat_pr     <- matrix(rep(pr_real, n), nrow = n, byrow = TRUE)

  result <- bioclim_xt(mat_tas, mat_tasmax, mat_tasmin, mat_pr, ncores = 2L)
  expect_equal(nrow(result), n)
  expect_equal(ncol(result), 19L)
  # Spot check: all rows should be identical
  expect_equal(result[1, ], result[500, ], tolerance = 1e-12)
  expect_equal(result[1, ], result[1000, ], tolerance = 1e-12)
})
test_that("L12: has_gdal returns logical", {
  expect_type(has_gdal(), "logical")
})
test_that("L12: message system works end-to-end", {
  clear_messages()
  expect_false(has_error())
  expect_false(has_warning())
  expect_equal(bioclim_errors(), character(0))
  expect_equal(bioclim_warnings(), character(0))
})
