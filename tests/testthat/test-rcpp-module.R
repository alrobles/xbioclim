# Tests for the Rcpp bioclim_mod module (ClimateBlock class)
#
# The module is loaded automatically when rxbioclim is attached.
# These tests exercise the C++ path and verify numerical equivalence
# with the pure-R bioclim() function.

# ── Shared test fixtures ──────────────────────────────────────────────────────

# Mock data: same convention as test-bioclim.R
mock_tas    <- 1:12
mock_tasmax <- 2:13
mock_tasmin <- 0:11
mock_pr     <- 1:12

to_mat <- function(x) matrix(as.double(x), nrow = 1L, ncol = 12L)

tol <- 1e-4

# ── ClimateBlock construction ─────────────────────────────────────────────────

test_that("ClimateBlock can be constructed from 1×12 matrices", {
  block <- new(ClimateBlock,
    to_mat(mock_tas), to_mat(mock_tasmax),
    to_mat(mock_tasmin), to_mat(mock_pr))
  expect_true(is(block, "Rcpp_ClimateBlock"))
  expect_equal(block$n_pixels(), 1L)
})

test_that("ClimateBlock constructor rejects wrong number of columns", {
  bad <- matrix(1:6, nrow = 1L, ncol = 6L)
  expect_error(
    new(ClimateBlock, bad, to_mat(mock_tasmax),
      to_mat(mock_tasmin), to_mat(mock_pr)),
    "12 columns"
  )
})

test_that("ClimateBlock constructor rejects mismatched row counts", {
  two_rows <- matrix(rep(mock_tas, 2), nrow = 2L, ncol = 12L)
  expect_error(
    new(ClimateBlock, two_rows, to_mat(mock_tasmax),
      to_mat(mock_tasmin), to_mat(mock_pr)),
    "same dimensions"
  )
})

# ── ClimateBlock$compute() — single pixel ────────────────────────────────────

test_that("ClimateBlock$compute() returns a 1×19 matrix with named columns", {
  block  <- new(ClimateBlock,
    to_mat(mock_tas), to_mat(mock_tasmax),
    to_mat(mock_tasmin), to_mat(mock_pr))
  result <- block$compute()

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1L)
  expect_equal(ncol(result), 19L)
  expect_equal(colnames(result), paste0("bio", sprintf("%02d", 1:19)))
})

test_that("ClimateBlock$compute() matches pure-R bioclim() for all 19 vars", {
  block   <- new(ClimateBlock,
    to_mat(mock_tas), to_mat(mock_tasmax),
    to_mat(mock_tasmin), to_mat(mock_pr))
  cpp_res <- block$compute()[1L, ]
  r_res   <- bioclim(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  # Compare element-by-element (NaN-safe via is.nan check)
  for (i in seq_along(r_res)) {
    if (is.nan(r_res[i])) {
      expect_true(is.nan(cpp_res[i]))
    } else {
      expect_equal(cpp_res[i], r_res[i], tolerance = tol)
    }
  }
})

test_that("ClimateBlock$compute() exact reference values (single pixel)", {
  block <- new(ClimateBlock,
    to_mat(mock_tas), to_mat(mock_tasmax),
    to_mat(mock_tasmin), to_mat(mock_pr))
  res   <- block$compute()[1L, ]

  expect_equal(res[["bio01"]],  6.5,      tolerance = tol)
  expect_equal(res[["bio02"]],  2.0,      tolerance = tol)
  expect_equal(res[["bio03"]],  15.3846,  tolerance = 1e-3)
  expect_equal(res[["bio04"]],  345.2053, tolerance = 1e-3)
  expect_equal(res[["bio05"]],  13.0,     tolerance = tol)
  expect_equal(res[["bio06"]],  0.0,      tolerance = tol)
  expect_equal(res[["bio07"]],  13.0,     tolerance = tol)
  expect_equal(res[["bio08"]],  11.0,     tolerance = tol)
  expect_equal(res[["bio09"]],  2.0,      tolerance = tol)
  expect_equal(res[["bio10"]],  11.0,     tolerance = tol)
  expect_equal(res[["bio11"]],  2.0,      tolerance = tol)
  expect_equal(res[["bio12"]],  78.0,     tolerance = tol)
  expect_equal(res[["bio13"]],  12.0,     tolerance = tol)
  expect_equal(res[["bio14"]],  1.0,      tolerance = tol)
  expect_equal(res[["bio15"]],  53.1085,  tolerance = 1e-3)
  expect_equal(res[["bio16"]],  33.0,     tolerance = tol)
  expect_equal(res[["bio17"]],  6.0,      tolerance = tol)
  expect_equal(res[["bio18"]],  33.0,     tolerance = tol)
  expect_equal(res[["bio19"]],  6.0,      tolerance = tol)
})

test_that("ClimateBlock$compute() returns NaN for bio03 when annual range is 0", {
  tasmax_eq <- rep(10.0, 12)
  tasmin_eq <- rep(10.0, 12)
  block <- new(ClimateBlock,
    to_mat(rep(10.0, 12)), to_mat(tasmax_eq),
    to_mat(tasmin_eq), to_mat(mock_pr))
  res   <- block$compute()[1L, ]
  expect_true(is.nan(res[["bio03"]]))
})

test_that("ClimateBlock$compute() returns NaN for bio15 when pr is all zero", {
  block <- new(ClimateBlock,
    to_mat(mock_tas), to_mat(mock_tasmax),
    to_mat(mock_tasmin), to_mat(rep(0.0, 12)))
  res   <- block$compute()[1L, ]
  expect_true(is.nan(res[["bio15"]]))
})

# ── ClimateBlock$compute() — multi-pixel ─────────────────────────────────────

test_that("ClimateBlock supports multiple pixels (n × 12 matrices)", {
  # Stack two copies of the same pixel
  tas2    <- matrix(rep(as.double(mock_tas),    2L), nrow = 2L, byrow = TRUE)
  tasmax2 <- matrix(rep(as.double(mock_tasmax), 2L), nrow = 2L, byrow = TRUE)
  tasmin2 <- matrix(rep(as.double(mock_tasmin), 2L), nrow = 2L, byrow = TRUE)
  pr2     <- matrix(rep(as.double(mock_pr),     2L), nrow = 2L, byrow = TRUE)

  block  <- new(ClimateBlock, tas2, tasmax2, tasmin2, pr2)
  expect_equal(block$n_pixels(), 2L)

  result <- block$compute()
  expect_equal(nrow(result), 2L)
  expect_equal(ncol(result), 19L)

  # Both rows must be identical and match the single-pixel reference
  expect_equal(result[1L, ], result[2L, ])

  r_ref <- bioclim(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  for (i in seq_along(r_ref)) {
    if (is.nan(r_ref[i])) {
      expect_true(is.nan(result[1L, i]))
    } else {
      expect_equal(result[1L, i], r_ref[i], tolerance = tol)
    }
  }
})

# ── bioclim_block() convenience wrapper ──────────────────────────────────────

test_that("bioclim_block() returns a named vector of length 19", {
  res <- bioclim_block(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_length(res, 19L)
  expect_named(res, paste0("bio", sprintf("%02d", 1:19)))
})

test_that("bioclim_block() matches bioclim() for all 19 variables", {
  cpp_res <- bioclim_block(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  r_res   <- bioclim(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  for (i in seq_along(r_res)) {
    if (is.nan(r_res[i])) {
      expect_true(is.nan(cpp_res[i]))
    } else {
      expect_equal(cpp_res[i], r_res[i], tolerance = tol)
    }
  }
})

test_that("bioclim_block() validates input length", {
  expect_error(bioclim_block(1:6, mock_tasmax, mock_tasmin, mock_pr),
               "must have length 12")
  expect_error(bioclim_block(mock_tas, 1:6, mock_tasmin, mock_pr),
               "must have length 12")
  expect_error(bioclim_block(mock_tas, mock_tasmax, 1:6, mock_pr),
               "must have length 12")
  expect_error(bioclim_block(mock_tas, mock_tasmax, mock_tasmin, 1:6),
               "must have length 12")
})

test_that("bioclim_block() validates input type", {
  expect_error(bioclim_block(letters[1:12], mock_tasmax, mock_tasmin, mock_pr),
               "must be numeric")
})
