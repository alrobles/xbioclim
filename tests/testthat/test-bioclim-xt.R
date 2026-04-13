# Numerical equivalence tests: bioclim_xt() vs bioclim_cpp()
# These cross-implementation tests validate that the zero-copy vectorized path
# (bioclim_xt) produces bit-identical results to the reference implementation
# (bioclim_cpp) across representative input cases.

# Standard linear monthly test pattern
tas_vec    <- 1:12
tasmax_vec <- 2:13
tasmin_vec <- 0:11
pr_vec     <- 1:12

# Helper: build 1-row matrices from vectors
make_mats <- function(tas, tasmax, tasmin, pr) {
  list(
    tas    = matrix(as.double(tas),    nrow = 1L),
    tasmax = matrix(as.double(tasmax), nrow = 1L),
    tasmin = matrix(as.double(tasmin), nrow = 1L),
    pr     = matrix(as.double(pr),     nrow = 1L)
  )
}

# Helper: build n-row matrices repeating the same pattern n times
make_mats_n <- function(n, tas, tasmax, tasmin, pr) {
  list(
    tas    = matrix(rep(as.double(tas),    n), nrow = n, ncol = 12L, byrow = TRUE),
    tasmax = matrix(rep(as.double(tasmax), n), nrow = n, ncol = 12L, byrow = TRUE),
    tasmin = matrix(rep(as.double(tasmin), n), nrow = n, ncol = 12L, byrow = TRUE),
    pr     = matrix(rep(as.double(pr),     n), nrow = n, ncol = 12L, byrow = TRUE)
  )
}

test_that("bioclim_xt matches bioclim_cpp on standard linear pattern", {
  m <- make_mats(tas_vec, tasmax_vec, tasmin_vec, pr_vec)
  cpp_out <- bioclim_cpp(m$tas, m$tasmax, m$tasmin, m$pr)
  xt_out  <- bioclim_xt(m$tas,  m$tasmax, m$tasmin, m$pr)
  # Strip column names for comparison (bioclim_xt strips them too via .compute_bioclim_block)
  colnames(cpp_out) <- NULL
  colnames(xt_out)  <- NULL
  expect_equal(xt_out, cpp_out)
})

test_that("bioclim_xt matches bioclim_cpp for multiple identical pixels", {
  n <- 50L
  m <- make_mats_n(n, tas_vec, tasmax_vec, tasmin_vec, pr_vec)
  cpp_out <- bioclim_cpp(m$tas, m$tasmax, m$tasmin, m$pr)
  xt_out  <- bioclim_xt(m$tas,  m$tasmax, m$tasmin, m$pr)
  colnames(cpp_out) <- NULL
  colnames(xt_out)  <- NULL
  expect_equal(xt_out, cpp_out)
})

test_that("bioclim_xt propagates NA rows identically to bioclim_cpp", {
  n <- 5L
  m <- make_mats_n(n, tas_vec, tasmax_vec, tasmin_vec, pr_vec)
  # Inject NA into rows 2 and 4
  m$tas[2, 3]    <- NA_real_
  m$tasmax[4, 1] <- NA_real_
  cpp_out <- bioclim_cpp(m$tas, m$tasmax, m$tasmin, m$pr)
  xt_out  <- bioclim_xt(m$tas,  m$tasmax, m$tasmin, m$pr)
  colnames(cpp_out) <- NULL
  colnames(xt_out)  <- NULL
  expect_equal(xt_out, cpp_out)
})

test_that("bioclim_xt handles all-uniform monthly values", {
  # All months same value: edge case for SD, argmax, argmin
  tas_u <- rep(15.0, 12)
  pr_u  <- rep(50.0, 12)
  m <- make_mats(tas_u, tas_u + 5, tas_u - 5, pr_u)
  cpp_out <- bioclim_cpp(m$tas, m$tasmax, m$tasmin, m$pr)
  xt_out  <- bioclim_xt(m$tas,  m$tasmax, m$tasmin, m$pr)
  colnames(cpp_out) <- NULL
  colnames(xt_out)  <- NULL
  expect_equal(xt_out, cpp_out)
})

test_that("bioclim_xt returns matrix with correct dimensions", {
  set.seed(42)
  n <- 10L
  tas_r    <- matrix(runif(n * 12, -10, 30), nrow = n)
  tasmax_r <- tas_r + runif(n * 12, 0, 10)
  tasmin_r <- tas_r - runif(n * 12, 0, 10)
  pr_r     <- matrix(runif(n * 12, 0, 200), nrow = n)
  out <- bioclim_xt(tas_r, tasmax_r, tasmin_r, pr_r)
  expect_equal(nrow(out), n)
  expect_equal(ncol(out), 19L)
  expect_equal(colnames(out),
               paste0("bio", formatC(1:19, width = 2, flag = "0")))
})
