# Test suite for all 19 bioclimatic variables
# Mock data follows xbioclim test convention:
#   tas[m] = m, tasmax[m] = m+1, tasmin[m] = m-1, pr[m] = m

# ── Mock data setup ──────────────────────────────────────────────────────────

# Test case 1: Linear increasing (matching xbioclim make_test_block)
mock_tas    <- 1:12
mock_tasmax <- 2:13
mock_tasmin <- 0:11
mock_pr     <- 1:12

# Test case 2: Reversed precipitation (temperature/precip anti-correlated)
mock_pr_rev <- 12:1

# Test case 3: Constant temperature, variable precipitation
mock_tas_const    <- rep(20, 12)
mock_tasmax_const <- rep(25, 12)
mock_tasmin_const <- rep(15, 12)
mock_pr_var       <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 50, 20)

tol <- 1e-4

# ── BIO01: Mean Annual Temperature ──────────────────────────────────────────

test_that("bio01 computes mean annual temperature", {
  expect_equal(bio01(mock_tas), 6.5, tolerance = tol)
})

test_that("bio01 with constant temperature", {
  expect_equal(bio01(mock_tas_const), 20.0, tolerance = tol)
})

# ── BIO02: Mean Diurnal Range ───────────────────────────────────────────────

test_that("bio02 computes mean diurnal range", {
  expect_equal(bio02(mock_tasmax, mock_tasmin), 2.0, tolerance = tol)
})

test_that("bio02 with constant diurnal range", {
  expect_equal(bio02(mock_tasmax_const, mock_tasmin_const), 10.0, tolerance = tol)
})

# ── BIO03: Isothermality ────────────────────────────────────────────────────

test_that("bio03 computes isothermality", {
  expect_equal(bio03(mock_tasmax, mock_tasmin), 15.3846, tolerance = 1e-3)
})

test_that("bio03 returns NaN when annual range is zero", {
  tasmax_eq <- rep(10, 12)
  tasmin_eq <- rep(10, 12)
  expect_true(is.nan(bio03(tasmax_eq, tasmin_eq)))
})

test_that("bio03 with constant diurnal range", {
  expect_equal(bio03(mock_tasmax_const, mock_tasmin_const), 100.0, tolerance = tol)
})

# ── BIO04: Temperature Seasonality ──────────────────────────────────────────

test_that("bio04 computes temperature seasonality", {
  expect_equal(bio04(mock_tas), 345.2053, tolerance = 1e-3)
})

test_that("bio04 is zero for constant temperature", {
  expect_equal(bio04(mock_tas_const), 0.0, tolerance = tol)
})

# ── BIO05: Max Temperature of Warmest Month ─────────────────────────────────

test_that("bio05 computes max of monthly max temperature", {
  expect_equal(bio05(mock_tasmax), 13.0, tolerance = tol)
})

test_that("bio05 with constant max temperature", {
  expect_equal(bio05(mock_tasmax_const), 25.0, tolerance = tol)
})

# ── BIO06: Min Temperature of Coldest Month ─────────────────────────────────

test_that("bio06 computes min of monthly min temperature", {
  expect_equal(bio06(mock_tasmin), 0.0, tolerance = tol)
})

test_that("bio06 with constant min temperature", {
  expect_equal(bio06(mock_tasmin_const), 15.0, tolerance = tol)
})

# ── BIO07: Temperature Annual Range ─────────────────────────────────────────

test_that("bio07 computes temperature annual range", {
  expect_equal(bio07(mock_tasmax, mock_tasmin), 13.0, tolerance = tol)
})

test_that("bio07 with constant temperatures", {
  expect_equal(bio07(mock_tasmax_const, mock_tasmin_const), 10.0, tolerance = tol)
})

# ── BIO08: Mean Temperature of Wettest Quarter ──────────────────────────────

test_that("bio08 computes mean temp of wettest quarter (correlated)", {
  # When pr = 1:12, wettest quarter starts at month 10 (Oct-Nov-Dec)
  expect_equal(bio08(mock_tas, mock_pr), 11.0, tolerance = tol)
})

test_that("bio08 computes mean temp of wettest quarter (anti-correlated)", {
  # When pr = 12:1, wettest quarter starts at month 1 (Jan-Feb-Mar)
  expect_equal(bio08(mock_tas, mock_pr_rev), 2.0, tolerance = tol)
})

# ── BIO09: Mean Temperature of Driest Quarter ───────────────────────────────

test_that("bio09 computes mean temp of driest quarter (correlated)", {
  # When pr = 1:12, driest quarter starts at month 1 (Jan-Feb-Mar)
  expect_equal(bio09(mock_tas, mock_pr), 2.0, tolerance = tol)
})

test_that("bio09 computes mean temp of driest quarter (anti-correlated)", {
  # When pr = 12:1, driest quarter starts at month 10 (Oct-Nov-Dec)
  expect_equal(bio09(mock_tas, mock_pr_rev), 11.0, tolerance = tol)
})

# ── BIO10: Mean Temperature of Warmest Quarter ──────────────────────────────

test_that("bio10 computes mean temp of warmest quarter", {
  # Warmest quarter starts at month 10 (Oct-Nov-Dec)
  expect_equal(bio10(mock_tas), 11.0, tolerance = tol)
})

test_that("bio10 with constant temperature returns that temperature", {
  expect_equal(bio10(mock_tas_const), 20.0, tolerance = tol)
})

# ── BIO11: Mean Temperature of Coldest Quarter ──────────────────────────────

test_that("bio11 computes mean temp of coldest quarter", {
  # Coldest quarter starts at month 1 (Jan-Feb-Mar)
  expect_equal(bio11(mock_tas), 2.0, tolerance = tol)
})

test_that("bio11 with constant temperature returns that temperature", {
  expect_equal(bio11(mock_tas_const), 20.0, tolerance = tol)
})

# ── BIO12: Annual Precipitation ─────────────────────────────────────────────

test_that("bio12 computes annual precipitation", {
  expect_equal(bio12(mock_pr), 78.0, tolerance = tol)
})

test_that("bio12 with reversed precipitation is the same sum", {
  expect_equal(bio12(mock_pr_rev), 78.0, tolerance = tol)
})

# ── BIO13: Precipitation of Wettest Month ───────────────────────────────────

test_that("bio13 computes max monthly precipitation", {
  expect_equal(bio13(mock_pr), 12.0, tolerance = tol)
})

test_that("bio13 with variable precipitation", {
  expect_equal(bio13(mock_pr_var), 100.0, tolerance = tol)
})

# ── BIO14: Precipitation of Driest Month ────────────────────────────────────

test_that("bio14 computes min monthly precipitation", {
  expect_equal(bio14(mock_pr), 1.0, tolerance = tol)
})

test_that("bio14 with variable precipitation", {
  expect_equal(bio14(mock_pr_var), 10.0, tolerance = tol)
})

# ── BIO15: Precipitation Seasonality (CV) ───────────────────────────────────

test_that("bio15 computes precipitation seasonality", {
  expect_equal(bio15(mock_pr), 53.1085, tolerance = 1e-3)
})

test_that("bio15 returns NaN when precipitation is zero", {
  expect_true(is.nan(bio15(rep(0, 12))))
})

test_that("bio15 is zero for constant precipitation", {
  expect_equal(bio15(rep(50, 12)), 0.0, tolerance = tol)
})

# ── BIO16: Precipitation of Wettest Quarter ─────────────────────────────────

test_that("bio16 computes precipitation of wettest quarter", {
  # Wettest quarter starts at month 10, sum = 10+11+12 = 33
  expect_equal(bio16(mock_pr), 33.0, tolerance = tol)
})

test_that("bio16 with reversed precipitation", {
  # Wettest quarter starts at month 1, sum = 12+11+10 = 33
  expect_equal(bio16(mock_pr_rev), 33.0, tolerance = tol)
})

# ── BIO17: Precipitation of Driest Quarter ──────────────────────────────────

test_that("bio17 computes precipitation of driest quarter", {
  # Driest quarter starts at month 1, sum = 1+2+3 = 6
  expect_equal(bio17(mock_pr), 6.0, tolerance = tol)
})

test_that("bio17 with reversed precipitation", {
  # Driest quarter starts at month 10, sum = 3+2+1 = 6
  expect_equal(bio17(mock_pr_rev), 6.0, tolerance = tol)
})

# ── BIO18: Precipitation of Warmest Quarter ─────────────────────────────────

test_that("bio18 computes precipitation of warmest quarter (correlated)", {
  # Warmest quarter starts at month 10, pr[10]+pr[11]+pr[12] = 10+11+12 = 33
  expect_equal(bio18(mock_tas, mock_pr), 33.0, tolerance = tol)
})

test_that("bio18 computes precipitation of warmest quarter (anti-correlated)", {
  # Warmest quarter starts at month 10, pr_rev[10]+pr_rev[11]+pr_rev[12] = 3+2+1 = 6
  expect_equal(bio18(mock_tas, mock_pr_rev), 6.0, tolerance = tol)
})

# ── BIO19: Precipitation of Coldest Quarter ─────────────────────────────────

test_that("bio19 computes precipitation of coldest quarter (correlated)", {
  # Coldest quarter starts at month 1, pr[1]+pr[2]+pr[3] = 1+2+3 = 6
  expect_equal(bio19(mock_tas, mock_pr), 6.0, tolerance = tol)
})

test_that("bio19 computes precipitation of coldest quarter (anti-correlated)", {
  # Coldest quarter starts at month 1, pr_rev[1]+pr_rev[2]+pr_rev[3] = 12+11+10 = 33
  expect_equal(bio19(mock_tas, mock_pr_rev), 33.0, tolerance = tol)
})

# ── bioclim() unified function ──────────────────────────────────────────────

test_that("bioclim() returns named vector of length 19", {
  result <- bioclim(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_length(result, 19)
  expect_named(result, paste0("bio", sprintf("%02d", 1:19)))
})

test_that("bioclim() matches individual functions", {
  result <- bioclim(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_equal(result[["bio01"]], bio01(mock_tas), tolerance = tol)
  expect_equal(result[["bio02"]], bio02(mock_tasmax, mock_tasmin), tolerance = tol)
  expect_equal(result[["bio03"]], bio03(mock_tasmax, mock_tasmin), tolerance = tol)
  expect_equal(result[["bio04"]], bio04(mock_tas), tolerance = tol)
  expect_equal(result[["bio05"]], bio05(mock_tasmax), tolerance = tol)
  expect_equal(result[["bio06"]], bio06(mock_tasmin), tolerance = tol)
  expect_equal(result[["bio07"]], bio07(mock_tasmax, mock_tasmin), tolerance = tol)
  expect_equal(result[["bio08"]], bio08(mock_tas, mock_pr), tolerance = tol)
  expect_equal(result[["bio09"]], bio09(mock_tas, mock_pr), tolerance = tol)
  expect_equal(result[["bio10"]], bio10(mock_tas), tolerance = tol)
  expect_equal(result[["bio11"]], bio11(mock_tas), tolerance = tol)
  expect_equal(result[["bio12"]], bio12(mock_pr), tolerance = tol)
  expect_equal(result[["bio13"]], bio13(mock_pr), tolerance = tol)
  expect_equal(result[["bio14"]], bio14(mock_pr), tolerance = tol)
  expect_equal(result[["bio15"]], bio15(mock_pr), tolerance = tol)
  expect_equal(result[["bio16"]], bio16(mock_pr), tolerance = tol)
  expect_equal(result[["bio17"]], bio17(mock_pr), tolerance = tol)
  expect_equal(result[["bio18"]], bio18(mock_tas, mock_pr), tolerance = tol)
  expect_equal(result[["bio19"]], bio19(mock_tas, mock_pr), tolerance = tol)
})

test_that("bioclim() exact values match xbioclim reference", {
  result <- bioclim(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_equal(result[["bio01"]],  6.5,       tolerance = tol)
  expect_equal(result[["bio02"]],  2.0,       tolerance = tol)
  expect_equal(result[["bio03"]],  15.3846,   tolerance = 1e-3)
  expect_equal(result[["bio04"]],  345.2053,  tolerance = 1e-3)
  expect_equal(result[["bio05"]],  13.0,      tolerance = tol)
  expect_equal(result[["bio06"]],  0.0,       tolerance = tol)
  expect_equal(result[["bio07"]],  13.0,      tolerance = tol)
  expect_equal(result[["bio08"]],  11.0,      tolerance = tol)
  expect_equal(result[["bio09"]],  2.0,       tolerance = tol)
  expect_equal(result[["bio10"]],  11.0,      tolerance = tol)
  expect_equal(result[["bio11"]],  2.0,       tolerance = tol)
  expect_equal(result[["bio12"]],  78.0,      tolerance = tol)
  expect_equal(result[["bio13"]],  12.0,      tolerance = tol)
  expect_equal(result[["bio14"]],  1.0,       tolerance = tol)
  expect_equal(result[["bio15"]],  53.1085,   tolerance = 1e-3)
  expect_equal(result[["bio16"]],  33.0,      tolerance = tol)
  expect_equal(result[["bio17"]],  6.0,       tolerance = tol)
  expect_equal(result[["bio18"]],  33.0,      tolerance = tol)
  expect_equal(result[["bio19"]],  6.0,       tolerance = tol)
})

# ── Input validation ────────────────────────────────────────────────────────

test_that("functions reject non-numeric input", {
  expect_error(bio01(letters[1:12]), "must be numeric")
  expect_error(bio02(letters[1:12], mock_tasmin), "must be numeric")
  expect_error(bio12(as.character(1:12)), "must be numeric")
})

test_that("functions reject wrong length input", {
  expect_error(bio01(1:6), "must have length 12")
  expect_error(bio02(1:12, 1:6), "must have length 12")
  expect_error(bio12(1:24), "must have length 12")
})

test_that("bioclim() validates all inputs", {
  expect_error(bioclim(1:6, mock_tasmax, mock_tasmin, mock_pr),
               "must have length 12")
  expect_error(bioclim(mock_tas, 1:6, mock_tasmin, mock_pr),
               "must have length 12")
  expect_error(bioclim(mock_tas, mock_tasmax, 1:6, mock_pr),
               "must have length 12")
  expect_error(bioclim(mock_tas, mock_tasmax, mock_tasmin, 1:6),
               "must have length 12")
})
