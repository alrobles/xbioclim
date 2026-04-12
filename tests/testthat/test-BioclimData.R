# Tests for the BioclimData S4 class and its C++-backed methods.
#
# Correspondence between R interface and xbioclim C++ implementation:
#   bio01_cpp   <->  bio01("BioclimData")  -  Mean Annual Temperature
#   bio02_cpp   <->  bio02("BioclimData")  -  Mean Diurnal Range
#   bio03_cpp   <->  bio03("BioclimData")  -  Isothermality
#   bio04_cpp   <->  bio04("BioclimData")  -  Temperature Seasonality
#   bio05_cpp   <->  bio05("BioclimData")  -  Max Temp Warmest Month
#   bio06_cpp   <->  bio06("BioclimData")  -  Min Temp Coldest Month
#   bio07_cpp   <->  bio07("BioclimData")  -  Temperature Annual Range
#   bio08_cpp   <->  bio08("BioclimData")  -  Mean Temp Wettest Quarter
#   bio09_cpp   <->  bio09("BioclimData")  -  Mean Temp Driest Quarter
#   bio10_cpp   <->  bio10("BioclimData")  -  Mean Temp Warmest Quarter
#   bio11_cpp   <->  bio11("BioclimData")  -  Mean Temp Coldest Quarter
#   bio12_cpp   <->  bio12("BioclimData")  -  Annual Precipitation
#   bio13_cpp   <->  bio13("BioclimData")  -  Precipitation Wettest Month
#   bio14_cpp   <->  bio14("BioclimData")  -  Precipitation Driest Month
#   bio15_cpp   <->  bio15("BioclimData")  -  Precipitation Seasonality
#   bio16_cpp   <->  bio16("BioclimData")  -  Precip Wettest Quarter
#   bio17_cpp   <->  bio17("BioclimData")  -  Precip Driest Quarter
#   bio18_cpp   <->  bio18("BioclimData")  -  Precip Warmest Quarter
#   bio19_cpp   <->  bio19("BioclimData")  -  Precip Coldest Quarter
#   bioclim_cpp <->  bioclim("BioclimData") - All 19 variables

# ── Mock data ─────────────────────────────────────────────────────────────────

# Same convention as test-bioclim.R: tas[m]=m, tasmax[m]=m+1, tasmin[m]=m-1
vec_tas    <- as.numeric(1:12)
vec_tasmax <- as.numeric(2:13)
vec_tasmin <- as.numeric(0:11)
vec_pr     <- as.numeric(1:12)
vec_pr_rev <- as.numeric(12:1)

tol <- 1e-4

# ── BioclimData constructor ───────────────────────────────────────────────────

test_that("BioclimData() accepts single-pixel vectors and promotes to matrix", {
  bd <- BioclimData(vec_tas, vec_tasmax, vec_tasmin, vec_pr)
  expect_s4_class(bd, "BioclimData")
  expect_true(is.matrix(bd@tas))
  expect_equal(nrow(bd@tas), 1L)
  expect_equal(ncol(bd@tas), 12L)
})

test_that("BioclimData() accepts multi-pixel matrices", {
  mat_tas    <- rbind(vec_tas,    vec_tas)
  mat_tasmax <- rbind(vec_tasmax, vec_tasmax)
  mat_tasmin <- rbind(vec_tasmin, vec_tasmin)
  mat_pr     <- rbind(vec_pr,     vec_pr)
  bd <- BioclimData(mat_tas, mat_tasmax, mat_tasmin, mat_pr)
  expect_s4_class(bd, "BioclimData")
  expect_equal(nrow(bd@tas), 2L)
})

test_that("BioclimData() rejects non-numeric vectors", {
  expect_error(BioclimData(letters[1:12], vec_tasmax, vec_tasmin, vec_pr),
               "must be numeric")
})

test_that("BioclimData() rejects vectors of wrong length", {
  expect_error(BioclimData(1:6, vec_tasmax, vec_tasmin, vec_pr),
               "must have length 12")
})

test_that("BioclimData() rejects matrices with wrong number of columns", {
  expect_error(
    BioclimData(matrix(1:13, nrow = 1), vec_tasmax, vec_tasmin, vec_pr),
    "must have 12 columns"
  )
})

test_that("BioclimData() rejects slots with mismatched row counts", {
  mat2 <- matrix(rep(vec_tas, 2), nrow = 2)
  expect_error(
    BioclimData(mat2, vec_tasmax, vec_tasmin, vec_pr),
    "must have 2 rows"
  )
})

# ── show method ───────────────────────────────────────────────────────────────

test_that("show() prints a one-line summary without error", {
  bd <- BioclimData(vec_tas, vec_tasmax, vec_tasmin, vec_pr)
  expect_output(show(bd), "BioclimData")
})

# ── Single-variable S4 methods – output shape and type ────────────────────────

test_that("bio01(BioclimData) returns numeric vector of length n_pixels", {
  bd <- BioclimData(vec_tas, vec_tasmax, vec_tasmin, vec_pr)
  r  <- bio01(bd)
  expect_type(r, "double")
  expect_length(r, 1L)
})

test_that("bio02(BioclimData) returns numeric vector of length n_pixels", {
  bd <- BioclimData(vec_tas, vec_tasmax, vec_tasmin, vec_pr)
  r  <- bio02(bd)
  expect_type(r, "double")
  expect_length(r, 1L)
})

# ── Single-variable S4 methods – values match R reference implementation ──────

bd1 <- BioclimData(vec_tas, vec_tasmax, vec_tasmin, vec_pr)

test_that("bio01(BioclimData) matches bio01(numeric)", {
  expect_equal(bio01(bd1)[[1]], bio01(vec_tas), tolerance = tol)
})

test_that("bio02(BioclimData) matches bio02(numeric)", {
  expect_equal(bio02(bd1)[[1]], bio02(vec_tasmax, vec_tasmin), tolerance = tol)
})

test_that("bio03(BioclimData) matches bio03(numeric)", {
  expect_equal(bio03(bd1)[[1]], bio03(vec_tasmax, vec_tasmin), tolerance = tol)
})

test_that("bio04(BioclimData) matches bio04(numeric)", {
  expect_equal(bio04(bd1)[[1]], bio04(vec_tas), tolerance = tol)
})

test_that("bio05(BioclimData) matches bio05(numeric)", {
  expect_equal(bio05(bd1)[[1]], bio05(vec_tasmax), tolerance = tol)
})

test_that("bio06(BioclimData) matches bio06(numeric)", {
  expect_equal(bio06(bd1)[[1]], bio06(vec_tasmin), tolerance = tol)
})

test_that("bio07(BioclimData) matches bio07(numeric)", {
  expect_equal(bio07(bd1)[[1]], bio07(vec_tasmax, vec_tasmin), tolerance = tol)
})

test_that("bio08(BioclimData) matches bio08(numeric)", {
  expect_equal(bio08(bd1)[[1]], bio08(vec_tas, vec_pr), tolerance = tol)
})

test_that("bio09(BioclimData) matches bio09(numeric)", {
  expect_equal(bio09(bd1)[[1]], bio09(vec_tas, vec_pr), tolerance = tol)
})

test_that("bio10(BioclimData) matches bio10(numeric)", {
  expect_equal(bio10(bd1)[[1]], bio10(vec_tas), tolerance = tol)
})

test_that("bio11(BioclimData) matches bio11(numeric)", {
  expect_equal(bio11(bd1)[[1]], bio11(vec_tas), tolerance = tol)
})

test_that("bio12(BioclimData) matches bio12(numeric)", {
  expect_equal(bio12(bd1)[[1]], bio12(vec_pr), tolerance = tol)
})

test_that("bio13(BioclimData) matches bio13(numeric)", {
  expect_equal(bio13(bd1)[[1]], bio13(vec_pr), tolerance = tol)
})

test_that("bio14(BioclimData) matches bio14(numeric)", {
  expect_equal(bio14(bd1)[[1]], bio14(vec_pr), tolerance = tol)
})

test_that("bio15(BioclimData) matches bio15(numeric)", {
  expect_equal(bio15(bd1)[[1]], bio15(vec_pr), tolerance = tol)
})

test_that("bio16(BioclimData) matches bio16(numeric)", {
  expect_equal(bio16(bd1)[[1]], bio16(vec_pr), tolerance = tol)
})

test_that("bio17(BioclimData) matches bio17(numeric)", {
  expect_equal(bio17(bd1)[[1]], bio17(vec_pr), tolerance = tol)
})

test_that("bio18(BioclimData) matches bio18(numeric)", {
  expect_equal(bio18(bd1)[[1]], bio18(vec_tas, vec_pr), tolerance = tol)
})

test_that("bio19(BioclimData) matches bio19(numeric)", {
  expect_equal(bio19(bd1)[[1]], bio19(vec_tas, vec_pr), tolerance = tol)
})

# ── bio03 NaN edge case via BioclimData ───────────────────────────────────────

test_that("bio03(BioclimData) returns NaN when annual range is zero", {
  bd_eq <- BioclimData(rep(10, 12), rep(10, 12), rep(10, 12), rep(1, 12))
  expect_true(is.nan(bio03(bd_eq)[[1]]))
})

# ── bio15 NaN edge case via BioclimData ───────────────────────────────────────

test_that("bio15(BioclimData) returns NaN when mean precipitation is zero", {
  bd_zero <- BioclimData(vec_tas, vec_tasmax, vec_tasmin, rep(0, 12))
  expect_true(is.nan(bio15(bd_zero)[[1]]))
})

# ── Batch bioclim(BioclimData) ────────────────────────────────────────────────

test_that("bioclim(BioclimData) returns a matrix with 19 named columns", {
  r <- bioclim(bd1)
  expect_true(is.matrix(r))
  expect_equal(ncol(r), 19L)
  expect_equal(nrow(r), 1L)
  expect_equal(colnames(r), paste0("bio", sprintf("%02d", 1:19)))
})

test_that("bioclim(BioclimData) column values match individual S4 methods", {
  r <- bioclim(bd1)
  expect_equal(r[1, "bio01"], bio01(bd1)[[1]], tolerance = tol, ignore_attr = TRUE)
  expect_equal(r[1, "bio02"], bio02(bd1)[[1]], tolerance = tol, ignore_attr = TRUE)
  expect_equal(r[1, "bio03"], bio03(bd1)[[1]], tolerance = tol, ignore_attr = TRUE)
  expect_equal(r[1, "bio04"], bio04(bd1)[[1]], tolerance = tol, ignore_attr = TRUE)
  expect_equal(r[1, "bio05"], bio05(bd1)[[1]], tolerance = tol, ignore_attr = TRUE)
  expect_equal(r[1, "bio06"], bio06(bd1)[[1]], tolerance = tol, ignore_attr = TRUE)
  expect_equal(r[1, "bio07"], bio07(bd1)[[1]], tolerance = tol, ignore_attr = TRUE)
  expect_equal(r[1, "bio08"], bio08(bd1)[[1]], tolerance = tol, ignore_attr = TRUE)
  expect_equal(r[1, "bio09"], bio09(bd1)[[1]], tolerance = tol, ignore_attr = TRUE)
  expect_equal(r[1, "bio10"], bio10(bd1)[[1]], tolerance = tol, ignore_attr = TRUE)
  expect_equal(r[1, "bio11"], bio11(bd1)[[1]], tolerance = tol, ignore_attr = TRUE)
  expect_equal(r[1, "bio12"], bio12(bd1)[[1]], tolerance = tol, ignore_attr = TRUE)
  expect_equal(r[1, "bio13"], bio13(bd1)[[1]], tolerance = tol, ignore_attr = TRUE)
  expect_equal(r[1, "bio14"], bio14(bd1)[[1]], tolerance = tol, ignore_attr = TRUE)
  expect_equal(r[1, "bio15"], bio15(bd1)[[1]], tolerance = tol, ignore_attr = TRUE)
  expect_equal(r[1, "bio16"], bio16(bd1)[[1]], tolerance = tol, ignore_attr = TRUE)
  expect_equal(r[1, "bio17"], bio17(bd1)[[1]], tolerance = tol, ignore_attr = TRUE)
  expect_equal(r[1, "bio18"], bio18(bd1)[[1]], tolerance = tol, ignore_attr = TRUE)
  expect_equal(r[1, "bio19"], bio19(bd1)[[1]], tolerance = tol, ignore_attr = TRUE)
})

test_that("bioclim(BioclimData) matches bioclim(numeric) values", {
  ref <- bioclim(vec_tas, vec_tasmax, vec_tasmin, vec_pr)
  r   <- bioclim(bd1)
  for (nm in names(ref)) {
    expect_equal(r[1, nm], ref[[nm]], tolerance = tol,
                 ignore_attr = TRUE,
                 label = paste0("bioclim() variable ", nm))
  }
})

# ── Multi-pixel (raster block) path ───────────────────────────────────────────

test_that("bio01(BioclimData) with n>1 pixels returns vector of length n", {
  mat_tas    <- rbind(vec_tas,    rev(vec_tas))
  mat_tasmax <- rbind(vec_tasmax, rev(vec_tasmax))
  mat_tasmin <- rbind(vec_tasmin, rev(vec_tasmin))
  mat_pr     <- rbind(vec_pr,     rev(vec_pr))
  bd2 <- BioclimData(mat_tas, mat_tasmax, mat_tasmin, mat_pr)

  r <- bio01(bd2)
  expect_type(r, "double")
  expect_length(r, 2L)
  expect_equal(r[[1]], bio01(vec_tas),       tolerance = tol)
  expect_equal(r[[2]], bio01(rev(vec_tas)),   tolerance = tol)
})

test_that("bioclim(BioclimData) with n>1 pixels returns matrix with n rows", {
  mat_tas    <- rbind(vec_tas,    rev(vec_tas))
  mat_tasmax <- rbind(vec_tasmax, rev(vec_tasmax))
  mat_tasmin <- rbind(vec_tasmin, rev(vec_tasmin))
  mat_pr     <- rbind(vec_pr,     vec_pr_rev)
  bd2 <- BioclimData(mat_tas, mat_tasmax, mat_tasmin, mat_pr)

  r <- bioclim(bd2)
  expect_true(is.matrix(r))
  expect_equal(nrow(r), 2L)
  expect_equal(ncol(r), 19L)

  ref1 <- bioclim(vec_tas, vec_tasmax, vec_tasmin, vec_pr)
  ref2 <- bioclim(rev(vec_tas), rev(vec_tasmax), rev(vec_tasmin), vec_pr_rev)
  for (nm in names(ref1)) {
    expect_equal(r[1, nm], ref1[[nm]], tolerance = tol,
                 ignore_attr = TRUE,
                 label = paste0("pixel 1 bioclim() variable ", nm))
    expect_equal(r[2, nm], ref2[[nm]], tolerance = tol,
                 ignore_attr = TRUE,
                 label = paste0("pixel 2 bioclim() variable ", nm))
  }
})

# ── xbioclim reference values (single pixel) ──────────────────────────────────

test_that("BioclimData S4 methods match xbioclim reference values", {
  r <- bioclim(bd1)
  expect_equal(r[1, "bio01"],  6.5,       tolerance = tol,  ignore_attr = TRUE)
  expect_equal(r[1, "bio02"],  2.0,       tolerance = tol,  ignore_attr = TRUE)
  expect_equal(r[1, "bio03"],  15.3846,   tolerance = 1e-3, ignore_attr = TRUE)
  expect_equal(r[1, "bio04"],  345.2053,  tolerance = 1e-3, ignore_attr = TRUE)
  expect_equal(r[1, "bio05"],  13.0,      tolerance = tol,  ignore_attr = TRUE)
  expect_equal(r[1, "bio06"],  0.0,       tolerance = tol,  ignore_attr = TRUE)
  expect_equal(r[1, "bio07"],  13.0,      tolerance = tol,  ignore_attr = TRUE)
  expect_equal(r[1, "bio08"],  11.0,      tolerance = tol,  ignore_attr = TRUE)
  expect_equal(r[1, "bio09"],  2.0,       tolerance = tol,  ignore_attr = TRUE)
  expect_equal(r[1, "bio10"],  11.0,      tolerance = tol,  ignore_attr = TRUE)
  expect_equal(r[1, "bio11"],  2.0,       tolerance = tol,  ignore_attr = TRUE)
  expect_equal(r[1, "bio12"],  78.0,      tolerance = tol,  ignore_attr = TRUE)
  expect_equal(r[1, "bio13"],  12.0,      tolerance = tol,  ignore_attr = TRUE)
  expect_equal(r[1, "bio14"],  1.0,       tolerance = tol,  ignore_attr = TRUE)
  expect_equal(r[1, "bio15"],  53.1085,   tolerance = 1e-3, ignore_attr = TRUE)
  expect_equal(r[1, "bio16"],  33.0,      tolerance = tol,  ignore_attr = TRUE)
  expect_equal(r[1, "bio17"],  6.0,       tolerance = tol,  ignore_attr = TRUE)
  expect_equal(r[1, "bio18"],  33.0,      tolerance = tol,  ignore_attr = TRUE)
  expect_equal(r[1, "bio19"],  6.0,       tolerance = tol,  ignore_attr = TRUE)
})
