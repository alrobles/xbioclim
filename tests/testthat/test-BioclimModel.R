# Tests for the BioclimModel S4 class with C++ pointer handle

# ── Mock data (same as test-bioclim.R) ──────────────────────────────────────

mock_tas    <- 1:12
mock_tasmax <- 2:13
mock_tasmin <- 0:11
mock_pr     <- 1:12
mock_pr_rev <- 12:1

tol <- 1e-4

# ── Construction ─────────────────────────────────────────────────────────────

test_that("BioclimModel() creates a valid object", {
  m <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_s4_class(m, "BioclimModel")
  expect_false(isVirtualClass("BioclimModel"))
})

test_that("BioclimModel pntr slot is an externalptr", {
  m <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_true(is(m@pntr, "externalptr"))
})

test_that("BioclimModel pntr is non-null", {
  m <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_false(bioclim_model_is_null(m@pntr))
})

test_that("BioclimModel() validates tas length", {
  expect_error(BioclimModel(1:6, mock_tasmax, mock_tasmin, mock_pr),
               "must have length 12")
})

test_that("BioclimModel() validates tasmax length", {
  expect_error(BioclimModel(mock_tas, 1:6, mock_tasmin, mock_pr),
               "must have length 12")
})

test_that("BioclimModel() validates tasmin length", {
  expect_error(BioclimModel(mock_tas, mock_tasmax, 1:6, mock_pr),
               "must have length 12")
})

test_that("BioclimModel() validates pr length", {
  expect_error(BioclimModel(mock_tas, mock_tasmax, mock_tasmin, 1:6),
               "must have length 12")
})

test_that("BioclimModel() validates non-numeric tas", {
  expect_error(BioclimModel(letters[1:12], mock_tasmax, mock_tasmin, mock_pr),
               "must be numeric")
})

# ── show method ───────────────────────────────────────────────────────────────

test_that("show() outputs class and pointer info", {
  m <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  out <- capture.output(show(m))
  expect_match(out[1], "BioclimModel")
  expect_match(out[2], "C\\+\\+ BioclimModel")
})

# ── S4 methods: single-argument generics ─────────────────────────────────────

test_that("bio01() dispatches correctly on BioclimModel", {
  m <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_equal(bio01(m), bio01(mock_tas), tolerance = tol)
  expect_equal(bio01(m), 6.5,            tolerance = tol)
})

test_that("bio04() dispatches correctly on BioclimModel", {
  m <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_equal(bio04(m), bio04(mock_tas), tolerance = tol)
  expect_equal(bio04(m), 345.2053,        tolerance = 1e-3)
})

test_that("bio05() dispatches correctly on BioclimModel", {
  m <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_equal(bio05(m), bio05(mock_tasmax), tolerance = tol)
  expect_equal(bio05(m), 13.0,              tolerance = tol)
})

test_that("bio06() dispatches correctly on BioclimModel", {
  m <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_equal(bio06(m), bio06(mock_tasmin), tolerance = tol)
  expect_equal(bio06(m), 0.0,               tolerance = tol)
})

test_that("bio10() dispatches correctly on BioclimModel", {
  m <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_equal(bio10(m), bio10(mock_tas), tolerance = tol)
  expect_equal(bio10(m), 11.0,            tolerance = tol)
})

test_that("bio11() dispatches correctly on BioclimModel", {
  m <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_equal(bio11(m), bio11(mock_tas), tolerance = tol)
  expect_equal(bio11(m), 2.0,             tolerance = tol)
})

test_that("bio12() dispatches correctly on BioclimModel", {
  m <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_equal(bio12(m), bio12(mock_pr), tolerance = tol)
  expect_equal(bio12(m), 78.0,           tolerance = tol)
})

test_that("bio13() dispatches correctly on BioclimModel", {
  m <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_equal(bio13(m), bio13(mock_pr), tolerance = tol)
  expect_equal(bio13(m), 12.0,           tolerance = tol)
})

test_that("bio14() dispatches correctly on BioclimModel", {
  m <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_equal(bio14(m), bio14(mock_pr), tolerance = tol)
  expect_equal(bio14(m), 1.0,            tolerance = tol)
})

test_that("bio15() dispatches correctly on BioclimModel", {
  m <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_equal(bio15(m), bio15(mock_pr), tolerance = 1e-3)
  expect_equal(bio15(m), 53.1085,        tolerance = 1e-3)
})

test_that("bio16() dispatches correctly on BioclimModel", {
  m <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_equal(bio16(m), bio16(mock_pr), tolerance = tol)
  expect_equal(bio16(m), 33.0,           tolerance = tol)
})

test_that("bio17() dispatches correctly on BioclimModel", {
  m <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_equal(bio17(m), bio17(mock_pr), tolerance = tol)
  expect_equal(bio17(m), 6.0,            tolerance = tol)
})

# ── S4 methods: two-argument generics (tasmax, tasmin) ───────────────────────

test_that("bio02() dispatches correctly on BioclimModel", {
  m <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_equal(bio02(m), bio02(mock_tasmax, mock_tasmin), tolerance = tol)
  expect_equal(bio02(m), 2.0, tolerance = tol)
})

test_that("bio03() dispatches correctly on BioclimModel", {
  m <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_equal(bio03(m), bio03(mock_tasmax, mock_tasmin), tolerance = 1e-3)
  expect_equal(bio03(m), 15.3846, tolerance = 1e-3)
})

test_that("bio07() dispatches correctly on BioclimModel", {
  m <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_equal(bio07(m), bio07(mock_tasmax, mock_tasmin), tolerance = tol)
  expect_equal(bio07(m), 13.0, tolerance = tol)
})

# ── S4 methods: two-argument generics (tas, pr) ──────────────────────────────

test_that("bio08() dispatches correctly on BioclimModel", {
  m <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_equal(bio08(m), bio08(mock_tas, mock_pr), tolerance = tol)
  expect_equal(bio08(m), 11.0, tolerance = tol)
})

test_that("bio09() dispatches correctly on BioclimModel", {
  m <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_equal(bio09(m), bio09(mock_tas, mock_pr), tolerance = tol)
  expect_equal(bio09(m), 2.0, tolerance = tol)
})

test_that("bio18() dispatches correctly on BioclimModel", {
  m <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_equal(bio18(m), bio18(mock_tas, mock_pr), tolerance = tol)
  expect_equal(bio18(m), 33.0, tolerance = tol)
})

test_that("bio19() dispatches correctly on BioclimModel", {
  m <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_equal(bio19(m), bio19(mock_tas, mock_pr), tolerance = tol)
  expect_equal(bio19(m), 6.0, tolerance = tol)
})

# ── bioclim() generic on BioclimModel ────────────────────────────────────────

test_that("bioclim() dispatches on BioclimModel and returns named vector", {
  m      <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  result <- bioclim(m)
  expect_length(result, 19)
  expect_named(result, paste0("bio", sprintf("%02d", 1:19)))
})

test_that("bioclim() on BioclimModel matches individual S4 methods", {
  m      <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  result <- bioclim(m)
  expect_equal(result[["bio01"]], bio01(m), tolerance = tol)
  expect_equal(result[["bio02"]], bio02(m), tolerance = tol)
  expect_equal(result[["bio03"]], bio03(m), tolerance = 1e-3)
  expect_equal(result[["bio04"]], bio04(m), tolerance = 1e-3)
  expect_equal(result[["bio05"]], bio05(m), tolerance = tol)
  expect_equal(result[["bio06"]], bio06(m), tolerance = tol)
  expect_equal(result[["bio07"]], bio07(m), tolerance = tol)
  expect_equal(result[["bio08"]], bio08(m), tolerance = tol)
  expect_equal(result[["bio09"]], bio09(m), tolerance = tol)
  expect_equal(result[["bio10"]], bio10(m), tolerance = tol)
  expect_equal(result[["bio11"]], bio11(m), tolerance = tol)
  expect_equal(result[["bio12"]], bio12(m), tolerance = tol)
  expect_equal(result[["bio13"]], bio13(m), tolerance = tol)
  expect_equal(result[["bio14"]], bio14(m), tolerance = tol)
  expect_equal(result[["bio15"]], bio15(m), tolerance = 1e-3)
  expect_equal(result[["bio16"]], bio16(m), tolerance = tol)
  expect_equal(result[["bio17"]], bio17(m), tolerance = tol)
  expect_equal(result[["bio18"]], bio18(m), tolerance = tol)
  expect_equal(result[["bio19"]], bio19(m), tolerance = tol)
})

test_that("bioclim() on BioclimModel matches plain bioclim() values", {
  m         <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  from_obj  <- bioclim(m)
  from_vecs <- bioclim(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_equal(from_obj, from_vecs, tolerance = tol)
})

test_that("bioclim() on BioclimModel matches xbioclim reference values", {
  m      <- BioclimModel(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  result <- bioclim(m)
  expect_equal(result[["bio01"]],  6.5,      tolerance = tol)
  expect_equal(result[["bio02"]],  2.0,      tolerance = tol)
  expect_equal(result[["bio03"]],  15.3846,  tolerance = 1e-3)
  expect_equal(result[["bio04"]],  345.2053, tolerance = 1e-3)
  expect_equal(result[["bio05"]],  13.0,     tolerance = tol)
  expect_equal(result[["bio06"]],  0.0,      tolerance = tol)
  expect_equal(result[["bio07"]],  13.0,     tolerance = tol)
  expect_equal(result[["bio08"]],  11.0,     tolerance = tol)
  expect_equal(result[["bio09"]],  2.0,      tolerance = tol)
  expect_equal(result[["bio10"]],  11.0,     tolerance = tol)
  expect_equal(result[["bio11"]],  2.0,      tolerance = tol)
  expect_equal(result[["bio12"]],  78.0,     tolerance = tol)
  expect_equal(result[["bio13"]],  12.0,     tolerance = tol)
  expect_equal(result[["bio14"]],  1.0,      tolerance = tol)
  expect_equal(result[["bio15"]],  53.1085,  tolerance = 1e-3)
  expect_equal(result[["bio16"]],  33.0,     tolerance = tol)
  expect_equal(result[["bio17"]],  6.0,      tolerance = tol)
  expect_equal(result[["bio18"]],  33.0,     tolerance = tol)
  expect_equal(result[["bio19"]],  6.0,      tolerance = tol)
})

# ── Backward compatibility: numeric vector calls still work ──────────────────

test_that("bio01() still works on numeric vectors after generic conversion", {
  expect_equal(bio01(mock_tas), 6.5, tolerance = tol)
})

test_that("bio02() still works on numeric vectors after generic conversion", {
  expect_equal(bio02(mock_tasmax, mock_tasmin), 2.0, tolerance = tol)
})

test_that("bio08() still works on numeric vectors after generic conversion", {
  expect_equal(bio08(mock_tas, mock_pr), 11.0, tolerance = tol)
})

test_that("bioclim() still works on numeric vectors after generic conversion", {
  result <- bioclim(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_length(result, 19)
  expect_equal(result[["bio01"]], 6.5, tolerance = tol)
})

# ── Edge cases ───────────────────────────────────────────────────────────────

test_that("bio03() returns NaN via BioclimModel when annual range is zero", {
  tasmax_eq <- rep(10, 12)
  tasmin_eq <- rep(10, 12)
  m <- BioclimModel(rep(10, 12), tasmax_eq, tasmin_eq, 1:12)
  expect_true(is.nan(bio03(m)))
})

test_that("bio15() returns NaN via BioclimModel when precipitation is zero", {
  m <- BioclimModel(1:12, 2:13, 0:11, rep(0, 12))
  expect_true(is.nan(bio15(m)))
})
