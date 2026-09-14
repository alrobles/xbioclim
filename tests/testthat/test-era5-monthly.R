# Test suite for ERA5-Land hourly-to-monthly aggregation functions.
#
# Mock data convention:
#   - Hourly temperature: sinusoidal diurnal cycle around a base temperature
#   - Hourly precipitation: small random accumulations in metres
#   - Tests validate against hand-computed expected values

# ── Mock data helpers ────────────────────────────────────────────────────────

# Generate a synthetic diurnal temperature cycle (K) for one day.
# Peak at hour 14 (2 PM), minimum at hour 4 (4 AM).
make_diurnal_cycle <- function(base_k = 285, amplitude = 5, n_hours = 24) {
  hours <- seq(0, n_hours - 1)
  base_k + amplitude * sin(2 * pi * (hours - 4) / 24)
}

# Generate hourly temperature for n_days with a known pattern.
make_hourly_t2m <- function(n_days, base_k = 285, amplitude = 5) {
  unlist(lapply(seq_len(n_days), function(d) {
    make_diurnal_cycle(base_k + (d - 1) * 0.1, amplitude)
  }))
}

# Generate hourly precipitation for n_days (m).
# Small random positive values with some zeros.
make_hourly_tp <- function(n_days, rate_m = 0.0001) {
  set.seed(123)
  pmax(0, rnorm(n_days * 24, rate_m, rate_m / 2))
}

# ── Constants ────────────────────────────────────────────────────────────────

tol <- 1e-6

# ── Test: R reference implementation ─────────────────────────────────────────

test_that("era5_t2m_to_monthly_r produces correct statistics for constant temp", {
  n_days <- 3L
  # Constant temperature: 285 K every hour
  hourly <- rep(285.0, 24 * n_days)
  result <- era5_t2m_to_monthly_r(hourly, n_days)

  expect_equal(result$tas,    285.0, tolerance = tol)
  expect_equal(result$tasmax, 285.0, tolerance = tol)
  expect_equal(result$tasmin, 285.0, tolerance = tol)
})

test_that("era5_t2m_to_monthly_r handles diurnal cycle correctly", {
  n_days <- 1L
  # Known diurnal cycle: base 285 K, amplitude 5 K
  hourly <- make_diurnal_cycle(285, 5, 24)

  result <- era5_t2m_to_monthly_r(hourly, n_days)

  # Daily mean ≈ base temperature (mean of sine = 0)
  expect_equal(result$tas, mean(hourly), tolerance = tol)
  # Daily max = base + amplitude (peak of sine)
  expect_equal(result$tasmax, max(hourly), tolerance = tol)
  # Daily min = base - amplitude
  expect_equal(result$tasmin, min(hourly), tolerance = tol)
})

test_that("era5_t2m_to_monthly_r Kelvin-to-Celsius conversion", {
  n_days <- 1L
  hourly <- rep(273.15, 24)  # 0°C in Kelvin
  result <- era5_t2m_to_monthly_r(hourly, n_days, to_celsius = TRUE)

  expect_equal(result$tas,    0.0, tolerance = tol)
  expect_equal(result$tasmax, 0.0, tolerance = tol)
  expect_equal(result$tasmin, 0.0, tolerance = tol)
})

test_that("era5_t2m_to_monthly_r multi-day aggregation", {
  n_days <- 3L
  # Day 1: constant 280 K, Day 2: constant 285 K, Day 3: constant 290 K
  hourly <- c(rep(280, 24), rep(285, 24), rep(290, 24))
  result <- era5_t2m_to_monthly_r(hourly, n_days)

  # Monthly mean of daily means = (280 + 285 + 290) / 3
  expect_equal(result$tas, 285.0, tolerance = tol)
  # Monthly mean of daily maxima = same (constant per day)
  expect_equal(result$tasmax, 285.0, tolerance = tol)
  expect_equal(result$tasmin, 285.0, tolerance = tol)
})

test_that("era5_t2m_to_monthly_r with varying daily extremes", {
  n_days <- 2L
  # Day 1: linear 270..293 (24 values)
  day1 <- seq(270, 293, length.out = 24)
  # Day 2: linear 275..298
  day2 <- seq(275, 298, length.out = 24)
  hourly <- c(day1, day2)
  result <- era5_t2m_to_monthly_r(hourly, n_days)

  exp_tas    <- mean(c(mean(day1), mean(day2)))
  exp_tasmax <- mean(c(max(day1), max(day2)))
  exp_tasmin <- mean(c(min(day1), min(day2)))

  expect_equal(result$tas,    exp_tas,    tolerance = tol)
  expect_equal(result$tasmax, exp_tasmax, tolerance = tol)
  expect_equal(result$tasmin, exp_tasmin, tolerance = tol)
})

test_that("era5_tp_to_monthly_r sums and converts correctly", {
  # 3 days of constant 0.001 m/hour precipitation
  n_hours <- 72
  hourly_tp <- rep(0.001, n_hours)
  pr <- era5_tp_to_monthly_r(hourly_tp)
  # Expected: 72 * 0.001 * 1000 = 72 mm
  expect_equal(pr, 72.0, tolerance = tol)
})

test_that("era5_tp_to_monthly_r handles zeros and negatives", {
  hourly_tp <- c(0.001, 0, -0.0001, 0.002, 0, 0)
  pr <- era5_tp_to_monthly_r(hourly_tp)
  # Only positive values: 0.001 + 0.002 = 0.003 m → 3 mm
  expect_equal(pr, 3.0, tolerance = tol)
})

test_that("era5_tp_to_monthly_r handles all-zero precipitation", {
  pr <- era5_tp_to_monthly_r(rep(0, 48))
  expect_equal(pr, 0.0, tolerance = tol)
})

test_that("era5_tp_to_monthly_r handles NAs", {
  hourly_tp <- c(0.001, NA, 0.002, NA, 0.003)
  pr <- era5_tp_to_monthly_r(hourly_tp)
  # 0.001 + 0.002 + 0.003 = 0.006 m → 6 mm
  expect_equal(pr, 6.0, tolerance = tol)
})

# ── Test: exported wrapper dispatching ─────────────────────────────────────

test_that("era5_t2m_to_monthly dispatches vector to R impl", {
  hourly <- rep(285, 48)
  result <- era5_t2m_to_monthly(hourly, n_days = 2L)
  expect_true(is.list(result))
  expect_named(result, c("tas", "tasmax", "tasmin"))
  expect_equal(result$tas, 285.0, tolerance = tol)
})

test_that("era5_tp_to_monthly dispatches vector to R impl", {
  hourly_tp <- rep(0.001, 48)
  pr <- era5_tp_to_monthly(hourly_tp)
  expect_equal(pr, 48.0, tolerance = tol)
})

# ── Test: C++ implementation via matrix input ──────────────────────────────

test_that("era5_t2m_to_monthly_cpp matches R reference for single pixel", {
  n_days <- 3L
  hourly_vec <- make_hourly_t2m(n_days, 285, 5)
  hourly_mat <- matrix(hourly_vec, nrow = 1L)

  r_result   <- era5_t2m_to_monthly_r(hourly_vec, n_days)
  cpp_result <- era5_t2m_to_monthly_cpp(hourly_mat, n_days)

  expect_equal(cpp_result$tas[1],    r_result$tas,    tolerance = tol)
  expect_equal(cpp_result$tasmax[1], r_result$tasmax, tolerance = tol)
  expect_equal(cpp_result$tasmin[1], r_result$tasmin, tolerance = tol)
})

test_that("era5_t2m_to_monthly_cpp matches R reference with Celsius conversion", {
  n_days <- 2L
  hourly_vec <- rep(273.15, 48)
  hourly_mat <- matrix(hourly_vec, nrow = 1L)

  r_result   <- era5_t2m_to_monthly_r(hourly_vec, n_days, to_celsius = TRUE)
  cpp_result <- era5_t2m_to_monthly_cpp(hourly_mat, n_days, to_celsius = TRUE)

  expect_equal(cpp_result$tas[1],    r_result$tas,    tolerance = tol)
  expect_equal(cpp_result$tasmax[1], r_result$tasmax, tolerance = tol)
  expect_equal(cpp_result$tasmin[1], r_result$tasmin, tolerance = tol)
})

test_that("era5_tp_to_monthly_cpp matches R reference for single pixel", {
  hourly_vec <- c(0.001, 0, 0.002, rep(0, 45))
  hourly_mat <- matrix(hourly_vec, nrow = 1L)

  r_result   <- era5_tp_to_monthly_r(hourly_vec)
  cpp_result <- era5_tp_to_monthly_cpp(hourly_mat)

  expect_equal(cpp_result[1], r_result, tolerance = tol)
})

test_that("era5_t2m_to_monthly_cpp handles multiple pixels", {
  n_days <- 2L
  n_pixels <- 5L
  n_hours <- 48L

  set.seed(42)
  hourly_mat <- matrix(
    280 + rnorm(n_pixels * n_hours, 0, 3),
    nrow = n_pixels, ncol = n_hours
  )

  cpp_result <- era5_t2m_to_monthly_cpp(hourly_mat, n_days)

  expect_length(cpp_result$tas,    n_pixels)
  expect_length(cpp_result$tasmax, n_pixels)
  expect_length(cpp_result$tasmin, n_pixels)

  # Cross-validate each pixel against R reference
  for (i in seq_len(n_pixels)) {
    r_result <- era5_t2m_to_monthly_r(hourly_mat[i, ], n_days)
    expect_equal(cpp_result$tas[i],    r_result$tas,    tolerance = tol)
    expect_equal(cpp_result$tasmax[i], r_result$tasmax, tolerance = tol)
    expect_equal(cpp_result$tasmin[i], r_result$tasmin, tolerance = tol)
  }
})

test_that("era5_tp_to_monthly_cpp handles multiple pixels", {
  n_pixels <- 5L
  n_hours <- 48L

  set.seed(42)
  hourly_mat <- matrix(
    pmax(0, rnorm(n_pixels * n_hours, 0.0001, 0.00005)),
    nrow = n_pixels, ncol = n_hours
  )

  cpp_result <- era5_tp_to_monthly_cpp(hourly_mat)
  expect_length(cpp_result, n_pixels)

  for (i in seq_len(n_pixels)) {
    r_result <- era5_tp_to_monthly_r(hourly_mat[i, ])
    expect_equal(cpp_result[i], r_result, tolerance = tol)
  }
})

test_that("era5_t2m_to_monthly_cpp rejects mismatched n_hours", {
  hourly_mat <- matrix(285, nrow = 1, ncol = 50)
  expect_error(era5_t2m_to_monthly_cpp(hourly_mat, n_days = 2L),
               "n_hours")
})

# ── Test: days_in_month helper ──────────────────────────────────────────────

test_that("days_in_month returns correct values", {
  # January 2020

  expect_equal(days_in_month(2020, 1), 31L)
  # February 2020 (leap year)
  expect_equal(days_in_month(2020, 2), 29L)
  # February 2021 (non-leap)
  expect_equal(days_in_month(2021, 2), 28L)
  # April 2020
  expect_equal(days_in_month(2020, 4), 30L)
})

# ── Test: full ERA5 → bioclim pipeline with mock data ──────────────────────

test_that("ERA5 monthly data feeds into bioclim correctly", {
  # Create 12 months of mock monthly data (already aggregated)
  # Simulate a temperate climate with seasonal variation
  tas    <- c(0, 2, 7, 12, 17, 22, 25, 24, 19, 13, 7, 2)
  tasmax <- tas + 5
  tasmin <- tas - 5
  pr     <- c(50, 45, 55, 60, 70, 40, 20, 25, 45, 65, 60, 55)

  # This should work with the existing bioclim function
  bio <- bioclim(tas, tasmax, tasmin, pr)

  expect_length(bio, 19L)
  expect_named(bio, paste0("bio", formatC(1:19, width = 2, flag = "0")))
  # BIO01 = mean annual temperature
  expect_equal(bio[["bio01"]], mean(tas), tolerance = 1e-4)
  # BIO12 = annual precipitation
  expect_equal(bio[["bio12"]], sum(pr), tolerance = 1e-4)
})
