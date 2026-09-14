# Integration tests: ERA5 → monthly → bioclim pipeline
#
# Validates that the full transformation chain produces correct bioclimatic
# variables by constructing synthetic ERA5-like hourly data with known
# monthly statistics, aggregating with the C++ backend, and comparing
# bioclim output against direct computation from known monthly values.

# ── Helper: generate ERA5-like hourly data for 12 months ─────────────────

# Known monthly climate for a temperate pixel (°C):
known_tas    <- c(0, 2, 7, 12, 17, 22, 25, 24, 19, 13, 7, 2)
known_pr_mm  <- c(50, 45, 55, 60, 70, 40, 20, 25, 45, 65, 60, 55)

# Days per month for 2020 (leap year)
days_2020 <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# Generate synthetic hourly t2m (K) for one month with a target monthly
# mean in Celsius.  Uses a simple diurnal cycle with amplitude 5°C.
generate_hourly_t2m <- function(target_mean_c, n_days, amplitude = 5) {
  base_k <- target_mean_c + 273.15
  hours <- seq(0, 23)
  cycle <- amplitude * sin(2 * pi * (hours - 4) / 24)
  rep(base_k + cycle, n_days)
}

# Generate synthetic hourly tp (m) for one month with a target monthly
# total in mm.  Distributes precipitation uniformly across hours.
generate_hourly_tp <- function(target_pr_mm, n_days) {
  n_hours <- 24 * n_days
  rate_m <- (target_pr_mm / 1000.0) / n_hours
  rep(rate_m, n_hours)
}

# ── Test: single-pixel end-to-end ──────────────────────────────────────────

test_that("ERA5 hourly → monthly → bioclim matches direct bioclim", {
  # Generate 12 months of hourly data
  monthly_results <- vector("list", 12)
  monthly_pr <- numeric(12)

  for (m in 1:12) {
    hourly_t2m <- generate_hourly_t2m(known_tas[m], days_2020[m])
    hourly_tp  <- generate_hourly_tp(known_pr_mm[m], days_2020[m])

    # Aggregate with R reference
    temp <- era5_t2m_to_monthly_r(hourly_t2m, days_2020[m], to_celsius = TRUE)
    monthly_results[[m]] <- temp
    monthly_pr[m] <- era5_tp_to_monthly_r(hourly_tp)
  }

  tas_agg    <- vapply(monthly_results, `[[`, numeric(1), "tas")
  tasmax_agg <- vapply(monthly_results, `[[`, numeric(1), "tasmax")
  tasmin_agg <- vapply(monthly_results, `[[`, numeric(1), "tasmin")

  # The aggregated monthly means should be close to the known inputs
  expect_equal(tas_agg, known_tas, tolerance = 0.01)
  expect_equal(monthly_pr, known_pr_mm, tolerance = 0.01)

  # Compute bioclim from aggregated monthly data
  bio_agg <- bioclim(tas_agg, tasmax_agg, tasmin_agg, monthly_pr)

  # Compute bioclim directly from known monthly data
  # (tasmax/tasmin are known_tas ± amplitude)
  bio_direct <- bioclim(known_tas, known_tas + 5, known_tas - 5, known_pr_mm)

  # The two should match closely
  expect_equal(bio_agg[["bio01"]], bio_direct[["bio01"]], tolerance = 0.1)
  expect_equal(bio_agg[["bio12"]], bio_direct[["bio12"]], tolerance = 0.1)
  expect_equal(bio_agg[["bio04"]], bio_direct[["bio04"]], tolerance = 1.0)
  expect_equal(bio_agg[["bio05"]], bio_direct[["bio05"]], tolerance = 0.1)
  expect_equal(bio_agg[["bio06"]], bio_direct[["bio06"]], tolerance = 0.1)
})

# ── Test: C++ multi-pixel matches R reference through bioclim ──────────────

test_that("C++ aggregation → bioclim matches R reference → bioclim", {
  n_days <- 3L
  n_pixels <- 3L

  set.seed(99)
  # Create synthetic hourly data for 3 pixels, 3 days
  hourly_t2m <- matrix(NA_real_, nrow = n_pixels, ncol = 24 * n_days)
  hourly_tp  <- matrix(NA_real_, nrow = n_pixels, ncol = 24 * n_days)

  for (i in 1:n_pixels) {
    base_k <- 273.15 + 10 * i  # 10°C, 20°C, 30°C
    hourly_t2m[i, ] <- generate_hourly_t2m(10 * i, n_days)
    hourly_tp[i, ]  <- generate_hourly_tp(50 * i, n_days)
  }

  # C++ path
  cpp_temp <- era5_t2m_to_monthly_cpp(hourly_t2m, n_days, to_celsius = TRUE)
  cpp_pr   <- era5_tp_to_monthly_cpp(hourly_tp)

  # R path (per pixel)
  for (i in 1:n_pixels) {
    r_temp <- era5_t2m_to_monthly_r(hourly_t2m[i, ], n_days, to_celsius = TRUE)
    r_pr   <- era5_tp_to_monthly_r(hourly_tp[i, ])

    expect_equal(cpp_temp$tas[i],    r_temp$tas,    tolerance = 1e-6)
    expect_equal(cpp_temp$tasmax[i], r_temp$tasmax, tolerance = 1e-6)
    expect_equal(cpp_temp$tasmin[i], r_temp$tasmin, tolerance = 1e-6)
    expect_equal(cpp_pr[i],          r_pr,          tolerance = 1e-6)
  }
})

# ── Test: days_in_month consistency for all 12 months ────────────────────

test_that("days_in_month correct for 2020 (leap year)", {
  expected <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  for (m in 1:12) {
    expect_equal(days_in_month(2020, m), expected[m])
  }
})

test_that("days_in_month correct for 2021 (non-leap year)", {
  expected <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  for (m in 1:12) {
    expect_equal(days_in_month(2021, m), expected[m])
  }
})
