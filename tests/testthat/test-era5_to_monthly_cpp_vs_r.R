# Cross-language consistency tests: era5_to_monthly_cpp() vs era5_to_monthly_r()
#
# Stage 3 QA gate (SSDLC): C++ implementation must match pure-R reference
# to within tolerance 1e-10 on identical inputs.

tol <- 1e-10

# ── Mock data helpers ────────────────────────────────────────────────────────

make_diurnal <- function(base_k = 285, amplitude = 5, n_hours = 24) {
  hours <- seq(0, n_hours - 1)
  base_k + amplitude * sin(2 * pi * (hours - 4) / 24)
}

make_hourly_t2m <- function(n_days, base_k = 285, amplitude = 5) {
  unlist(lapply(seq_len(n_days), function(d) {
    make_diurnal(base_k + (d - 1) * 0.1, amplitude)
  }))
}

make_hourly_tp <- function(n_days, rate_m = 0.0001) {
  set.seed(42)
  pmax(0, rnorm(n_days * 24, rate_m, rate_m / 2))
}

# ── Test: single pixel, standard pattern ─────────────────────────────────

test_that("era5_to_monthly_cpp matches _r on standard diurnal pattern", {
  n_days <- 31L
  t2m_vec <- make_hourly_t2m(n_days, 285, 5)
  tp_vec  <- make_hourly_tp(n_days, 0.0001)

  r_result   <- era5_to_monthly_r(t2m_vec, tp_vec, n_days)
  t2m_mat    <- matrix(t2m_vec, nrow = 1L)
  tp_mat     <- matrix(tp_vec,  nrow = 1L)
  cpp_result <- era5_to_monthly_cpp(t2m_mat, tp_mat, n_days)

  expect_equal(cpp_result$tas[1],    r_result$tas,    tolerance = tol)
  expect_equal(cpp_result$tasmax[1], r_result$tasmax, tolerance = tol)
  expect_equal(cpp_result$tasmin[1], r_result$tasmin, tolerance = tol)
  expect_equal(cpp_result$pr[1],     r_result$pr,     tolerance = tol)
})

# ── Test: single pixel, Celsius conversion ───────────────────────────────

test_that("era5_to_monthly_cpp matches _r with Celsius conversion", {
  n_days <- 28L
  t2m_vec <- make_hourly_t2m(n_days, 273.15, 3)
  tp_vec  <- make_hourly_tp(n_days, 0.0002)

  r_result   <- era5_to_monthly_r(t2m_vec, tp_vec, n_days, to_celsius = TRUE)
  t2m_mat    <- matrix(t2m_vec, nrow = 1L)
  tp_mat     <- matrix(tp_vec,  nrow = 1L)
  cpp_result <- era5_to_monthly_cpp(t2m_mat, tp_mat, n_days, to_celsius = TRUE)

  expect_equal(cpp_result$tas[1],    r_result$tas,    tolerance = tol)
  expect_equal(cpp_result$tasmax[1], r_result$tasmax, tolerance = tol)
  expect_equal(cpp_result$tasmin[1], r_result$tasmin, tolerance = tol)
  expect_equal(cpp_result$pr[1],     r_result$pr,     tolerance = tol)
})

# ── Test: multiple pixels ────────────────────────────────────────────────

test_that("era5_to_monthly_cpp matches _r for 10 pixels", {
  n_days   <- 30L
  n_pixels <- 10L
  n_hours  <- 24L * n_days

  set.seed(99)
  t2m_mat <- matrix(NA_real_, nrow = n_pixels, ncol = n_hours)
  tp_mat  <- matrix(NA_real_, nrow = n_pixels, ncol = n_hours)

  for (i in seq_len(n_pixels)) {
    t2m_mat[i, ] <- make_hourly_t2m(n_days, 270 + i * 3, 4 + i * 0.5)
    set.seed(i * 7)
    tp_mat[i, ]  <- pmax(0, rnorm(n_hours, 0.0001 * i, 0.00005))
  }

  cpp_result <- era5_to_monthly_cpp(t2m_mat, tp_mat, n_days)

  for (i in seq_len(n_pixels)) {
    r_result <- era5_to_monthly_r(t2m_mat[i, ], tp_mat[i, ], n_days)
    expect_equal(cpp_result$tas[i],    r_result$tas,    tolerance = tol,
                 label = sprintf("pixel %d tas", i))
    expect_equal(cpp_result$tasmax[i], r_result$tasmax, tolerance = tol,
                 label = sprintf("pixel %d tasmax", i))
    expect_equal(cpp_result$tasmin[i], r_result$tasmin, tolerance = tol,
                 label = sprintf("pixel %d tasmin", i))
    expect_equal(cpp_result$pr[i],     r_result$pr,     tolerance = tol,
                 label = sprintf("pixel %d pr", i))
  }
})

# ── Test: constant temperature ───────────────────────────────────────────

test_that("era5_to_monthly_cpp matches _r with constant temperature", {
  n_days  <- 5L
  t2m_vec <- rep(290.0, 24 * n_days)
  tp_vec  <- rep(0.001, 24 * n_days)

  r_result   <- era5_to_monthly_r(t2m_vec, tp_vec, n_days)
  cpp_result <- era5_to_monthly_cpp(
    matrix(t2m_vec, nrow = 1L),
    matrix(tp_vec,  nrow = 1L),
    n_days
  )

  expect_equal(cpp_result$tas[1],    r_result$tas,    tolerance = tol)
  expect_equal(cpp_result$tasmax[1], r_result$tasmax, tolerance = tol)
  expect_equal(cpp_result$tasmin[1], r_result$tasmin, tolerance = tol)
  expect_equal(cpp_result$pr[1],     r_result$pr,     tolerance = tol)
})

# ── Test: zero precipitation ─────────────────────────────────────────────

test_that("era5_to_monthly_cpp matches _r with zero precipitation", {
  n_days  <- 3L
  t2m_vec <- make_hourly_t2m(n_days, 300, 8)
  tp_vec  <- rep(0, 24 * n_days)

  r_result   <- era5_to_monthly_r(t2m_vec, tp_vec, n_days)
  cpp_result <- era5_to_monthly_cpp(
    matrix(t2m_vec, nrow = 1L),
    matrix(tp_vec,  nrow = 1L),
    n_days
  )

  expect_equal(cpp_result$pr[1], r_result$pr, tolerance = tol)
  expect_equal(cpp_result$pr[1], 0.0, tolerance = tol)
})

# ── Test: exported wrapper dispatches correctly ──────────────────────────

test_that("era5_to_monthly wrapper returns all 4 variables", {
  n_days  <- 2L
  t2m_vec <- make_hourly_t2m(n_days)
  tp_vec  <- make_hourly_tp(n_days)

  result <- era5_to_monthly(t2m_vec, tp_vec, n_days)
  expect_named(result, c("tas", "tasmax", "tasmin", "pr"))
  expect_length(result$tas, 1L)
  expect_length(result$pr, 1L)
})

test_that("era5_to_monthly matrix path matches vector path", {
  n_days  <- 3L
  t2m_vec <- make_hourly_t2m(n_days, 288, 6)
  tp_vec  <- make_hourly_tp(n_days, 0.0003)

  vec_result <- era5_to_monthly(t2m_vec, tp_vec, n_days)
  mat_result <- era5_to_monthly(
    matrix(t2m_vec, nrow = 1L),
    matrix(tp_vec,  nrow = 1L),
    n_days
  )

  expect_equal(mat_result$tas[1],    vec_result$tas,    tolerance = tol)
  expect_equal(mat_result$tasmax[1], vec_result$tasmax, tolerance = tol)
  expect_equal(mat_result$tasmin[1], vec_result$tasmin, tolerance = tol)
  expect_equal(mat_result$pr[1],     vec_result$pr,     tolerance = tol)
})

# ── Test: validation errors ──────────────────────────────────────────────

test_that("era5_to_monthly_cpp rejects mismatched pixel count", {
  t2m_mat <- matrix(285, nrow = 3, ncol = 48)
  tp_mat  <- matrix(0.001, nrow = 5, ncol = 48)
  expect_error(era5_to_monthly_cpp(t2m_mat, tp_mat, 2L), "same number of pixels")
})

test_that("era5_to_monthly_cpp rejects mismatched n_hours vs n_days", {
  t2m_mat <- matrix(285, nrow = 1, ncol = 50)
  tp_mat  <- matrix(0.001, nrow = 1, ncol = 50)
  expect_error(era5_to_monthly_cpp(t2m_mat, tp_mat, 2L), "columns")
})

# ── Test: full 12-month cross-validation ─────────────────────────────────

test_that("era5_to_monthly_cpp matches _r for all 12 months of 2020", {
  days_2020 <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  base_temps <- c(265, 268, 275, 282, 290, 298, 302, 300, 293, 284, 275, 268)

  for (m in seq_along(days_2020)) {
    nd <- days_2020[m]
    t2m_vec <- make_hourly_t2m(nd, base_temps[m], 5)
    set.seed(m * 13)
    tp_vec <- pmax(0, rnorm(24 * nd, 0.00015, 0.00008))

    r_result   <- era5_to_monthly_r(t2m_vec, tp_vec, nd)
    cpp_result <- era5_to_monthly_cpp(
      matrix(t2m_vec, nrow = 1L),
      matrix(tp_vec,  nrow = 1L),
      nd
    )

    expect_equal(cpp_result$tas[1],    r_result$tas,    tolerance = tol,
                 label = sprintf("month %d tas", m))
    expect_equal(cpp_result$tasmax[1], r_result$tasmax, tolerance = tol,
                 label = sprintf("month %d tasmax", m))
    expect_equal(cpp_result$tasmin[1], r_result$tasmin, tolerance = tol,
                 label = sprintf("month %d tasmin", m))
    expect_equal(cpp_result$pr[1],     r_result$pr,     tolerance = tol,
                 label = sprintf("month %d pr", m))
  }
})
