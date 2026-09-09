# Test suite for quarterly/seasonal climate variables

mock_tas    <- 1:12
mock_tasmax <- 2:13
mock_tasmin <- 0:11
mock_pr     <- 1:12

tol <- 1e-4

test_that("quarterly_variables works for a single pixel", {
  q <- quarterly_variables(mock_tas, mock_tasmax, mock_tasmin, mock_pr,
                           months = 1:3)
  expect_length(q, 6)
  expect_equal(names(q), c("tmean_s", "tmax_max", "tmin_min",
                           "trange", "pr_tot", "pr_cv"))
  expect_equal(q[["tmean_s"]], mean(mock_tas[1:3]), tolerance = tol)
  expect_equal(q[["tmax_max"]], max(mock_tasmax[1:3]), tolerance = tol)
  expect_equal(q[["tmin_min"]], min(mock_tasmin[1:3]), tolerance = tol)
  expect_equal(q[["trange"]], q[["tmax_max"]] - q[["tmin_min"]], tolerance = tol)
  expect_equal(q[["pr_tot"]], sum(mock_pr[1:3]), tolerance = tol)
})

test_that("quarterly_fixed maps meteorological quarters correctly", {
  q1 <- quarterly_fixed(mock_tas, mock_tasmax, mock_tasmin, mock_pr,
                        quarter = 1, type = "meteorological")
  q1_dir <- quarterly_variables(mock_tas, mock_tasmax, mock_tasmin, mock_pr,
                                months = c(12, 1, 2))
  expect_equal(q1, q1_dir, tolerance = tol)

  q3 <- quarterly_fixed(mock_tas, mock_tasmax, mock_tasmin, mock_pr,
                        quarter = 3, type = "meteorological")
  q3_dir <- quarterly_variables(mock_tas, mock_tasmax, mock_tasmin, mock_pr,
                                months = 6:8)
  expect_equal(q3, q3_dir, tolerance = tol)
})

test_that("quarterly_fixed maps calendar quarters correctly", {
  q4 <- quarterly_fixed(mock_tas, mock_tasmax, mock_tasmin, mock_pr,
                        quarter = 4, type = "calendar")
  q4_dir <- quarterly_variables(mock_tas, mock_tasmax, mock_tasmin, mock_pr,
                                months = 10:12)
  expect_equal(q4, q4_dir, tolerance = tol)
})

test_that("quarterly_rolling wraps around the year", {
  r12 <- quarterly_rolling(mock_tas, mock_tasmax, mock_tasmin, mock_pr,
                           start = 12)
  r12_dir <- quarterly_variables(mock_tas, mock_tasmax, mock_tasmin, mock_pr,
                                 months = c(12, 1, 2))
  expect_equal(r12, r12_dir, tolerance = tol)
})

test_that("quarterly_variables handles NA with default na.rm = FALSE", {
  tas <- 1:12
  tas[2] <- NA
  q <- quarterly_variables(tas, mock_tasmax, mock_tasmin, mock_pr,
                           months = 1:3)
  expect_true(all(is.na(q)))
})

test_that("quarterly_variables handles NA with na.rm = TRUE", {
  tas <- 1:12
  tas[2] <- NA
  q <- quarterly_variables(tas, mock_tasmax, mock_tasmin, mock_pr,
                           months = 1:3, na.rm = TRUE)
  expect_equal(q[["tmean_s"]], mean(c(1, 3)), tolerance = tol)
})

test_that("quarterly_variables works on matrices", {
  n <- 10
  m_tas    <- matrix(rep(1:12, n), n, 12, byrow = TRUE)
  m_tasmax <- m_tas + 2
  m_tasmin <- m_tas - 2
  m_pr     <- matrix(rep(1:12, n), n, 12, byrow = TRUE)

  out <- quarterly_variables(m_tas, m_tasmax, m_tasmin, m_pr, months = 4:6)
  expect_equal(dim(out), c(n, 6))
  expect_equal(colnames(out), c("tmean_s", "tmax_max", "tmin_min",
                                "trange", "pr_tot", "pr_cv"))
  expect_equal(unname(out[1, "tmean_s"]), mean(4:6), tolerance = tol)
  expect_equal(unname(out[1, "pr_tot"]), sum(4:6), tolerance = tol)
})

test_that("quarterly_variables works on BioclimData", {
  bd <- BioclimData(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  q <- quarterly_variables(bd, months = 1:3)
  expect_equal(dim(q), c(1, 6))
  expect_equal(unname(q[1, "tmean_s"]), mean(mock_tas[1:3]), tolerance = tol)
})

test_that("quarterly_variables errors on invalid months", {
  expect_error(quarterly_variables(mock_tas, mock_tasmax, mock_tasmin, mock_pr,
                                   months = c(0, 1, 2)))
  expect_error(quarterly_variables(mock_tas, mock_tasmax, mock_tasmin, mock_pr,
                                   months = c(10, 11, 13)))
})
