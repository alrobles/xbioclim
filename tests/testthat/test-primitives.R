# Tests for primitive helper functions

test_that("sd_pop computes population standard deviation", {
  # For 1:12, population SD = sqrt(sum((1:12 - 6.5)^2) / 12)
  expected <- sqrt(sum((1:12 - 6.5)^2) / 12)
  expect_equal(sd_pop(1:12), expected, tolerance = 1e-6)
})

test_that("sd_pop is zero for constant values", {
  expect_equal(sd_pop(rep(5, 12)), 0.0)
})

test_that("sd_pop differs from sample sd", {
  x <- 1:12
  # Population SD uses N, sample SD uses N-1
  expect_true(sd_pop(x) < sd(x))
})

test_that("sd_pop returns NaN for empty input", {
  expect_true(is.nan(sd_pop(numeric(0))))
})

test_that("rolling_quarter_sum computes correct sums", {
  x <- 1:12
  qs <- rolling_quarter_sum(x)
  expect_length(qs, 12)
  # First quarter: 1+2+3 = 6

expect_equal(qs[1], 6)
  # Last quarter (Dec-Jan-Feb): 12+1+2 = 15
  expect_equal(qs[12], 15)
  # Middle quarter (month 6): 6+7+8 = 21
  expect_equal(qs[6], 21)
})

test_that("rolling_quarter_sum handles circular wrapping", {
  x <- c(rep(0, 11), 100)
  qs <- rolling_quarter_sum(x)
  # Quarter starting at month 11: 0+100+0 = 100
  expect_equal(qs[11], 100)
  # Quarter starting at month 12: 100+0+0 = 100
  expect_equal(qs[12], 100)
  # Quarter starting at month 10: 0+0+100 = 100
  expect_equal(qs[10], 100)
})

test_that("rolling_quarter_mean computes correct means", {
  x <- 1:12
  qm <- rolling_quarter_mean(x)
  expect_length(qm, 12)
  # First quarter: (1+2+3)/3 = 2
  expect_equal(qm[1], 2)
  # Quarter at month 10: (10+11+12)/3 = 11
  expect_equal(qm[10], 11)
})

test_that("quarter_argmax returns correct starting month", {
  x <- 1:12
  # Maximum sum quarter starts at month 10 (10+11+12 = 33)
  expect_equal(quarter_argmax(x), 10)
})

test_that("quarter_argmin returns correct starting month", {
  x <- 1:12
  # Minimum sum quarter starts at month 1 (1+2+3 = 6)
  expect_equal(quarter_argmin(x), 1)
})

test_that("quarter_values extracts correct months", {
  x <- 1:12
  # Quarter starting at month 1: months 1, 2, 3
  expect_equal(quarter_values(x, 1), c(1, 2, 3))
  # Quarter starting at month 10: months 10, 11, 12
  expect_equal(quarter_values(x, 10), c(10, 11, 12))
  # Quarter starting at month 12 (wraps): months 12, 1, 2
  expect_equal(quarter_values(x, 12), c(12, 1, 2))
})

test_that("validate_monthly rejects non-numeric", {
  expect_error(validate_monthly(letters[1:12], "test"), "must be numeric")
})

test_that("validate_monthly rejects wrong length", {
  expect_error(validate_monthly(1:6, "test"), "must have length 12")
  expect_error(validate_monthly(1:24, "test"), "must have length 12")
})

test_that("validate_monthly accepts valid input", {
  expect_silent(validate_monthly(1:12, "test"))
  expect_silent(validate_monthly(as.numeric(1:12), "test"))
})
