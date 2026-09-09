# Test suite for bioclim_window and bioclim_rolling

mock_tas    <- 1:12
mock_tasmax <- 2:13
mock_tasmin <- 0:11
mock_pr     <- 1:12
mock_pr_rev <- 12:1

tol <- 1e-4

test_that("bioclim_window matches full bioclim over 12 months", {
  w <- bioclim_window(mock_tas, mock_tasmax, mock_tasmin, mock_pr,
                      months = 1:12)
  b <- bioclim(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_equal(unname(w), unname(b), tolerance = tol)
})

test_that("bioclim_rolling with window=3 matches full bioclim", {
  r <- bioclim_rolling(mock_tas, mock_tasmax, mock_tasmin, mock_pr,
                       window = 3)
  b <- bioclim(mock_tas, mock_tasmax, mock_tasmin, mock_pr)
  expect_equal(unname(r), unname(b), tolerance = tol)
})

test_that("bioclim_window base stats are over selected months", {
  w <- bioclim_window(mock_tas, mock_tasmax, mock_tasmin, mock_pr,
                      months = 1:2)
  expect_equal(unname(w["bio01"]), mean(mock_tas[1:2]), tolerance = tol)
  expect_equal(unname(w["bio05"]), max(mock_tasmax[1:2]), tolerance = tol)
  expect_equal(unname(w["bio06"]), min(mock_tasmin[1:2]), tolerance = tol)
  expect_equal(unname(w["bio12"]), sum(mock_pr[1:2]), tolerance = tol)
})

test_that("bioclim_window quarter vars are NA when window > length(months)", {
  w <- bioclim_window(mock_tas, mock_tasmax, mock_tasmin, mock_pr,
                      months = 1:2)
  expect_true(all(is.na(w[8:11])))
  expect_true(all(is.na(w[16:19])))
})

test_that("bioclim_window works with wrap-around months", {
  w <- bioclim_window(mock_tas, mock_tasmax, mock_tasmin, mock_pr,
                      months = c(12, 1, 2))
  expect_equal(unname(w["bio01"]), mean(c(mock_tas[12], mock_tas[1:2])),
               tolerance = tol)
  expect_equal(unname(w["bio12"]), sum(c(mock_pr[12], mock_pr[1:2])),
               tolerance = tol)
  expect_false(is.na(w["bio08"]))
})

test_that("bioclim_window resolves months from start/end strings", {
  w <- bioclim_window(mock_tas, mock_tasmax, mock_tasmin, mock_pr,
                      start = "06-07", end = "07-19")
  expect_equal(unname(w["bio01"]), mean(mock_tas[6:7]), tolerance = tol)
})

test_that("bioclim_rolling finds best 3-month windows", {
  r <- bioclim_rolling(mock_tas, mock_tasmax, mock_tasmin, mock_pr,
                       window = 3)
  expect_equal(unname(r["bio10"]), 11.0, tolerance = tol) # warmest quarter = Aug-Sep-Oct? tas mean 9+10+11=30/3=10? Wait
})

test_that("bioclim_window works on matrices", {
  n <- 5
  m_tas    <- matrix(rep(1:12, n), n, 12, byrow = TRUE)
  m_tasmax <- m_tas + 2
  m_tasmin <- m_tas - 2
  m_pr     <- matrix(rep(1:12, n), n, 12, byrow = TRUE)
  out <- bioclim_window(m_tas, m_tasmax, m_tasmin, m_pr, months = 6:9)
  expect_equal(dim(out), c(n, 19))
  expect_equal(colnames(out), c(paste0("bio", sprintf("%02d", 1:19))))
})

test_that("bioclim_window errors on invalid months", {
  expect_error(bioclim_window(mock_tas, mock_tasmax, mock_tasmin, mock_pr,
                              months = c(0, 1, 2)))
})

test_that("bioclim_rolling errors on invalid window", {
  expect_error(bioclim_rolling(mock_tas, mock_tasmax, mock_tasmin, mock_pr,
                               window = 1))
  expect_error(bioclim_rolling(mock_tas, mock_tasmax, mock_tasmin, mock_pr,
                               window = 12))
})
