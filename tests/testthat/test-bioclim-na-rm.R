test_that("na.rm = TRUE skips missing precipitation but keeps temperature", {
  tas  <- as.numeric(1:12)
  tmax <- tas + 1
  tmin <- tas - 1
  pr   <- as.numeric(c(2, 5, 3, 8, 15, 30, 60, 45, 20, 10, 5, 1))

  # Mask out the wettest month (July, position 7)
  pr_na <- pr
  pr_na[7] <- NA

  bd <- BioclimData(tas, tmax, tmin, pr_na)

  res <- bioclim(bd, na.rm = TRUE)[1, ]

  # bio01 depends only on tas, so unchanged
  expect_equal(res[["bio01"]], mean(tas), tolerance = 1e-10)

  # bio12 is the sum of non-NA precipitation
  expect_equal(res[["bio12"]], sum(pr_na, na.rm = TRUE), tolerance = 1e-10)

  # bio13 is the max of non-NA precipitation
  expect_equal(res[["bio13"]], max(pr_na, na.rm = TRUE), tolerance = 1e-10)

  # bio15 is 100 * CV of non-NA precipitation (population sd)
  pr_valid <- pr_na[!is.na(pr_na)]
  pr_mean  <- mean(pr_valid)
  pr_sd    <- sqrt(mean((pr_valid - pr_mean)^2))
  expect_equal(res[["bio15"]], 100 * pr_sd / pr_mean, tolerance = 1e-10)

  # With na.rm = FALSE the whole pixel should be NA
  res_strict <- bioclim(bd, na.rm = FALSE)[1, ]
  expect_true(all(is.na(res_strict)))
})

test_that("na.rm works for single-pixel vector bioclim()", {
  tas  <- as.numeric(1:12)
  tmax <- tas + 1
  tmin <- tas - 1
  pr   <- as.numeric(c(2, 5, 3, 8, 15, 30, 60, 45, 20, 10, 5, 1))

  pr_na <- pr
  pr_na[7] <- NA

  res <- bioclim(tas, tmax, tmin, pr_na, na.rm = TRUE)

  expect_equal(res[["bio01"]], mean(tas), tolerance = 1e-10)
  expect_equal(res[["bio12"]], sum(pr_na, na.rm = TRUE), tolerance = 1e-10)
  expect_equal(res[["bio13"]], max(pr_na, na.rm = TRUE), tolerance = 1e-10)
})

test_that("na.rm = TRUE handles missing temperature as well", {
  tas  <- as.numeric(1:12)
  tmax <- tas + 1
  tmin <- tas - 1
  pr   <- as.numeric(c(2, 5, 3, 8, 15, 30, 60, 45, 20, 10, 5, 1))

  tas_na <- tas
  tas_na[3] <- NA

  bd <- BioclimData(tas_na, tmax, tmin, pr)
  res <- bioclim(bd, na.rm = TRUE)[1, ]

  # bio01 is mean of non-NA tas
  expect_equal(res[["bio01"]], mean(tas_na, na.rm = TRUE), tolerance = 1e-10)

  # bio04 is 100 * population sd of non-NA tas
  tas_valid <- tas_na[!is.na(tas_na)]
  tas_mean  <- mean(tas_valid)
  tas_sd    <- sqrt(mean((tas_valid - tas_mean)^2))
  expect_equal(res[["bio04"]], 100 * tas_sd, tolerance = 1e-10)

  # bio12 (pr) should be unchanged because pr has no NA
  expect_equal(res[["bio12"]], sum(pr), tolerance = 1e-10)
})
