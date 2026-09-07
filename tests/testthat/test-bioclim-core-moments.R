test_that("bioclim core moments: bio01, bio04 and bio15 match R mean/sd", {
  # Controlled 12-month ramp that has non-trivial mean and variance.
  tas  <- as.numeric(1:12)
  tmax <- tas + 1
  tmin <- tas - 1
  pr   <- as.numeric(c(2, 5, 3, 8, 15, 30, 60, 45, 20, 10, 5, 1))

  result <- bioclim(tas, tmax, tmin, pr)

  # bio01 = annual mean temperature
  expect_equal(result[["bio01"]], mean(tas), tolerance = 1e-10,
               label = "bio01 == mean(tas)")

  # bio04 = 100 * population standard deviation of monthly mean temperature
  sd_pop <- sqrt(mean((tas - mean(tas))^2))
  expect_equal(result[["bio04"]], 100 * sd_pop, tolerance = 1e-10,
               label = "bio04 == 100 * sd_pop(tas)")

  # bio15 = 100 * cv of precipitation (population sd / mean)
  pr_mean <- mean(pr)
  pr_sd   <- sqrt(mean((pr - pr_mean)^2))
  expect_equal(result[["bio15"]], 100 * pr_sd / pr_mean, tolerance = 1e-10,
               label = "bio15 == 100 * cv(pr)")
})
