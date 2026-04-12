#' Compute Bioclimatic Variables from Monthly Climate Data
#'
#' @description
#' Computes the 19 standard bioclimatic variables (BIO01-BIO19) from monthly
#' climate data following the WorldClim specification. This is an R
#' implementation of the xbioclim C++ library.
#'
#' @details
#' The 19 bioclimatic variables are:
#' \itemize{
#'   \item BIO01: Mean Annual Temperature
#'   \item BIO02: Mean Diurnal Range (mean of monthly (tasmax - tasmin))
#'   \item BIO03: Isothermality (100 * BIO02 / BIO07)
#'   \item BIO04: Temperature Seasonality (100 * population SD of monthly tas)
#'   \item BIO05: Max Temperature of Warmest Month
#'   \item BIO06: Min Temperature of Coldest Month
#'   \item BIO07: Temperature Annual Range (BIO05 - BIO06)
#'   \item BIO08: Mean Temperature of Wettest Quarter
#'   \item BIO09: Mean Temperature of Driest Quarter
#'   \item BIO10: Mean Temperature of Warmest Quarter
#'   \item BIO11: Mean Temperature of Coldest Quarter
#'   \item BIO12: Annual Precipitation
#'   \item BIO13: Precipitation of Wettest Month
#'   \item BIO14: Precipitation of Driest Month
#'   \item BIO15: Precipitation Seasonality (CV)
#'   \item BIO16: Precipitation of Wettest Quarter
#'   \item BIO17: Precipitation of Driest Quarter
#'   \item BIO18: Precipitation of Warmest Quarter
#'   \item BIO19: Precipitation of Coldest Quarter
#' }
#'
#' @name bioclim-variables
NULL

#' @describeIn bioclim-variables BIO01 - Mean Annual Temperature
#' @param tas Numeric vector of length 12: monthly mean temperature.
#' @return A single numeric value.
#' @export
#' @examples
#' tas <- 1:12
#' bio01(tas)
bio01 <- function(tas) {
  validate_monthly(tas, "tas")
  mean(tas)
}

#' @describeIn bioclim-variables BIO02 - Mean Diurnal Range
#' @param tasmax Numeric vector of length 12: monthly maximum temperature.
#' @param tasmin Numeric vector of length 12: monthly minimum temperature.
#' @return A single numeric value.
#' @export
#' @examples
#' tasmax <- 3:14
#' tasmin <- 0:11
#' bio02(tasmax, tasmin)
bio02 <- function(tasmax, tasmin) {
  validate_monthly(tasmax, "tasmax")
  validate_monthly(tasmin, "tasmin")
  mean(tasmax - tasmin)
}

#' @describeIn bioclim-variables BIO03 - Isothermality (100 * BIO02 / BIO07)
#' @param tasmax Numeric vector of length 12: monthly maximum temperature.
#' @param tasmin Numeric vector of length 12: monthly minimum temperature.
#' @return A single numeric value.
#' @export
#' @examples
#' tasmax <- 3:14
#' tasmin <- 0:11
#' bio03(tasmax, tasmin)
bio03 <- function(tasmax, tasmin) {
  b02 <- bio02(tasmax, tasmin)
  b07 <- bio07(tasmax, tasmin)
  if (b07 == 0) return(NaN)
  100 * b02 / b07
}

#' @describeIn bioclim-variables BIO04 - Temperature Seasonality
#' @param tas Numeric vector of length 12: monthly mean temperature.
#' @return A single numeric value.
#' @export
#' @examples
#' tas <- 1:12
#' bio04(tas)
bio04 <- function(tas) {
  validate_monthly(tas, "tas")
  100 * sd_pop(tas)
}

#' @describeIn bioclim-variables BIO05 - Max Temperature of Warmest Month
#' @param tasmax Numeric vector of length 12: monthly maximum temperature.
#' @return A single numeric value.
#' @export
#' @examples
#' tasmax <- 3:14
#' bio05(tasmax)
bio05 <- function(tasmax) {
  validate_monthly(tasmax, "tasmax")
  max(tasmax)
}

#' @describeIn bioclim-variables BIO06 - Min Temperature of Coldest Month
#' @param tasmin Numeric vector of length 12: monthly minimum temperature.
#' @return A single numeric value.
#' @export
#' @examples
#' tasmin <- 0:11
#' bio06(tasmin)
bio06 <- function(tasmin) {
  validate_monthly(tasmin, "tasmin")
  min(tasmin)
}

#' @describeIn bioclim-variables BIO07 - Temperature Annual Range (BIO05 - BIO06)
#' @param tasmax Numeric vector of length 12: monthly maximum temperature.
#' @param tasmin Numeric vector of length 12: monthly minimum temperature.
#' @return A single numeric value.
#' @export
#' @examples
#' tasmax <- 3:14
#' tasmin <- 0:11
#' bio07(tasmax, tasmin)
bio07 <- function(tasmax, tasmin) {
  bio05(tasmax) - bio06(tasmin)
}

#' @describeIn bioclim-variables BIO08 - Mean Temperature of Wettest Quarter
#' @param tas Numeric vector of length 12: monthly mean temperature.
#' @param pr Numeric vector of length 12: monthly precipitation.
#' @return A single numeric value.
#' @export
#' @examples
#' tas <- 1:12
#' pr <- 1:12
#' bio08(tas, pr)
bio08 <- function(tas, pr) {
  validate_monthly(tas, "tas")
  validate_monthly(pr, "pr")
  wet_start <- quarter_argmax(pr)
  mean(quarter_values(tas, wet_start))
}

#' @describeIn bioclim-variables BIO09 - Mean Temperature of Driest Quarter
#' @param tas Numeric vector of length 12: monthly mean temperature.
#' @param pr Numeric vector of length 12: monthly precipitation.
#' @return A single numeric value.
#' @export
#' @examples
#' tas <- 1:12
#' pr <- 1:12
#' bio09(tas, pr)
bio09 <- function(tas, pr) {
  validate_monthly(tas, "tas")
  validate_monthly(pr, "pr")
  dry_start <- quarter_argmin(pr)
  mean(quarter_values(tas, dry_start))
}

#' @describeIn bioclim-variables BIO10 - Mean Temperature of Warmest Quarter
#' @param tas Numeric vector of length 12: monthly mean temperature.
#' @return A single numeric value.
#' @export
#' @examples
#' tas <- 1:12
#' bio10(tas)
bio10 <- function(tas) {
  validate_monthly(tas, "tas")
  warm_start <- quarter_argmax(tas)
  mean(quarter_values(tas, warm_start))
}

#' @describeIn bioclim-variables BIO11 - Mean Temperature of Coldest Quarter
#' @param tas Numeric vector of length 12: monthly mean temperature.
#' @return A single numeric value.
#' @export
#' @examples
#' tas <- 1:12
#' bio11(tas)
bio11 <- function(tas) {
  validate_monthly(tas, "tas")
  cold_start <- quarter_argmin(tas)
  mean(quarter_values(tas, cold_start))
}

#' @describeIn bioclim-variables BIO12 - Annual Precipitation
#' @param pr Numeric vector of length 12: monthly precipitation.
#' @return A single numeric value.
#' @export
#' @examples
#' pr <- 1:12
#' bio12(pr)
bio12 <- function(pr) {
  validate_monthly(pr, "pr")
  sum(pr)
}

#' @describeIn bioclim-variables BIO13 - Precipitation of Wettest Month
#' @param pr Numeric vector of length 12: monthly precipitation.
#' @return A single numeric value.
#' @export
#' @examples
#' pr <- 1:12
#' bio13(pr)
bio13 <- function(pr) {
  validate_monthly(pr, "pr")
  max(pr)
}

#' @describeIn bioclim-variables BIO14 - Precipitation of Driest Month
#' @param pr Numeric vector of length 12: monthly precipitation.
#' @return A single numeric value.
#' @export
#' @examples
#' pr <- 1:12
#' bio14(pr)
bio14 <- function(pr) {
  validate_monthly(pr, "pr")
  min(pr)
}

#' @describeIn bioclim-variables BIO15 - Precipitation Seasonality (CV)
#' @param pr Numeric vector of length 12: monthly precipitation.
#' @return A single numeric value.
#' @export
#' @examples
#' pr <- 1:12
#' bio15(pr)
bio15 <- function(pr) {
  validate_monthly(pr, "pr")
  pr_mean <- mean(pr)
  if (pr_mean == 0) return(NaN)
  100 * sd_pop(pr) / pr_mean
}

#' @describeIn bioclim-variables BIO16 - Precipitation of Wettest Quarter
#' @param pr Numeric vector of length 12: monthly precipitation.
#' @return A single numeric value.
#' @export
#' @examples
#' pr <- 1:12
#' bio16(pr)
bio16 <- function(pr) {
  validate_monthly(pr, "pr")
  wet_start <- quarter_argmax(pr)
  sum(quarter_values(pr, wet_start))
}

#' @describeIn bioclim-variables BIO17 - Precipitation of Driest Quarter
#' @param pr Numeric vector of length 12: monthly precipitation.
#' @return A single numeric value.
#' @export
#' @examples
#' pr <- 1:12
#' bio17(pr)
bio17 <- function(pr) {
  validate_monthly(pr, "pr")
  dry_start <- quarter_argmin(pr)
  sum(quarter_values(pr, dry_start))
}

#' @describeIn bioclim-variables BIO18 - Precipitation of Warmest Quarter
#' @param tas Numeric vector of length 12: monthly mean temperature.
#' @param pr Numeric vector of length 12: monthly precipitation.
#' @return A single numeric value.
#' @export
#' @examples
#' tas <- 1:12
#' pr <- 1:12
#' bio18(tas, pr)
bio18 <- function(tas, pr) {
  validate_monthly(tas, "tas")
  validate_monthly(pr, "pr")
  warm_start <- quarter_argmax(tas)
  sum(quarter_values(pr, warm_start))
}

#' @describeIn bioclim-variables BIO19 - Precipitation of Coldest Quarter
#' @param tas Numeric vector of length 12: monthly mean temperature.
#' @param pr Numeric vector of length 12: monthly precipitation.
#' @return A single numeric value.
#' @export
#' @examples
#' tas <- 1:12
#' pr <- 1:12
#' bio19(tas, pr)
bio19 <- function(tas, pr) {
  validate_monthly(tas, "tas")
  validate_monthly(pr, "pr")
  cold_start <- quarter_argmin(tas)
  sum(quarter_values(pr, cold_start))
}

#' Compute All 19 Bioclimatic Variables
#'
#' Computes all 19 bioclimatic variables at once from monthly climate data.
#'
#' @param tas Numeric vector of length 12: monthly mean temperature.
#' @param tasmax Numeric vector of length 12: monthly maximum temperature.
#' @param tasmin Numeric vector of length 12: monthly minimum temperature.
#' @param pr Numeric vector of length 12: monthly precipitation.
#'
#' @return A named numeric vector of length 19 (bio01 through bio19).
#'
#' @export
#' @examples
#' tas <- 1:12
#' tasmax <- 2:13
#' tasmin <- 0:11
#' pr <- 1:12
#' bioclim(tas, tasmax, tasmin, pr)
bioclim <- function(tas, tasmax, tasmin, pr) {
  validate_monthly(tas, "tas")
  validate_monthly(tasmax, "tasmax")
  validate_monthly(tasmin, "tasmin")
  validate_monthly(pr, "pr")

  # Temperature basics
  b01 <- mean(tas)
  b02 <- mean(tasmax - tasmin)
  b05 <- max(tasmax)
  b06 <- min(tasmin)
  b07 <- b05 - b06
  b03 <- if (b07 == 0) NaN else 100 * b02 / b07
  b04 <- 100 * sd_pop(tas)

  # Quarter indices
  wet_start  <- quarter_argmax(pr)
  dry_start  <- quarter_argmin(pr)
  warm_start <- quarter_argmax(tas)
  cold_start <- quarter_argmin(tas)

  # Temperature quarter means
  b08 <- mean(quarter_values(tas, wet_start))
  b09 <- mean(quarter_values(tas, dry_start))
  b10 <- mean(quarter_values(tas, warm_start))
  b11 <- mean(quarter_values(tas, cold_start))

  # Precipitation basics
  b12 <- sum(pr)
  b13 <- max(pr)
  b14 <- min(pr)
  pr_mean <- mean(pr)
  b15 <- if (pr_mean == 0) NaN else 100 * sd_pop(pr) / pr_mean

  # Precipitation quarter sums
  b16 <- sum(quarter_values(pr, wet_start))
  b17 <- sum(quarter_values(pr, dry_start))
  b18 <- sum(quarter_values(pr, warm_start))
  b19 <- sum(quarter_values(pr, cold_start))

  c(bio01 = b01, bio02 = b02, bio03 = b03, bio04 = b04,
    bio05 = b05, bio06 = b06, bio07 = b07,
    bio08 = b08, bio09 = b09, bio10 = b10, bio11 = b11,
    bio12 = b12, bio13 = b13, bio14 = b14, bio15 = b15,
    bio16 = b16, bio17 = b17, bio18 = b18, bio19 = b19)
}
