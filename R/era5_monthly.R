#' Convert ERA5-Land Hourly Data to CHELSA-Compatible Monthly Variables
#'
#' @description
#' Aggregates ERA5-Land hourly reanalysis data to CHELSA-compatible monthly
#' climate variables.  The four output variables (\code{tas}, \code{tasmax},
#' \code{tasmin}, \code{pr}) can be fed directly into
#' \code{\link{bioclim}} or \code{\link{bioclim_raster}} to compute
#' bioclimatic variables BIO01–BIO19.
#'
#' @details
#' \strong{Temperature aggregation.}
#' ERA5-Land provides instantaneous 2-m temperature (\code{t2m}) at hourly
#' resolution in Kelvin.
#' The hourly values are first grouped into calendar days (24 hours each):
#' \itemize{
#'   \item \code{tas}    — monthly mean of daily means
#'   \item \code{tasmax} — monthly mean of daily maxima
#'   \item \code{tasmin} — monthly mean of daily minima
#' }
#'
#' \strong{Precipitation aggregation.}
#' ERA5-Land provides total precipitation (\code{tp}) as hourly accumulations
#' in metres of water equivalent.  The hourly values are summed over the month
#' and converted to \eqn{\mathrm{kg\,m^{-2}\,month^{-1}}}{kg m-2 month-1}
#' (= mm) by multiplying by 1000.
#'
#' \strong{Unit conventions.}
#' By default, temperatures are returned in Kelvin to match the CHELSA
#' convention.
#' Set \code{to_celsius = TRUE} to obtain degrees Celsius instead (common for
#' WorldClim-style bioclimatic variables).
#'
#' @name era5-monthly
NULL

# ── Reference implementations (pure R, @keywords internal) ──────────────────

#' @describeIn era5-monthly Reference implementation: hourly t2m → monthly
#'   temperature statistics.  Pure-R, single-pixel.
#' @param hourly_t2m Numeric vector of hourly 2-m temperatures (K).
#'   Length must equal \code{24 * n_days}.
#' @param n_days Integer: number of days in the month.
#' @param to_celsius Logical: convert K → °C? Default \code{FALSE}.
#' @return Named list: \code{tas}, \code{tasmax}, \code{tasmin}.
#' @keywords internal
era5_t2m_to_monthly_r <- function(hourly_t2m, n_days, to_celsius = FALSE) {
  n_hours <- length(hourly_t2m)
  stopifnot(n_hours == 24L * n_days)

  offset <- if (to_celsius) -273.15 else 0.0

  daily_mean <- numeric(n_days)
  daily_max  <- numeric(n_days)
  daily_min  <- numeric(n_days)

  for (d in seq_len(n_days)) {
    idx <- ((d - 1L) * 24L + 1L):(d * 24L)
    vals <- hourly_t2m[idx]
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0L) {
      daily_mean[d] <- NaN
      daily_max[d]  <- NaN
      daily_min[d]  <- NaN
    } else {
      daily_mean[d] <- mean(vals)
      daily_max[d]  <- max(vals)
      daily_min[d]  <- min(vals)
    }
  }

  list(
    tas    = mean(daily_mean) + offset,
    tasmax = mean(daily_max)  + offset,
    tasmin = mean(daily_min)  + offset
  )
}

#' @describeIn era5-monthly Reference implementation: hourly tp → monthly
#'   precipitation.  Pure-R, single-pixel.
#' @param hourly_tp Numeric vector of hourly total precipitation (m).
#' @return Numeric scalar: monthly precipitation in mm (kg m-2).
#' @keywords internal
era5_tp_to_monthly_r <- function(hourly_tp) {
  vals <- hourly_tp[!is.na(hourly_tp)]
  vals <- vals[vals > 0]
  sum(vals) * 1000.0
}

# ── Unified reference implementation ──────────────────────────────────────

#' @describeIn era5-monthly Unified reference implementation: hourly t2m + tp →
#'   monthly tas, tasmax, tasmin, pr.  Pure-R, single-pixel.
#' @param hourly_t2m Numeric vector of hourly 2-m temperatures (K).
#'   Length must equal \code{24 * n_days}.
#' @param hourly_tp Numeric vector of hourly total precipitation (m).
#' @param n_days Integer: number of days in the month.
#' @param to_celsius Logical: convert K → °C? Default \code{FALSE}.
#' @return Named list: \code{tas}, \code{tasmax}, \code{tasmin}, \code{pr}.
#' @keywords internal
era5_to_monthly_r <- function(hourly_t2m, hourly_tp, n_days,
                              to_celsius = FALSE) {
  temp <- era5_t2m_to_monthly_r(hourly_t2m, n_days, to_celsius)
  pr   <- era5_tp_to_monthly_r(hourly_tp)
  list(
    tas    = temp$tas,
    tasmax = temp$tasmax,
    tasmin = temp$tasmin,
    pr     = pr
  )
}

# ── Exported wrappers ──────────────────────────────────────────────────────

#' @describeIn era5-monthly Convert hourly 2-m temperature to monthly
#'   statistics (C++ backend).
#'
#' @param hourly_t2m Numeric matrix (n_pixels × n_hours) of hourly 2-m
#'   temperatures (K).  Or a numeric vector for a single pixel.
#' @param n_days Integer: number of days in the month (e.g. 31 for January).
#' @param to_celsius Logical: convert output from Kelvin to Celsius?
#'   Default \code{FALSE}.
#' @param ncores Integer: OpenMP thread count. Default \code{1L}.
#'
#' @return Named list with \code{tas}, \code{tasmax}, \code{tasmin} — each a
#'   numeric vector of length n_pixels.
#'
#' @examples
#' # Single pixel: 3 days of hourly data (72 hours) at ~285 K
#' set.seed(42)
#' hourly <- 285 + cumsum(rnorm(72, 0, 0.5))
#' result <- era5_t2m_to_monthly(hourly, n_days = 3L)
#' result$tas     # monthly mean temperature (K)
#' result$tasmax  # monthly mean of daily maxima (K)
#' result$tasmin  # monthly mean of daily minima (K)
#'
#' @export
era5_t2m_to_monthly <- function(hourly_t2m, n_days,
                                to_celsius = FALSE,
                                ncores = 1L) {
  if (is.numeric(hourly_t2m) && is.null(dim(hourly_t2m))) {
    # Single pixel: use reference implementation
    return(era5_t2m_to_monthly_r(hourly_t2m, n_days, to_celsius))
  }
  if (!is.matrix(hourly_t2m)) {
    hourly_t2m <- as.matrix(hourly_t2m)
  }
  stopifnot(ncol(hourly_t2m) == 24L * n_days)
  era5_t2m_to_monthly_cpp(hourly_t2m, n_days, to_celsius, ncores)
}

#' @describeIn era5-monthly Convert hourly total precipitation to monthly
#'   total (C++ backend).
#'
#' @param hourly_tp Numeric matrix (n_pixels × n_hours) of hourly total
#'   precipitation (m).  Or a numeric vector for a single pixel.
#' @param ncores Integer: OpenMP thread count. Default \code{1L}.
#'
#' @return Numeric vector of length n_pixels: monthly precipitation in
#'   kg m-2 (mm).
#'
#' @examples
#' # Single pixel: 3 days of hourly precipitation (72 hours)
#' set.seed(42)
#' hourly_tp <- pmax(0, rnorm(72, 0.0001, 0.00005))
#' era5_tp_to_monthly(hourly_tp)  # total in mm
#'
#' @export
era5_tp_to_monthly <- function(hourly_tp, ncores = 1L) {
  if (is.numeric(hourly_tp) && is.null(dim(hourly_tp))) {
    return(era5_tp_to_monthly_r(hourly_tp))
  }
  if (!is.matrix(hourly_tp)) {
    hourly_tp <- as.matrix(hourly_tp)
  }
  era5_tp_to_monthly_cpp(hourly_tp, ncores)
}

# ── Unified exported wrapper ───────────────────────────────────────────────

#' @describeIn era5-monthly Convert ERA5-Land hourly t2m and tp to the four
#'   CHELSA-compatible monthly climate variables in a single call.
#'
#' @param hourly_t2m Numeric matrix (n_pixels × n_hours) of hourly 2-m
#'   temperatures (K), or a numeric vector for a single pixel.
#' @param hourly_tp Numeric matrix (n_pixels × n_hours) of hourly total
#'   precipitation (m), or a numeric vector for a single pixel.
#' @param n_days Integer: number of days in the month.
#' @param to_celsius Logical: convert temperatures from Kelvin to Celsius?
#'   Default \code{FALSE}.
#' @param ncores Integer: OpenMP thread count. Default \code{1L}.
#'
#' @return Named list with \code{tas}, \code{tasmax}, \code{tasmin}, \code{pr}
#'   — each a numeric vector of length n_pixels (or scalar for single pixel).
#'
#' @examples
#' # Single pixel: 3 days of synthetic hourly data
#' n_days <- 3L
#' set.seed(42)
#' hourly_t2m <- 285 + 5 * sin(2 * pi * (seq(0, 71) - 4) / 24)
#' hourly_tp  <- pmax(0, rnorm(72, 0.0001, 0.00005))
#' result <- era5_to_monthly(hourly_t2m, hourly_tp, n_days)
#' result$tas     # monthly mean temperature (K)
#' result$tasmax  # monthly mean of daily maxima (K)
#' result$tasmin  # monthly mean of daily minima (K)
#' result$pr      # monthly precipitation (mm)
#'
#' @export
era5_to_monthly <- function(hourly_t2m, hourly_tp, n_days,
                            to_celsius = FALSE,
                            ncores = 1L) {
  # Single-pixel vector path: use R reference
  if (is.numeric(hourly_t2m) && is.null(dim(hourly_t2m))) {
    return(era5_to_monthly_r(hourly_t2m, hourly_tp, n_days, to_celsius))
  }
  # Matrix path: use C++ backend
  if (!is.matrix(hourly_t2m)) hourly_t2m <- as.matrix(hourly_t2m)
  if (!is.matrix(hourly_tp))  hourly_tp  <- as.matrix(hourly_tp)
  stopifnot(ncol(hourly_t2m) == 24L * n_days)
  era5_to_monthly_cpp(hourly_t2m, hourly_tp, n_days, to_celsius, ncores)
}

# ── Helper: days in month ──────────────────────────────────────────────────

#' Number of days in a given month/year
#'
#' @param year  Integer year.
#' @param month Integer month (1–12).
#' @return Integer number of days.
#' @keywords internal
days_in_month <- function(year, month) {
  if (month == 12L) {
    next_year  <- year + 1L
    next_month <- 1L
  } else {
    next_year  <- year
    next_month <- month + 1L
  }
  as.integer(as.Date(
    paste(next_year, next_month, 1L, sep = "-"),
    format = "%Y-%m-%d"
  ) - as.Date(
    paste(year, month, 1L, sep = "-"),
    format = "%Y-%m-%d"
  ))
}
