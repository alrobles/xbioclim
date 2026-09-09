#' Quarterly / seasonal climate variables
#'
#' Computes six quarterly/seasonal variables for an arbitrary set of months.
#'
#' @param tas    Numeric vector (length 12), matrix (pixels x 12),
#'   [SpatRaster-class] with 12 layers, or [BioclimData-class].
#' @param tasmax Same form as `tas`: monthly maximum temperature.
#' @param tasmin Same form as `tas`: monthly minimum temperature.
#' @param pr     Same form as `tas`: monthly precipitation.
#' @param months Integer vector of 1-based month indices to include (e.g.
#'   `c(1, 2, 3)` for a fixed quarter, or `c(12, 1, 2)` for DJF).
#' @param na.rm  Logical. If `FALSE` (default), a single `NA` among the
#'   selected months makes all outputs `NA` for that pixel. If `TRUE`,
#'   missing months are skipped.
#' @param ...    Not currently used.
#'
#' @return
#' * Vector input: a named numeric vector of length 6.
#' * Matrix / [BioclimData-class] input: a numeric matrix (pixels x 6) with
#'   columns `tmean_s`, `tmax_max`, `tmin_min`, `trange`, `pr_tot`, `pr_cv`.
#' * [SpatRaster-class] input: a 6-layer [SpatRaster-class] with the same
#'   geometry and the same names.
#'
#' @details
#' The output variables are:
#' \itemize{
#'   \item `tmean_s` — mean `tas` over the selected months.
#'   \item `tmax_max` — maximum monthly `tasmax` among the selected months.
#'   \item `tmin_min` — minimum monthly `tasmin` among the selected months.
#'   \item `trange` — `tmax_max - tmin_min`.
#'   \item `pr_tot` — total `pr` over the selected months.
#'   \item `pr_cv` — coefficient of variation of monthly `pr`
#'     (`100 * sd_pop / mean`), matching the BIO15 scaling.
#' }
#'
#' Use [quarterly_fixed()] for the four standard fixed quarters and
#' [quarterly_rolling()] for rolling 3-month windows.
#'
#' @export
quarterly_variables <- function(tas, tasmax, tasmin, pr, months, na.rm = FALSE, ...) {
  months <- as.integer(months)
  if (any(months < 1L | months > 12L)) {
    stop("'months' must contain values between 1 and 12")
  }

  UseMethod("quarterly_variables", tas)
}

#' @rdname quarterly_variables
#' @export
quarterly_variables.default <- function(tas, tasmax, tasmin, pr, months, na.rm = FALSE, ...) {
  # Dispatch based on the actual class of tas.
  if (is.vector(tas)) {
    quarterly_variables.vector(tas, tasmax, tasmin, pr, months, na.rm, ...)
  } else if (is.matrix(tas)) {
    quarterly_variables.matrix(tas, tasmax, tasmin, pr, months, na.rm, ...)
  } else if (inherits(tas, "SpatRaster")) {
    quarterly_variables.SpatRaster(tas, tasmax, tasmin, pr, months, na.rm, ...)
  } else if (inherits(tas, "BioclimData")) {
    quarterly_variables.BioclimData(tas, tasmax, tasmin, pr, months, na.rm, ...)
  } else {
    stop("Unsupported input type for 'tas'")
  }
}

#' @exportS3Method quarterly_variables vector
quarterly_variables.vector <- function(tas, tasmax, tasmin, pr, months, na.rm = FALSE, ...) {
  m <- BioclimData(tas, tasmax, tasmin, pr)
  res <- quarterly_variables_cpp(m@tas, m@tasmax, m@tasmin, m@pr, months, na.rm)
  out <- res[1L, ]
  out
}

#' @exportS3Method quarterly_variables matrix
quarterly_variables.matrix <- function(tas, tasmax, tasmin, pr, months, na.rm = FALSE, ...) {
  if (!is.numeric(tas) || !is.numeric(tasmax) ||
      !is.numeric(tasmin) || !is.numeric(pr)) {
    stop("All inputs must be numeric")
  }
  if (ncol(tas) != 12L || ncol(tasmax) != 12L ||
      ncol(tasmin) != 12L || ncol(pr) != 12L) {
    stop("All inputs must have 12 columns (one per month)")
  }
  if (nrow(tas) != nrow(tasmax) || nrow(tas) != nrow(tasmin) ||
      nrow(tas) != nrow(pr)) {
    stop("All matrices must have the same number of rows")
  }
  quarterly_variables_cpp(tas, tasmax, tasmin, pr, months, na.rm)
}

#' @exportS3Method quarterly_variables BioclimData
quarterly_variables.BioclimData <- function(tas, tasmax, tasmin, pr, months, na.rm = FALSE, ...) {
  quarterly_variables_cpp(tas@tas, tas@tasmax, tas@tasmin, tas@pr, months, na.rm)
}

#' @exportS3Method quarterly_variables SpatRaster
quarterly_variables.SpatRaster <- function(tas, tasmax, tasmin, pr, months, na.rm = FALSE, ...) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required for SpatRaster inputs")
  }
  if (terra::nlyr(tas) != 12L || terra::nlyr(tasmax) != 12L ||
      terra::nlyr(tasmin) != 12L || terra::nlyr(pr) != 12L) {
    stop("All SpatRaster inputs must have 12 layers (one per month)")
  }

  vals_tas    <- terra::values(tas,    mat = TRUE)
  vals_tasmax <- terra::values(tasmax, mat = TRUE)
  vals_tasmin <- terra::values(tasmin, mat = TRUE)
  vals_pr     <- terra::values(pr,     mat = TRUE)

  out <- quarterly_variables_cpp(vals_tas, vals_tasmax, vals_tasmin,
                                 vals_pr, months, na.rm)

  template <- tas[[1]]
  r <- terra::rast(nrows = terra::nrow(template),
                   ncols = terra::ncol(template),
                   nlyrs = 6L,
                   extent = terra::ext(template),
                   crs = terra::crs(template))
  names(r) <- c("tmean_s", "tmax_max", "tmin_min", "trange", "pr_tot", "pr_cv")
  terra::values(r) <- out
  r
}

#' Fixed-quarter seasonal variables
#'
#' Convenience wrapper around [quarterly_variables()] for the four standard
#' fixed quarters.
#'
#' @param quarter Integer 1-4.
#' @param type    Character: `"meteorological"` (DJF, MAM, JJA, SON) or
#'   `"calendar"` (JFM, AMJ, JAS, OND).
#' @inheritParams quarterly_variables
#' @export
quarterly_fixed <- function(tas, tasmax, tasmin, pr, quarter,
                            type = c("meteorological", "calendar"),
                            na.rm = FALSE, ...) {
  type <- match.arg(type)
  if (quarter < 1L || quarter > 4L) {
    stop("'quarter' must be between 1 and 4")
  }

  if (type == "meteorological") {
    months <- switch(quarter,
                     `1` = c(12L, 1L, 2L),
                     `2` = c(3L, 4L, 5L),
                     `3` = c(6L, 7L, 8L),
                     `4` = c(9L, 10L, 11L))
  } else {
    months <- switch(quarter,
                     `1` = c(1L, 2L, 3L),
                     `2` = c(4L, 5L, 6L),
                     `3` = c(7L, 8L, 9L),
                     `4` = c(10L, 11L, 12L))
  }

  quarterly_variables(tas, tasmax, tasmin, pr, months, na.rm, ...)
}

#' Rolling-quarter seasonal variables
#'
#' Convenience wrapper around [quarterly_variables()] for a rolling 3-month
#' quarter starting at `start`.
#'
#' @param start Integer 1-12: starting month.
#' @inheritParams quarterly_variables
#' @export
quarterly_rolling <- function(tas, tasmax, tasmin, pr, start, na.rm = FALSE, ...) {
  if (start < 1L || start > 12L) {
    stop("'start' must be between 1 and 12")
  }
  months <- ((start - 1L):(start + 1L)) %% 12L + 1L
  quarterly_variables(tas, tasmax, tasmin, pr, months, na.rm, ...)
}
