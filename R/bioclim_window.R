#' Bioclimatic variables over an arbitrary window of months
#'
#' Computes the 19 standard bioclimatic variables over a user-defined
#' subset of months.  This is useful when a species' relevant season
#' does not coincide with fixed quarters (e.g. eBird Status & Trends
#' breeding windows).
#'
#' @param tas    Numeric vector (length 12), matrix (pixels x 12),
#'   [SpatRaster-class] with 12 layers, or [BioclimData-class].
#' @param tasmax Same form as `tas`: monthly maximum temperature.
#' @param tasmin Same form as `tas`: monthly minimum temperature.
#' @param pr     Same form as `tas`: monthly precipitation.
#' @param months Integer vector of 1-based months in the window, e.g.
#'   `c(6, 7)` for Jun-Jul or `c(12, 1, 2)` for DJF.  If `start` and
#'   `end` are supplied, `months` is derived automatically.
#' @param start,end Optional character strings in `"MM-DD"` format.
#'   When both are supplied, the months that overlap the day range are
#'   selected.  Wrap-around across the year end is supported.
#' @param window Integer: length of the internal rolling sub-window used
#'   to compute `bio08`..`bio19`.  Defaults to 3 (quarterly).  Must be
#'   `<= length(months)` for the rolling variables to be non-`NA`.
#' @param na.rm  Logical. If `FALSE` (default), a single `NA` among the
#'   selected months makes all outputs `NA` for that pixel. If `TRUE`,
#'   missing months are skipped.
#' @param ...    Not currently used.
#'
#' @return
#' * Vector input: a named numeric vector of length 19 (`bio01`..`bio19`).
#' * Matrix / [BioclimData-class] input: a numeric matrix (pixels x 19).
#' * [SpatRaster-class] input: a 19-layer [SpatRaster-class].
#'
#' @details
#' The base variables (`bio01`..`bio07`, `bio12`..`bio15`) are computed
#' directly over the selected months.  The rolling variables
#' (`bio08`..`bio11`, `bio16`..`bio19`) are computed over the best
#' contiguous `window`-month period that is fully contained in `months`.
#' If `length(months) < window` or no such contiguous period exists,
#' those columns are `NA`.
#'
#' @seealso [bioclim_rolling()] for rolling optimal windows over the full
#'   year with a free `window` length.
#'
#' @export
bioclim_window <- function(tas, tasmax, tasmin, pr,
                           months = NULL,
                           start = NULL, end = NULL,
                           window = 3L, na.rm = FALSE, ...) {
  months <- resolve_months(tas, months, start, end)
  months <- as.integer(months)
  if (any(months < 1L | months > 12L)) {
    stop("'months' must contain values between 1 and 12")
  }

  if (is.vector(tas)) {
    m <- BioclimData(tas, tasmax, tasmin, pr)
    res <- bioclim_window_cpp(m@tas, m@tasmax, m@tasmin, m@pr,
                              months, window, na.rm)
    res[1L, ]
  } else if (is.matrix(tas)) {
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
    bioclim_window_cpp(tas, tasmax, tasmin, pr, months, window, na.rm)
  } else if (inherits(tas, "SpatRaster")) {
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

    out <- bioclim_window_cpp(vals_tas, vals_tasmax, vals_tasmin,
                              vals_pr, months, window, na.rm)

    template <- tas[[1]]
    r <- terra::rast(nrows = terra::nrow(template),
                     ncols = terra::ncol(template),
                     nlyrs = 19L,
                     extent = terra::ext(template),
                     crs = terra::crs(template))
    names(r) <- c("bio01", "bio02", "bio03", "bio04", "bio05",
                  "bio06", "bio07", "bio08", "bio09", "bio10",
                  "bio11", "bio12", "bio13", "bio14", "bio15",
                  "bio16", "bio17", "bio18", "bio19")
    terra::values(r) <- out
    r
  } else if (inherits(tas, "BioclimData")) {
    bioclim_window_cpp(tas@tas, tas@tasmax, tas@tasmin, tas@pr,
                       months, window, na.rm)
  } else {
    stop("Unsupported input type for 'tas'")
  }
}

#' Rolling optimal-window bioclimatic variables
#'
#' Computes the 19 bioclimatic variables using a rolling window of arbitrary
#' length over the full 12 months.  The base variables (`bio01`..`bio07`,
#' `bio12`..`bio15`) are annual, while the rolling variables
#' (`bio08`..`bio11`, `bio16`..`bio19`) use the best `window`-month period.
#' With `window = 3` this is equivalent to the standard BIO08-BIO19.
#'
#' @param tas    Numeric vector (length 12), matrix (pixels x 12),
#'   [SpatRaster-class] with 12 layers, or [BioclimData-class].
#' @param tasmax Same form as `tas`: monthly maximum temperature.
#' @param tasmin Same form as `tas`: monthly minimum temperature.
#' @param pr     Same form as `tas`: monthly precipitation.
#' @param window Integer: length of the rolling window in months (2-11).
#'   Default is 3.
#' @param na.rm  Logical. If `FALSE` (default), a single `NA` makes all
#'   outputs `NA` for that pixel. If `TRUE`, missing months are skipped.
#' @param ...    Not currently used.
#'
#' @return Named vector (length 19), matrix (pixels x 19), or 19-layer
#'   [SpatRaster-class].
#'
#' @export
bioclim_rolling <- function(tas, tasmax, tasmin, pr,
                            window = 3L, na.rm = FALSE, ...) {
  window <- as.integer(window)
  if (window < 2L || window > 11L) {
    stop("'window' must be between 2 and 11")
  }

  if (is.vector(tas)) {
    m <- BioclimData(tas, tasmax, tasmin, pr)
    res <- bioclim_rolling_cpp(m@tas, m@tasmax, m@tasmin, m@pr,
                               window, na.rm)
    res[1L, ]
  } else if (is.matrix(tas)) {
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
    bioclim_rolling_cpp(tas, tasmax, tasmin, pr, window, na.rm)
  } else if (inherits(tas, "SpatRaster")) {
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

    out <- bioclim_rolling_cpp(vals_tas, vals_tasmax, vals_tasmin,
                               vals_pr, window, na.rm)

    template <- tas[[1]]
    r <- terra::rast(nrows = terra::nrow(template),
                     ncols = terra::ncol(template),
                     nlyrs = 19L,
                     extent = terra::ext(template),
                     crs = terra::crs(template))
    names(r) <- c("bio01", "bio02", "bio03", "bio04", "bio05",
                  "bio06", "bio07", "bio08", "bio09", "bio10",
                  "bio11", "bio12", "bio13", "bio14", "bio15",
                  "bio16", "bio17", "bio18", "bio19")
    terra::values(r) <- out
    r
  } else if (inherits(tas, "BioclimData")) {
    bioclim_rolling_cpp(tas@tas, tas@tasmax, tas@tasmin, tas@pr,
                        window, na.rm)
  } else {
    stop("Unsupported input type for 'tas'")
  }
}

# ── helper: resolve months from dates or explicit vector ─────────────────────

resolve_months <- function(tas, months, start, end) {
  if (!is.null(months)) {
    return(as.integer(months))
  }
  if (is.null(start) || is.null(end)) {
    stop("Either 'months' or both 'start' and 'end' must be supplied")
  }

  parse_month <- function(s) {
    if (!grepl("^[0-9]{2}-[0-9]{2}$", s)) {
      stop("Date strings must be in 'MM-DD' format")
    }
    as.integer(substr(s, 1, 2))
  }

  m_start <- parse_month(start)
  m_end   <- parse_month(end)

  if (m_start < 1L || m_start > 12L || m_end < 1L || m_end > 12L) {
    stop("Parsed months must be between 1 and 12")
  }

  if (m_start <= m_end) {
    return(m_start:m_end)
  } else {
    # Wrap across year end (e.g. Dec-Feb)
    return(c(m_start:12L, 1L:m_end))
  }
}
