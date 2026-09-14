#' Compute Bioclimatic Variables from ERA5-Land Data
#'
#' @description
#' End-to-end pipeline that reads ERA5-Land hourly GRIB/NetCDF files,
#' aggregates them to CHELSA-compatible monthly climate variables, and computes
#' the 19 standard bioclimatic variables (BIO01–BIO19).
#'
#' @details
#' The pipeline proceeds in three stages:
#' \enumerate{
#'   \item \strong{Monthly aggregation}: For each of the 12 calendar months,
#'     hourly \code{t2m} is aggregated to \code{tas}, \code{tasmax}, and
#'     \code{tasmin}; hourly \code{tp} is summed to \code{pr}.
#'   \item \strong{Stack}: The 12 monthly layers are assembled into
#'     \code{terra::SpatRaster} objects with 12 bands each.
#'   \item \strong{Bioclim}: \code{\link{bioclim_raster}} computes BIO01–BIO19.
#' }
#'
#' @param t2m_files Character vector of 12 file paths to monthly ERA5-Land
#'   hourly 2-m temperature files (one per calendar month, January–December).
#'   Each file may be GRIB or NetCDF.
#' @param tp_files Character vector of 12 file paths to monthly ERA5-Land
#'   hourly total precipitation files (same order as \code{t2m_files}).
#' @param year Integer: the calendar year (used to determine days per month).
#' @param output Character path to an output directory for GeoTIFF files.
#'   Defaults to a temporary directory.
#' @param to_celsius Logical: convert temperatures to Celsius?
#'   Default \code{TRUE} for WorldClim-convention bioclimatic variables.
#' @param variables Integer vector of bioclimatic variables to compute
#'   (1–19).
#'   Default \code{1:19} (all).
#' @param ncores Integer: OpenMP threads for aggregation. Default \code{1L}.
#' @param save_monthly Logical: write intermediate monthly GeoTIFFs?
#'   Default \code{FALSE}.
#'
#' @return A \code{terra::SpatRaster} with one layer per bioclimatic variable.
#'
#' @seealso \code{\link{era5_t2m_to_monthly}}, \code{\link{era5_tp_to_monthly}},
#'   \code{\link{bioclim_raster}}
#'
#' @examples
#' \dontrun{
#' # Paths to ERA5-Land GRIB files on the HPC cluster
#' t2m_files <- sprintf("era5land_t2m_hourly_2020_%02d.grib", 1:12)
#' tp_files  <- sprintf("era5land_tp_hourly_2020_%02d.grib", 1:12)
#'
#' bio <- era5_bioclim(t2m_files, tp_files, year = 2020L, ncores = 4L)
#' terra::plot(bio[[1]])  # BIO01
#' }
#'
#' @export
era5_bioclim <- function(t2m_files,
                         tp_files,
                         year,
                         output = tempdir(),
                         to_celsius = TRUE,
                         variables = 1:19,
                         ncores = 1L,
                         save_monthly = FALSE) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required for era5_bioclim(). ",
         "Install it with: install.packages('terra')", call. = FALSE)
  }

  stopifnot(length(t2m_files) == 12L)
  stopifnot(length(tp_files)  == 12L)

  if (!dir.exists(output)) dir.create(output, recursive = TRUE)

  # Days per month for the given year
  n_days <- vapply(1:12, function(m) days_in_month(year, m), integer(1))

  # Allocate lists for monthly rasters
  tas_list    <- vector("list", 12L)
  tasmax_list <- vector("list", 12L)
  tasmin_list <- vector("list", 12L)
  pr_list     <- vector("list", 12L)

  for (m in seq_len(12L)) {
    message(sprintf("[xbioclim] Processing month %02d/%d (%d days)...",
                    m, year, n_days[m]))

    # Read hourly t2m
    r_t2m <- terra::rast(t2m_files[m])
    expected_hours <- 24L * n_days[m]
    if (terra::nlyr(r_t2m) != expected_hours) {
      warning(sprintf(
        "Month %d: expected %d layers (24 x %d days), got %d",
        m, expected_hours, n_days[m], terra::nlyr(r_t2m)
      ))
    }

    # Read hourly tp
    r_tp <- terra::rast(tp_files[m])

    # Extract values as matrices [n_pixels x n_hours]
    vals_t2m <- terra::values(r_t2m)
    vals_tp  <- terra::values(r_tp)

    # Aggregate using C++ backend
    temp_result <- era5_t2m_to_monthly_cpp(
      vals_t2m, n_days[m], to_celsius, ncores
    )
    pr_vals <- era5_tp_to_monthly_cpp(vals_tp, ncores)

    # Create single-band rasters from results
    template <- r_t2m[[1]]
    tas_list[[m]]    <- terra::setValues(terra::rast(template), temp_result$tas)
    tasmax_list[[m]] <- terra::setValues(terra::rast(template), temp_result$tasmax)
    tasmin_list[[m]] <- terra::setValues(terra::rast(template), temp_result$tasmin)
    pr_list[[m]]     <- terra::setValues(terra::rast(template), pr_vals)

    if (save_monthly) {
      terra::writeRaster(
        tas_list[[m]],
        file.path(output, sprintf("tas_%d_%02d.tif", year, m)),
        overwrite = TRUE
      )
      terra::writeRaster(
        tasmax_list[[m]],
        file.path(output, sprintf("tasmax_%d_%02d.tif", year, m)),
        overwrite = TRUE
      )
      terra::writeRaster(
        tasmin_list[[m]],
        file.path(output, sprintf("tasmin_%d_%02d.tif", year, m)),
        overwrite = TRUE
      )
      terra::writeRaster(
        pr_list[[m]],
        file.path(output, sprintf("pr_%d_%02d.tif", year, m)),
        overwrite = TRUE
      )
    }
  }

  # Stack monthly layers into 12-band rasters
  tas_stack    <- terra::rast(tas_list)
  tasmax_stack <- terra::rast(tasmax_list)
  tasmin_stack <- terra::rast(tasmin_list)
  pr_stack     <- terra::rast(pr_list)

  names(tas_stack)    <- sprintf("tas_%02d", 1:12)
  names(tasmax_stack) <- sprintf("tasmax_%02d", 1:12)
  names(tasmin_stack) <- sprintf("tasmin_%02d", 1:12)
  names(pr_stack)     <- sprintf("pr_%02d", 1:12)

  message("[xbioclim] Computing bioclimatic variables...")

  bio <- bioclim_raster(tas_stack, tasmax_stack, tasmin_stack, pr_stack,
                        ncores = ncores)
  if (!identical(variables, 1:19)) {
    bio <- bio[[variables]]
  }

  message("[xbioclim] Done.")
  bio
}

#' Generate Bioclimatic Variables for Multiple Years
#'
#' Batch-processes multiple years of ERA5-Land data through the
#' \code{\link{era5_bioclim}} pipeline.
#'
#' @param base_dir Character: base directory containing ERA5-Land raw data.
#'   Expected structure: \code{base_dir/t2m/YYYY/era5land_t2m_hourly_YYYY_MM.grib}
#'   and \code{base_dir/tp/YYYY/era5land_tp_hourly_YYYY_MM.grib}.
#' @param years Integer vector of years to process.
#' @param output_dir Character: output directory.
#' @param t2m_pattern Character: filename pattern for t2m files.
#'   Must contain \code{\%d} for year and \code{\%02d} for month.
#'   Default: \code{"era5land_t2m_hourly_\%d_\%02d.grib"}.
#' @param tp_pattern Character: filename pattern for tp files.
#'   Default: \code{"era5land_tp_hourly_\%d_\%02d.grib"}.
#' @param ... Additional arguments passed to \code{\link{era5_bioclim}}.
#'
#' @return A named list of \code{terra::SpatRaster} objects, one per year.
#'
#' @examples
#' \dontrun{
#' bio_all <- era5_bioclim_years(
#'   base_dir   = "/scratch/era5-land/raw",
#'   years      = 1980:2020,
#'   output_dir = "/scratch/era5-land/bioclim",
#'   ncores     = 8L
#' )
#' }
#'
#' @export
era5_bioclim_years <- function(base_dir,
                               years,
                               output_dir,
                               t2m_pattern = "era5land_t2m_hourly_%d_%02d.grib",
                               tp_pattern  = "era5land_tp_hourly_%d_%02d.grib",
                               ...) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  results <- vector("list", length(years))
  names(results) <- as.character(years)

  for (i in seq_along(years)) {
    y <- years[i]
    message(sprintf("[xbioclim] === Year %d (%d/%d) ===", y, i, length(years)))

    t2m_files <- file.path(
      base_dir, "t2m", y,
      sprintf(t2m_pattern, y, 1:12)
    )
    tp_files <- file.path(
      base_dir, "tp", y,
      sprintf(tp_pattern, y, 1:12)
    )

    # Verify files exist
    missing_t2m <- !file.exists(t2m_files)
    missing_tp  <- !file.exists(tp_files)
    if (any(missing_t2m) || any(missing_tp)) {
      warning(sprintf(
        "Year %d: missing %d t2m and %d tp files, skipping.",
        y, sum(missing_t2m), sum(missing_tp)
      ))
      next
    }

    year_output <- file.path(output_dir, as.character(y))
    results[[as.character(y)]] <- era5_bioclim(
      t2m_files  = t2m_files,
      tp_files   = tp_files,
      year       = y,
      output     = year_output,
      ...
    )
  }

  results
}
