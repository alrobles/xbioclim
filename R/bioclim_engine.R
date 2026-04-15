#' Compute Bioclimatic Variables via the Native GDAL-Tiled Engine
#'
#' @description
#' High-level R interface to the \code{BioclimEngine} C++ tiled computation
#' pipeline.  Reads four sets of monthly climate rasters (mean temperature,
#' maximum temperature, minimum temperature, and precipitation) from disk,
#' computes the requested bioclimatic variables, and writes each variable to
#' a separate single-band GeoTIFF inside a directory — one tile at a time so
#' that peak memory is proportional to \code{tile_size}, not the full raster
#' extent.
#'
#' @details
#' **GDAL requirement.** This function requires the package to have been
#' compiled with GDAL support (see \code{\link{has_gdal}}).  If GDAL is not
#' available an informative error is raised immediately.
#'
#' **Variable selection.** By default all 19 standard bioclimatic variables
#' (BIO01–BIO19) are computed and written.
#' Pass \code{variables} as an integer vector (e.g. \code{c(1, 12, 15)}) to
#' restrict the output to a subset.
#' The engine always computes all 19 internally (they share intermediate
#' values), but only the selected ones are written to disk, reducing I/O and
#' file size.
#'
#' **Separate output files.** Each selected variable is written to its own
#' single-band GeoTIFF named \code{bio01.tif} … \code{bio19.tif} inside the
#' \code{output} directory.  This keeps individual files small and allows
#' downstream tools to load only the variables they need.
#'
#' **Tiled processing.** The engine reads and writes rasters in square tiles of
#' \code{tile_size} × \code{tile_size} pixels.  Choosing a large tile improves
#' I/O efficiency; a small tile reduces peak RAM.  The default (256) is a
#' good balance for most use cases.
#'
#' **Multi-band vs. single-band inputs.** Each climate variable can be
#' supplied either as twelve single-band files (one per calendar month) or as
#' one multi-band file with exactly 12 bands.  A \code{terra::SpatRaster} with
#' 12 layers is also accepted; its on-disk source paths are extracted
#' automatically via \code{terra::sources()}.
#'
#' **Mask support.** An optional raster or vector mask can be used to restrict
#' computation to a specific region.  Pixels outside the mask are written as
#' \code{NaN}.  Accepted formats:
#' \itemize{
#'   \item \code{character} — path to any GDAL-readable raster or OGR vector.
#'   \item \code{sf} object — written to a temporary GeoJSON and rasterized.
#'   \item \code{terra::SpatRaster} — written to a temporary GeoTIFF.
#' }
#'
#' **Output.** If \pkg{terra} is installed the function returns a
#' \code{terra::SpatRaster} whose layers correspond to the selected variables.
#' Otherwise it returns a character vector of the output file paths.
#'
#' @param tas Character vector of length 1 (12-band file) or 12 (one file per
#'   month), or a \code{terra::SpatRaster} with 12 layers: monthly mean
#'   temperature.
#' @param tasmax Like \code{tas} but for monthly maximum temperature.
#' @param tasmin Like \code{tas} but for monthly minimum temperature.
#' @param pr Like \code{tas} but for monthly precipitation.
#' @param output Character string: path to the output directory where
#'   individual GeoTIFF files will be written (one per variable, named
#'   \code{bio01.tif} … \code{bio19.tif}).  The directory is created
#'   automatically if it does not exist.  Defaults to a temporary directory.
#' @param variables Integer vector of variable numbers to compute, with
#'   values in \code{1:19}.  Default is \code{1:19} (all 19 variables).
#'   For example, \code{c(1, 12)} writes only BIO01 and BIO12.
#' @param mask Optional mask: a character file path, an \code{sf} object
#'   (polygon), or a \code{terra::SpatRaster}.  \code{NULL} (default) means
#'   no masking.
#' @param threads Positive integer: number of OpenMP threads to use.  Default
#'   is \code{1L}.
#' @param tile_size Positive integer: tile dimension (pixels) for tiled I/O.
#'   Default is \code{256L}.
#' @param overwrite Logical: whether to overwrite existing output files inside
#'   \code{output}.  Default is \code{FALSE}.
#' @param device Character scalar: compute device to use.  One of
#'   \code{"auto"} (default), \code{"cpu"}, or \code{"gpu"}.  \code{"auto"}
#'   selects the GPU when a CUDA device is available, otherwise falls back to
#'   the CPU.  \code{"gpu"} on a system without CUDA emits a warning and falls
#'   back to the CPU.  When \code{device = "gpu"} and \code{tile_size} is left
#'   at its default, the tile size is automatically scaled to match the
#'   detected GPU memory (4096 for high-memory GPUs such as the A100, 1024
#'   otherwise).
#'
#' @return If \pkg{terra} is installed, a \code{terra::SpatRaster} with one
#'   layer per selected variable (named \code{bio01} … \code{bio19}).
#'   Otherwise a character vector of the output file paths.
#'
#' @seealso \code{\link{bioclim_raster}} for the in-memory R/terra path,
#'   \code{\link{has_gdal}} to check GDAL availability,
#'   \code{\link{engine_create}} for the low-level XPtr interface.
#'
#' @export
#' @examples
#' \donttest{
#' if (has_gdal() && requireNamespace("terra", quietly = TRUE)) {
#'   library(terra)
#'
#'   # Create tiny synthetic climate rasters (10x10 pixels, 12 layers each)
#'   make_rast <- function(vals, file) {
#'     r <- rast(nrows = 10, ncols = 10, nlyrs = 12,
#'               xmin = 0, xmax = 1, ymin = 0, ymax = 1, crs = "EPSG:4326")
#'     for (m in seq_len(12)) values(r[[m]]) <- vals[m]
#'     writeRaster(r, file, overwrite = TRUE)
#'     file
#'   }
#'
#'   tmp <- tempdir()
#'   tas_file    <- make_rast(c(5,7,10,14,18,22,25,24,20,15,10,6),
#'                             file.path(tmp, "tas.tif"))
#'   tasmax_file <- make_rast(c(8,10,14,18,23,28,32,31,26,19,13,9),
#'                             file.path(tmp, "tasmax.tif"))
#'   tasmin_file <- make_rast(c(1,3,6,10,13,17,20,19,15,10,6,2),
#'                             file.path(tmp, "tasmin.tif"))
#'   pr_file     <- make_rast(c(60,55,48,35,28,22,18,20,35,55,65,68),
#'                             file.path(tmp, "pr.tif"))
#'
#'   # Compute all 19 variables (one file each in a directory)
#'   out_dir <- file.path(tmp, "bioclim_out")
#'   result <- bioclim_engine(tas_file, tasmax_file, tasmin_file, pr_file,
#'                             output = out_dir, overwrite = TRUE)
#'   nlyr(result)   # 19
#'   list.files(out_dir, pattern = "[.]tif$")
#'
#'   # Compute only BIO01 and BIO12
#'   out_dir2 <- file.path(tmp, "bioclim_subset")
#'   result2 <- bioclim_engine(tas_file, tasmax_file, tasmin_file, pr_file,
#'                              output = out_dir2, variables = c(1L, 12L),
#'                              overwrite = TRUE)
#'   nlyr(result2)  # 2
#'   names(result2) # "bio01" "bio12"
#' }
#' }
bioclim_engine <- function(
    tas,
    tasmax,
    tasmin,
    pr,
    output     = tempfile("bioclim_"),
    variables  = 1:19,
    mask       = NULL,
    threads    = 1L,
    tile_size  = 256L,
    overwrite  = FALSE,
    device     = c("auto", "cpu", "gpu")
) {

  device <- match.arg(device)

  # Auto-detect: use GPU if available, otherwise CPU.
  if (device == "auto") {
    device <- if (has_cuda()) "gpu" else "cpu"
  }

  # GPU requested but not available -> warn and fall back to CPU.
  if (device == "gpu" && !has_cuda()) {
    warning(
      "CUDA GPU requested but not available. Falling back to CPU.",
      call. = FALSE
    )
    device <- "cpu"
  }

  # Auto-scale tile size for GPU when the user did not override it.
  if (device == "gpu" && identical(tile_size, 256L)) {
    info <- cuda_info()
    if (length(info) > 0L && isTRUE(info$memory_gb >= 40)) {
      tile_size <- 4096L
      message(
        "GPU detected (", info$name, ", ",
        round(info$memory_gb, 1), " GB). Using tile_size=4096."
      )
    } else if (length(info) > 0L) {
      tile_size <- 1024L
    }
  }

  # -- 0. GDAL availability ---------------------------------------------------
  if (!has_gdal()) {
    stop(
      "bioclim_engine() requires GDAL support. ",
      "The package was built without GDAL. ",
      "Please reinstall rxbioclim on a system with libgdal-dev.",
      call. = FALSE
    )
  }

  # -- 1. Validate scalar arguments -------------------------------------------

  # output: must be a single character string (directory path)
  if (!is.character(output) || length(output) != 1L || is.na(output)) {
    stop("'output' must be a single non-NA character string.", call. = FALSE)
  }

  # variables: integer vector with values in 1:19
  variables <- suppressWarnings(as.integer(variables))
  if (length(variables) == 0L || any(is.na(variables))) {
    stop("'variables' must be a non-empty integer vector.", call. = FALSE)
  }
  if (any(variables < 1L) || any(variables > 19L)) {
    stop("'variables' must contain values between 1 and 19.", call. = FALSE)
  }
  if (anyDuplicated(variables)) {
    stop("'variables' must not contain duplicates.", call. = FALSE)
  }

  threads_int <- suppressWarnings(as.integer(threads))
  if (length(threads_int) != 1L || is.na(threads_int) || threads_int < 1L) {
    stop("'threads' must be a finite scalar integer >= 1.", call. = FALSE)
  }
  threads <- threads_int

  tile_size_int <- suppressWarnings(as.integer(tile_size))
  if (length(tile_size_int) != 1L || is.na(tile_size_int) ||
      tile_size_int < 1L) {
    stop("'tile_size' must be a finite scalar integer >= 1.", call. = FALSE)
  }
  tile_size <- tile_size_int

  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("'overwrite' must be a single non-NA logical value.", call. = FALSE)
  }

  # -- 2. Prepare output directory --------------------------------------------
  if (file.exists(output) && !file.info(output)$isdir) {
    stop(
      "'output' exists but is not a directory: ", output, "\n",
      "Please supply a directory path.",
      call. = FALSE
    )
  }
  dir.create(output, recursive = TRUE, showWarnings = FALSE)

  # Check for existing output files when overwrite = FALSE
  expected_files <- sprintf("bio%02d.tif", variables)
  existing <- expected_files[file.exists(file.path(output, expected_files))]
  if (length(existing) > 0L && !isTRUE(overwrite)) {
    stop(
      "Output files already exist in '", output, "':\n  ",
      paste(existing, collapse = ", "), "\n",
      "Set overwrite = TRUE to allow overwriting.",
      call. = FALSE
    )
  }

  # -- 3. Resolve climate file paths ------------------------------------------
  tas_files    <- .resolve_climate_input(tas,    "tas")
  tasmax_files <- .resolve_climate_input(tasmax, "tasmax")
  tasmin_files <- .resolve_climate_input(tasmin, "tasmin")
  pr_files     <- .resolve_climate_input(pr,     "pr")

  # -- 4. Resolve mask path (NULL means no mask) ------------------------------
  mask_tmp  <- NULL  # temp file handle for cleanup
  mask_path <- ""    # empty string => no mask in the engine

  if (!is.null(mask)) {
    resolved <- .resolve_engine_mask(mask)
    mask_path <- resolved$path
    mask_tmp  <- resolved$tmp   # non-NULL if we created a temp file
  }

  on.exit({
    if (!is.null(mask_tmp)) unlink(mask_tmp)
  }, add = TRUE)

  # -- 5. Call engine pipeline ------------------------------------------------
  eng <- engine_create()
  engine_open(eng, tas_files, tasmax_files, tasmin_files, pr_files)
  engine_set_output(eng, output)
  engine_set_variables(eng, as.integer(variables))
  if (nzchar(mask_path)) engine_set_mask(eng, mask_path)
  engine_set_threads(eng, threads)
  engine_set_tile_size(eng, tile_size)
  engine_set_device(eng, device)
  engine_compute(eng)

  # -- 6. Return result -------------------------------------------------------
  out_files <- file.path(output, sprintf("bio%02d.tif", variables))

  if (requireNamespace("terra", quietly = TRUE)) {
    layers <- lapply(out_files, terra::rast)
    result <- terra::rast(layers)
    names(result) <- sprintf("bio%02d", variables)
    result
  } else {
    out_files
  }
}

# -- Internal helpers ----------------------------------------------------------

#' Resolve a climate input to a character vector of file paths
#'
#' Accepts a character vector (length 1 or 12) or a terra::SpatRaster (12
#' layers).  Returns a character vector of length 1 or 12 and validates that
#' all files exist.
#'
#' @param x The input to resolve.
#' @param name Variable name for error messages.
#' @return Character vector of length 1 or 12.
#' @keywords internal
.resolve_climate_input <- function(x, name) {

  # terra::SpatRaster -> extract source paths
  if (inherits(x, "SpatRaster")) {
    if (!requireNamespace("terra", quietly = TRUE)) {
      stop(
        "Package 'terra' is required to use a SpatRaster as '", name, "'. ",
        "Install it with: install.packages('terra')",
        call. = FALSE
      )
    }
    if (terra::nlyr(x) != 12L) {
      stop(
        sprintf(
          "'%s' SpatRaster must have 12 layers (one per month), got %d.",
          name, terra::nlyr(x)
        ),
        call. = FALSE
      )
    }
    srcs <- terra::sources(x)
    # terra::sources() returns one path per layer; for a single multi-band
    # file all 12 entries are identical -- collapse to the unique path.
    u <- unique(srcs[nzchar(srcs)])
    if (length(u) == 0L) {
      stop(
        "'", name, "' SpatRaster has no on-disk source. ",
        "Write it to disk first with terra::writeRaster().",
        call. = FALSE
      )
    }
    if (length(u) == 1L) {
      x <- u          # single multi-band file
    } else if (length(u) == 12L) {
      x <- srcs       # 12 single-band files
    } else {
      stop(
        sprintf(
          "'%s' SpatRaster sources resolved to %d unique files; ",
          name, length(u)
        ),
        "expected 1 (multi-band) or 12 (one per month).",
        call. = FALSE
      )
    }
  }

  # character validation
  if (!is.character(x) || !length(x) %in% c(1L, 12L)) {
    stop(
      sprintf(
        "'%s' must be a character vector of length 1 or 12, or a ",
        name
      ),
      "terra::SpatRaster with 12 layers.",
      call. = FALSE
    )
  }

  missing_files <- x[!file.exists(x)]
  if (length(missing_files) > 0L) {
    stop(
      sprintf(
        "The following '%s' files do not exist:\n  %s",
        name, paste(missing_files, collapse = "\n  ")
      ),
      call. = FALSE
    )
  }

  x
}

#' Resolve a mask argument to a file path
#'
#' Accepts a character string, an sf object, or a terra::SpatRaster.
#' Returns a list with \code{path} (character) and \code{tmp} (path to clean
#' up, or \code{NULL}).
#'
#' @param mask The mask argument passed by the user.
#' @return Named list with \code{path} and \code{tmp}.
#' @keywords internal
.resolve_engine_mask <- function(mask) {

  # terra::SpatRaster -> write to temp GeoTIFF
  if (inherits(mask, "SpatRaster")) {
    if (!requireNamespace("terra", quietly = TRUE)) {
      stop(
        "Package 'terra' is required to use a SpatRaster as 'mask'. ",
        "Install it with: install.packages('terra')",
        call. = FALSE
      )
    }
    tmp <- tempfile(fileext = ".tif")
    terra::writeRaster(mask, tmp, overwrite = TRUE)
    return(list(path = tmp, tmp = tmp))
  }

  # sf / sfc -> write to temp GeoJSON
  if (inherits(mask, "sf") || inherits(mask, "sfc")) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop(
        "Package 'sf' is required to use an sf object as 'mask'. ",
        "Install it with: install.packages('sf')",
        call. = FALSE
      )
    }
    tmp <- tempfile(fileext = ".geojson")
    sf::st_write(mask, tmp, quiet = TRUE)
    return(list(path = tmp, tmp = tmp))
  }

  # character file path
  if (!is.character(mask) || length(mask) != 1L || is.na(mask)) {
    stop(
      "'mask' must be a character file path, an sf object, or a ",
      "terra::SpatRaster.",
      call. = FALSE
    )
  }
  if (!file.exists(mask)) {
    stop("'mask' file does not exist: ", mask, call. = FALSE)
  }

  list(path = mask, tmp = NULL)
}
