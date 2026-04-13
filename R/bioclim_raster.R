#' Block-Based Raster Processing for Bioclimatic Variables
#'
#' @description
#' Functions for computing bioclimatic variables from monthly climate rasters
#' (SpatRaster objects) using terra's block-loop architecture for
#' memory-efficient processing of large rasters.
#'
#' @name bioclim-raster
NULL

#' Validate SpatRaster Input
#'
#' Checks that input is a SpatRaster with 12 layers.
#'
#' @param x The input to validate.
#' @param name Name of the variable for error messages.
#' @return Invisible NULL. Throws an error if validation fails.
#' @keywords internal
validate_spatraster <- function(x, name = "input") {
  if (!inherits(x, "SpatRaster")) {
    stop(sprintf("'%s' must be a SpatRaster", name), call. = FALSE)
  }
  if (terra::nlyr(x) != 12L) {
    stop(
      sprintf("'%s' must have 12 layers (one per month), got %d",
              name, terra::nlyr(x)),
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Compute Bioclimatic Variables for a Block of Cells
#'
#' Processes a matrix of monthly climate values (one row per cell, 12 columns
#' per month) and returns a matrix of 19 bioclimatic variable values.
#' Cells with any NA input values are returned as all-NA rows.
#'
#' @param v_tas    Numeric matrix (n_cells x 12): monthly mean temperature.
#' @param v_tasmax Numeric matrix (n_cells x 12): monthly max temperature.
#' @param v_tasmin Numeric matrix (n_cells x 12): monthly min temperature.
#' @param v_pr     Numeric matrix (n_cells x 12): monthly precipitation.
#'
#' @return A numeric matrix (n_cells x 19) of bioclimatic variable values.
#' @keywords internal
bioclim_block <- function(v_tas, v_tasmax, v_tasmin, v_pr) {
  n_cells <- nrow(v_tas)
  result  <- matrix(NA_real_, nrow = n_cells, ncol = 19L)
  for (j in seq_len(n_cells)) {
    t_row  <- v_tas[j, ]
    tx_row <- v_tasmax[j, ]
    tn_row <- v_tasmin[j, ]
    p_row  <- v_pr[j, ]
    if (anyNA(t_row) || anyNA(tx_row) || anyNA(tn_row) || anyNA(p_row)) {
      next
    }
    result[j, ] <- bioclim(t_row, tx_row, tn_row, p_row)
  }
  result
}

#' Compute Bioclimatic Variables from Monthly Climate Rasters
#'
#' @description
#' Computes the 19 standard bioclimatic variables (BIO01-BIO19) from monthly
#' climate raster data (SpatRaster objects with 12 layers, one per month) using
#' terra's block-loop architecture for memory-efficient processing.
#'
#' @details
#' This function follows terra's block-loop pattern: the input rasters are read
#' one horizontal block at a time (using `terra::readStart()`,
#' `terra::readValues()`, and `terra::readStop()`), keeping peak memory use
#' proportional to the block size rather than the full raster extent. Results
#' are written to the output raster block by block (using
#' `terra::writeStart()`, `terra::writeValues()`, and `terra::writeStop()`).
#'
#' When `ncores > 1`, cells within each block are split across workers using
#' `parallel::makeCluster()` / `parallel::parLapply()`, enabling multi-core
#' computation. Only one block is held in memory at a time regardless of the
#' number of cores.
#'
#' @param tas     SpatRaster with 12 layers: monthly mean temperature.
#' @param tasmax  SpatRaster with 12 layers: monthly maximum temperature.
#' @param tasmin  SpatRaster with 12 layers: monthly minimum temperature.
#' @param pr      SpatRaster with 12 layers: monthly precipitation.
#' @param filename Character string: output file path. Pass `""` (default) to
#'   keep the result in memory.
#' @param n_blocks Integer: target number of row blocks. If `NULL` (default),
#'   terra selects an appropriate number based on available memory.
#' @param ncores  Integer: number of CPU cores for within-block parallel
#'   processing. Default is `1` (sequential). Values greater than 1 require
#'   the \pkg{parallel} package.
#' @param overwrite Logical: whether to overwrite an existing output file.
#'   Default is `FALSE`.
#' @param ...     Additional arguments passed to `terra::writeStart()`.
#'
#' @return A \link[terra]{SpatRaster} with 19 layers named `bio01` through
#'   `bio19`, sharing the spatial extent, resolution and CRS of `tas`.
#'
#' @seealso [bioclim()] for single-pixel (vector) computation.
#'
#' @export
#' @examplesIf requireNamespace("terra", quietly = TRUE)
#' library(terra)
#'
#' # Create small test rasters: 4 rows x 3 cols, 12 monthly layers
#' make_rast <- function(vals) {
#'   r <- rast(nrows = 4, ncols = 3, nlyr = 12)
#'   values(r) <- vals
#'   r
#' }
#'
#' n <- 4 * 3  # number of cells
#' tas    <- make_rast(matrix(rep(1:12, each = n), nrow = n, ncol = 12))
#' tasmax <- make_rast(matrix(rep(2:13, each = n), nrow = n, ncol = 12))
#' tasmin <- make_rast(matrix(rep(0:11, each = n), nrow = n, ncol = 12))
#' pr     <- make_rast(matrix(rep(1:12, each = n), nrow = n, ncol = 12))
#'
#' # Compute all 19 bioclimatic variables in a single pass
#' result <- bioclim_raster(tas, tasmax, tasmin, pr)
#' nlyr(result)  # 19
bioclim_raster <- function(tas, tasmax, tasmin, pr,
                            filename = "", n_blocks = NULL,
                            ncores = 1L, overwrite = FALSE, ...) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop(
      "Package 'terra' is required for bioclim_raster(). ",
      "Install it with: install.packages('terra')",
      call. = FALSE
    )
  }

  validate_spatraster(tas,    "tas")
  validate_spatraster(tasmax, "tasmax")
  validate_spatraster(tasmin, "tasmin")
  validate_spatraster(pr,     "pr")

  ncores_int <- suppressWarnings(as.integer(ncores))
  if (length(ncores_int) != 1L || is.na(ncores_int) || ncores_int < 1L) {
    stop("'ncores' must be a finite scalar integer >= 1", call. = FALSE)
  }
  ncores <- ncores_int

  # Create output raster: 19 layers, same footprint as tas
  out <- terra::rast(tas[[1L]], nlyrs = 19L)
  names(out) <- paste0("bio", sprintf("%02d", 1:19))

  # Determine blocks
  bk <- if (is.null(n_blocks)) {
    terra::blocks(out)
  } else {
    terra::blocks(out, n = as.integer(n_blocks))
  }

  # State variables for cleanup tracking
  cl            <- NULL
  write_started <- FALSE

  # Single comprehensive on.exit: runs on both normal exit and errors.
  # try() prevents cascading failures during cleanup.
  on.exit(
    {
      try(terra::readStop(tas),    silent = TRUE)
      try(terra::readStop(tasmax), silent = TRUE)
      try(terra::readStop(tasmin), silent = TRUE)
      try(terra::readStop(pr),     silent = TRUE)
      if (write_started) try(terra::writeStop(out), silent = TRUE)
      if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE)
    },
    add = TRUE
  )

  # Setup parallel cluster if requested
  if (ncores > 1L) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      warning(
        "Package 'parallel' is not available; falling back to ncores = 1.",
        call. = FALSE
      )
      ncores <- 1L
    } else {
      cl <- parallel::makeCluster(ncores)
      # Export all package functions needed by workers, sourcing explicitly
      # from the package namespace so closures resolve correctly on workers.
      parallel::clusterExport(
        cl,
        varlist = c(
          "bioclim", "bioclim_block", "sd_pop",
          "rolling_quarter_sum", "rolling_quarter_mean",
          "quarter_argmax", "quarter_argmin",
          "quarter_values", "validate_monthly"
        ),
        envir = asNamespace("rxbioclim")
      )
    }
  }

  # Open input rasters for block-reading
  terra::readStart(tas)
  terra::readStart(tasmax)
  terra::readStart(tasmin)
  terra::readStart(pr)

  # Open output raster for block-writing
  terra::writeStart(out, filename = filename, overwrite = overwrite, ...)
  write_started <- TRUE

  # Block loop: read -> compute -> write
  for (i in seq_len(bk$n)) {
    row_start <- bk$row[i]
    n_rows    <- bk$nrows[i]

    v_tas    <- terra::readValues(tas,    row = row_start,
                                  nrows = n_rows, mat = TRUE)
    v_tasmax <- terra::readValues(tasmax, row = row_start,
                                  nrows = n_rows, mat = TRUE)
    v_tasmin <- terra::readValues(tasmin, row = row_start,
                                  nrows = n_rows, mat = TRUE)
    v_pr     <- terra::readValues(pr,     row = row_start,
                                  nrows = n_rows, mat = TRUE)

    n_cells <- nrow(v_tas)

    if (!is.null(cl)) {
      # Parallel: split cells into chunks across workers.
      # `bioclim_block` is available on workers via clusterExport above.
      chunks       <- parallel::splitIndices(n_cells, ncores)
      result_parts <- parallel::parLapply(cl, chunks, function(idx) {
        bioclim_block(
          v_tas[idx,    , drop = FALSE],
          v_tasmax[idx, , drop = FALSE],
          v_tasmin[idx, , drop = FALSE],
          v_pr[idx,     , drop = FALSE]
        )
      })
      result_mat <- do.call(rbind, result_parts)
    } else {
      result_mat <- bioclim_block(v_tas, v_tasmax, v_tasmin, v_pr)
    }

    terra::writeValues(out, result_mat, start = row_start, nrows = n_rows)
  }

  # Explicit finalization on the happy path; prevents double-call from on.exit
  write_started <- FALSE
  terra::readStop(tas)
  terra::readStop(tasmax)
  terra::readStop(tasmin)
  terra::readStop(pr)
  terra::writeStop(out)
  if (!is.null(cl)) {
    parallel::stopCluster(cl)
    cl <- NULL
  }
  out
}
