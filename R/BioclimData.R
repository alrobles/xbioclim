#' S4 Class for Bioclimatic Input Data
#'
#' @description
#' Holds monthly climate data for one or more pixels (a raster block) and
#' exposes S4 methods for computing each of the 19 standard bioclimatic
#' variables as well as the full batch computation.  Each slot is a numeric
#' matrix with **12 columns** (one per calendar month) and **one row per
#' pixel**; single-pixel inputs (plain numeric vectors of length 12) are
#' automatically promoted to 1-row matrices by the constructor.
#'
#' The S4 methods delegate to the compiled C++ routines exported by the
#' package (e.g. `bio01_cpp`, `bioclim_cpp`), mirroring the xbioclim
#' convention for all 19 variables.
#'
#' @slot tas    Numeric matrix (pixels × 12): monthly mean temperature.
#' @slot tasmax Numeric matrix (pixels × 12): monthly maximum temperature.
#' @slot tasmin Numeric matrix (pixels × 12): monthly minimum temperature.
#' @slot pr     Numeric matrix (pixels × 12): monthly precipitation.
#'
#' @export
setClass(
  "BioclimData",
  representation(
    tas    = "matrix",
    tasmax = "matrix",
    tasmin = "matrix",
    pr     = "matrix"
  ),
  validity = function(object) {
    msgs <- character(0)

    check_slot <- function(x, nm) {
      if (!is.numeric(x))
        msgs <<- c(msgs, sprintf("'%s' must be numeric", nm))
      if (ncol(x) != 12L)
        msgs <<- c(msgs, sprintf("'%s' must have 12 columns (one per month), got %d",
                                 nm, ncol(x)))
    }

    check_slot(object@tas,    "tas")
    check_slot(object@tasmax, "tasmax")
    check_slot(object@tasmin, "tasmin")
    check_slot(object@pr,     "pr")

    n_rows <- nrow(object@tas)
    for (nm in c("tasmax", "tasmin", "pr")) {
      sl <- slot(object, nm)
      if (nrow(sl) != n_rows)
        msgs <- c(msgs, sprintf("'%s' must have %d rows (same as 'tas'), got %d",
                                nm, n_rows, nrow(sl)))
    }

    if (length(msgs)) msgs else TRUE
  }
)

# ── Constructor ────────────────────────────────────────────────────────────────

#' Create a BioclimData Object
#'
#' Constructs a [BioclimData-class] object from monthly climate data.
#' Plain numeric vectors of length 12 are automatically coerced to 1-row
#' matrices so that single-pixel and raster-block inputs are handled uniformly.
#'
#' @param tas    Numeric vector (length 12) **or** matrix (pixels × 12):
#'   monthly mean temperature.
#' @param tasmax Numeric vector (length 12) **or** matrix (pixels × 12):
#'   monthly maximum temperature.
#' @param tasmin Numeric vector (length 12) **or** matrix (pixels × 12):
#'   monthly minimum temperature.
#' @param pr     Numeric vector (length 12) **or** matrix (pixels × 12):
#'   monthly precipitation.
#'
#' @return A [BioclimData-class] object.
#'
#' @export
#' @examples
#' tas    <- 1:12
#' tasmax <- 2:13
#' tasmin <- 0:11
#' pr     <- 1:12
#' bd <- BioclimData(tas, tasmax, tasmin, pr)
#' bd
BioclimData <- function(tas, tasmax, tasmin, pr) {
  to_matrix <- function(x, nm) {
    if (is.vector(x)) {
      if (!is.numeric(x))
        stop(sprintf("'%s' must be numeric", nm), call. = FALSE)
      if (length(x) != 12L)
        stop(sprintf("'%s' must have length 12 (one value per month), got %d",
                     nm, length(x)), call. = FALSE)
      matrix(x, nrow = 1L)
    } else if (is.matrix(x)) {
      if (!is.numeric(x))
        stop(sprintf("'%s' must be numeric", nm), call. = FALSE)
      storage.mode(x) <- "double"
      x
    } else {
      stop(sprintf("'%s' must be a numeric vector or matrix", nm), call. = FALSE)
    }
  }
  new("BioclimData",
      tas    = to_matrix(tas,    "tas"),
      tasmax = to_matrix(tasmax, "tasmax"),
      tasmin = to_matrix(tasmin, "tasmin"),
      pr     = to_matrix(pr,     "pr"))
}

# ── show ───────────────────────────────────────────────────────────────────────

setMethod("show", "BioclimData", function(object) {
  n <- nrow(object@tas)
  cat(sprintf("BioclimData: %d pixel%s, 12 months\n", n, if (n == 1) "" else "s"))
  cat("  Slots: tas, tasmax, tasmin, pr  (each a ", n, " x 12 matrix)\n", sep = "")
  invisible(object)
})

# ── S4 generics ────────────────────────────────────────────────────────────────
# We promote the existing S4-unaware R functions to proper S4 generics,
# preserving backwards compatibility for plain numeric-vector callers while
# adding efficient C++-backed methods for BioclimData objects.

#' @rdname bioclim-variables
#' @export
setGeneric("bio01",   function(tas, ...)                   standardGeneric("bio01"))
#' @rdname bioclim-variables
#' @export
setGeneric("bio02",   function(tasmax, tasmin, ...)         standardGeneric("bio02"))
#' @rdname bioclim-variables
#' @export
setGeneric("bio03",   function(tasmax, tasmin, ...)         standardGeneric("bio03"))
#' @rdname bioclim-variables
#' @export
setGeneric("bio04",   function(tas, ...)                   standardGeneric("bio04"))
#' @rdname bioclim-variables
#' @export
setGeneric("bio05",   function(tasmax, ...)                standardGeneric("bio05"))
#' @rdname bioclim-variables
#' @export
setGeneric("bio06",   function(tasmin, ...)                standardGeneric("bio06"))
#' @rdname bioclim-variables
#' @export
setGeneric("bio07",   function(tasmax, tasmin, ...)         standardGeneric("bio07"))
#' @rdname bioclim-variables
#' @export
setGeneric("bio08",   function(tas, pr, ...)               standardGeneric("bio08"))
#' @rdname bioclim-variables
#' @export
setGeneric("bio09",   function(tas, pr, ...)               standardGeneric("bio09"))
#' @rdname bioclim-variables
#' @export
setGeneric("bio10",   function(tas, ...)                   standardGeneric("bio10"))
#' @rdname bioclim-variables
#' @export
setGeneric("bio11",   function(tas, ...)                   standardGeneric("bio11"))
#' @rdname bioclim-variables
#' @export
setGeneric("bio12",   function(pr, ...)                    standardGeneric("bio12"))
#' @rdname bioclim-variables
#' @export
setGeneric("bio13",   function(pr, ...)                    standardGeneric("bio13"))
#' @rdname bioclim-variables
#' @export
setGeneric("bio14",   function(pr, ...)                    standardGeneric("bio14"))
#' @rdname bioclim-variables
#' @export
setGeneric("bio15",   function(pr, ...)                    standardGeneric("bio15"))
#' @rdname bioclim-variables
#' @export
setGeneric("bio16",   function(pr, ...)                    standardGeneric("bio16"))
#' @rdname bioclim-variables
#' @export
setGeneric("bio17",   function(pr, ...)                    standardGeneric("bio17"))
#' @rdname bioclim-variables
#' @export
setGeneric("bio18",   function(tas, pr, ...)               standardGeneric("bio18"))
#' @rdname bioclim-variables
#' @export
setGeneric("bio19",   function(tas, pr, ...)               standardGeneric("bio19"))
#' @rdname bioclim-variables
#' @export
setGeneric("bioclim", function(tas, tasmax, tasmin, pr, ...)
  standardGeneric("bioclim"))

# ── numeric/ANY methods (existing single-pixel R behaviour) ───────────────────

setMethod("bio01",   "ANY",
  function(tas, ...) { validate_monthly(tas, "tas"); mean(tas) })

setMethod("bio02",   signature("ANY", "ANY"),
  function(tasmax, tasmin, ...) {
    validate_monthly(tasmax, "tasmax"); validate_monthly(tasmin, "tasmin")
    mean(tasmax - tasmin) })

setMethod("bio03",   signature("ANY", "ANY"),
  function(tasmax, tasmin, ...) {
    b02 <- bio02(tasmax, tasmin)
    b07 <- bio07(tasmax, tasmin)
    if (b07 == 0) return(NaN)
    100 * b02 / b07 })

setMethod("bio04",   "ANY",
  function(tas, ...) { validate_monthly(tas, "tas"); 100 * sd_pop(tas) })

setMethod("bio05",   "ANY",
  function(tasmax, ...) { validate_monthly(tasmax, "tasmax"); max(tasmax) })

setMethod("bio06",   "ANY",
  function(tasmin, ...) { validate_monthly(tasmin, "tasmin"); min(tasmin) })

setMethod("bio07",   signature("ANY", "ANY"),
  function(tasmax, tasmin, ...) { bio05(tasmax) - bio06(tasmin) })

setMethod("bio08",   signature("ANY", "ANY"),
  function(tas, pr, ...) {
    validate_monthly(tas, "tas"); validate_monthly(pr, "pr")
    mean(quarter_values(tas, quarter_argmax(pr))) })

setMethod("bio09",   signature("ANY", "ANY"),
  function(tas, pr, ...) {
    validate_monthly(tas, "tas"); validate_monthly(pr, "pr")
    mean(quarter_values(tas, quarter_argmin(pr))) })

setMethod("bio10",   "ANY",
  function(tas, ...) {
    validate_monthly(tas, "tas")
    mean(quarter_values(tas, quarter_argmax(tas))) })

setMethod("bio11",   "ANY",
  function(tas, ...) {
    validate_monthly(tas, "tas")
    mean(quarter_values(tas, quarter_argmin(tas))) })

setMethod("bio12",   "ANY",
  function(pr, ...) { validate_monthly(pr, "pr"); sum(pr) })

setMethod("bio13",   "ANY",
  function(pr, ...) { validate_monthly(pr, "pr"); max(pr) })

setMethod("bio14",   "ANY",
  function(pr, ...) { validate_monthly(pr, "pr"); min(pr) })

setMethod("bio15",   "ANY",
  function(pr, ...) {
    validate_monthly(pr, "pr")
    pr_mean <- mean(pr)
    if (pr_mean == 0) return(NaN)
    100 * sd_pop(pr) / pr_mean })

setMethod("bio16",   "ANY",
  function(pr, ...) {
    validate_monthly(pr, "pr")
    sum(quarter_values(pr, quarter_argmax(pr))) })

setMethod("bio17",   "ANY",
  function(pr, ...) {
    validate_monthly(pr, "pr")
    sum(quarter_values(pr, quarter_argmin(pr))) })

setMethod("bio18",   signature("ANY", "ANY"),
  function(tas, pr, ...) {
    validate_monthly(tas, "tas"); validate_monthly(pr, "pr")
    sum(quarter_values(pr, quarter_argmax(tas))) })

setMethod("bio19",   signature("ANY", "ANY"),
  function(tas, pr, ...) {
    validate_monthly(tas, "tas"); validate_monthly(pr, "pr")
    sum(quarter_values(pr, quarter_argmin(tas))) })

setMethod("bioclim", signature("ANY", "ANY", "ANY", "ANY"),
  function(tas, tasmax, tasmin, pr, ...) {
    validate_monthly(tas, "tas")
    validate_monthly(tasmax, "tasmax")
    validate_monthly(tasmin, "tasmin")
    validate_monthly(pr, "pr")

    b01 <- mean(tas)
    b02 <- mean(tasmax - tasmin)
    b05 <- max(tasmax)
    b06 <- min(tasmin)
    b07 <- b05 - b06
    b03 <- if (b07 == 0) NaN else 100 * b02 / b07
    b04 <- 100 * sd_pop(tas)

    wet_start  <- quarter_argmax(pr)
    dry_start  <- quarter_argmin(pr)
    warm_start <- quarter_argmax(tas)
    cold_start <- quarter_argmin(tas)

    b08 <- mean(quarter_values(tas, wet_start))
    b09 <- mean(quarter_values(tas, dry_start))
    b10 <- mean(quarter_values(tas, warm_start))
    b11 <- mean(quarter_values(tas, cold_start))

    b12 <- sum(pr)
    b13 <- max(pr)
    b14 <- min(pr)
    pr_mean <- mean(pr)
    b15 <- if (pr_mean == 0) NaN else 100 * sd_pop(pr) / pr_mean

    b16 <- sum(quarter_values(pr, wet_start))
    b17 <- sum(quarter_values(pr, dry_start))
    b18 <- sum(quarter_values(pr, warm_start))
    b19 <- sum(quarter_values(pr, cold_start))

    c(bio01 = b01, bio02 = b02, bio03 = b03, bio04 = b04,
      bio05 = b05, bio06 = b06, bio07 = b07,
      bio08 = b08, bio09 = b09, bio10 = b10, bio11 = b11,
      bio12 = b12, bio13 = b13, bio14 = b14, bio15 = b15,
      bio16 = b16, bio17 = b17, bio18 = b18, bio19 = b19)
  })

# ── BioclimData methods (C++ raster-block paths) ──────────────────────────────

setMethod("bio01",   "BioclimData",
  function(tas, ...) bio01_cpp(tas@tas))

setMethod("bio02",   signature("BioclimData", "missing"),
  function(tasmax, tasmin, ...) bio02_cpp(tasmax@tasmax, tasmax@tasmin))

setMethod("bio03",   signature("BioclimData", "missing"),
  function(tasmax, tasmin, ...) bio03_cpp(tasmax@tasmax, tasmax@tasmin))

setMethod("bio04",   "BioclimData",
  function(tas, ...) bio04_cpp(tas@tas))

setMethod("bio05",   "BioclimData",
  function(tasmax, ...) bio05_cpp(tasmax@tasmax))

setMethod("bio06",   "BioclimData",
  function(tasmin, ...) bio06_cpp(tasmin@tasmin))

setMethod("bio07",   signature("BioclimData", "missing"),
  function(tasmax, tasmin, ...) bio07_cpp(tasmax@tasmax, tasmax@tasmin))

setMethod("bio08",   signature("BioclimData", "missing"),
  function(tas, pr, ...) bio08_cpp(tas@tas, tas@pr))

setMethod("bio09",   signature("BioclimData", "missing"),
  function(tas, pr, ...) bio09_cpp(tas@tas, tas@pr))

setMethod("bio10",   "BioclimData",
  function(tas, ...) bio10_cpp(tas@tas))

setMethod("bio11",   "BioclimData",
  function(tas, ...) bio11_cpp(tas@tas))

setMethod("bio12",   "BioclimData",
  function(pr, ...) bio12_cpp(pr@pr))

setMethod("bio13",   "BioclimData",
  function(pr, ...) bio13_cpp(pr@pr))

setMethod("bio14",   "BioclimData",
  function(pr, ...) bio14_cpp(pr@pr))

setMethod("bio15",   "BioclimData",
  function(pr, ...) bio15_cpp(pr@pr))

setMethod("bio16",   "BioclimData",
  function(pr, ...) bio16_cpp(pr@pr))

setMethod("bio17",   "BioclimData",
  function(pr, ...) bio17_cpp(pr@pr))

setMethod("bio18",   signature("BioclimData", "missing"),
  function(tas, pr, ...) bio18_cpp(tas@tas, tas@pr))

setMethod("bio19",   signature("BioclimData", "missing"),
  function(tas, pr, ...) bio19_cpp(tas@tas, tas@pr))

setMethod("bioclim", signature("BioclimData", "missing", "missing", "missing"),
  function(tas, tasmax, tasmin, pr, ...) {
    bioclim_cpp(tas@tas, tas@tasmax, tas@tasmin, tas@pr)
  })
