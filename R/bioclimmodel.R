#' BioclimModel S4 Class
#'
#' An S4 class that wraps a C++ \code{BioclimModel} object via an opaque
#' external pointer handle, following the terra package pattern for C++ object
#' handles. All 19 bioclimatic variable computations are delegated to the
#' underlying C++ object via Rcpp.
#'
#' @slot pntr An \code{externalptr} to the underlying C++ \code{BioclimModel}
#'   object.
#'
#' @name BioclimModel-class
#' @exportClass BioclimModel
setClass("BioclimModel", representation(pntr = "externalptr"))

# Validity function: ensures the slot holds a live, non-null C++ pointer.
setValidity("BioclimModel", function(object) {
  if (is.null(object@pntr) || bioclim_model_is_null(object@pntr)) {
    "BioclimModel contains a null or invalid C++ pointer"
  } else {
    TRUE
  }
})

# ── Constructor ──────────────────────────────────────────────────────────────

#' Create a BioclimModel Object
#'
#' Constructs a \code{\link{BioclimModel-class}} S4 object backed by a C++
#' \code{BioclimModel} instance. The four monthly climate arrays are validated
#' and passed to the C++ object; all bioclimatic variable computations delegate
#' to that object via Rcpp.
#'
#' @param tas    Numeric vector of length 12: monthly mean temperature.
#' @param tasmax Numeric vector of length 12: monthly maximum temperature.
#' @param tasmin Numeric vector of length 12: monthly minimum temperature.
#' @param pr     Numeric vector of length 12: monthly precipitation.
#'
#' @return A \code{\link{BioclimModel-class}} object.
#'
#' @export
#' @examples
#' tas    <- 1:12
#' tasmax <- 2:13
#' tasmin <- 0:11
#' pr     <- 1:12
#' m <- BioclimModel(tas, tasmax, tasmin, pr)
#' bio01(m)
#' bioclim(m)
BioclimModel <- function(tas, tasmax, tasmin, pr) {
  validate_monthly(tas,    "tas")
  validate_monthly(tasmax, "tasmax")
  validate_monthly(tasmin, "tasmin")
  validate_monthly(pr,     "pr")
  ptr <- bioclim_model_new(
    as.double(tas), as.double(tasmax),
    as.double(tasmin), as.double(pr)
  )
  obj <- methods::new("BioclimModel", pntr = ptr)
  methods::validObject(obj)
  obj
}

# ── show method ──────────────────────────────────────────────────────────────

setMethod("show", "BioclimModel", function(object) {
  cat("class       : BioclimModel\n")
  cat("pntr        : <C++ BioclimModel>\n")
  invisible(object)
})

# ── BioclimModel S4 methods ───────────────────────────────────────────────────
#
# The generics are defined in BioclimData.R (which loads first in the Collate
# order).  Only the BioclimModel-specific dispatch paths are registered here.
# All other dispatch (ANY, BioclimData) remains as set up in BioclimData.R.

#' S4 Methods for BioclimModel Objects
#'
#' S4 method implementations for all 19 bioclimatic variable functions and
#' \code{bioclim} that dispatch to the underlying C++ object via Rcpp when
#' the first argument is a \code{\link{BioclimModel-class}} instance.
#'
#' @name BioclimModel-methods
#' @aliases bio01,BioclimModel-method bio04,BioclimModel-method bio05,BioclimModel-method
#' @aliases bio06,BioclimModel-method bio10,BioclimModel-method bio11,BioclimModel-method
#' @aliases bio12,BioclimModel-method bio13,BioclimModel-method bio14,BioclimModel-method
#' @aliases bio15,BioclimModel-method bio16,BioclimModel-method bio17,BioclimModel-method
#' @aliases bio02,BioclimModel,missing-method bio03,BioclimModel,missing-method
#' @aliases bio07,BioclimModel,missing-method bio08,BioclimModel,missing-method
#' @aliases bio09,BioclimModel,missing-method bio18,BioclimModel,missing-method
#' @aliases bio19,BioclimModel,missing-method
#' @aliases bioclim,BioclimModel,ANY,ANY,ANY-method
NULL

# ── Single-argument BioclimModel methods ──────────────────────────────────────

setMethod("bio01", "BioclimModel", function(tas, ...)    bioclim_model_bio01(tas@pntr))
setMethod("bio04", "BioclimModel", function(tas, ...)    bioclim_model_bio04(tas@pntr))
setMethod("bio05", "BioclimModel", function(tasmax, ...) bioclim_model_bio05(tasmax@pntr))
setMethod("bio06", "BioclimModel", function(tasmin, ...) bioclim_model_bio06(tasmin@pntr))
setMethod("bio10", "BioclimModel", function(tas, ...)    bioclim_model_bio10(tas@pntr))
setMethod("bio11", "BioclimModel", function(tas, ...)    bioclim_model_bio11(tas@pntr))
setMethod("bio12", "BioclimModel", function(pr, ...)     bioclim_model_bio12(pr@pntr))
setMethod("bio13", "BioclimModel", function(pr, ...)     bioclim_model_bio13(pr@pntr))
setMethod("bio14", "BioclimModel", function(pr, ...)     bioclim_model_bio14(pr@pntr))
setMethod("bio15", "BioclimModel", function(pr, ...)     bioclim_model_bio15(pr@pntr))
setMethod("bio16", "BioclimModel", function(pr, ...)     bioclim_model_bio16(pr@pntr))
setMethod("bio17", "BioclimModel", function(pr, ...)     bioclim_model_bio17(pr@pntr))

# ── Two-argument BioclimModel methods: (tasmax, tasmin) ───────────────────────

setMethod("bio02", signature("BioclimModel", "missing"),
          function(tasmax, tasmin, ...) bioclim_model_bio02(tasmax@pntr))

setMethod("bio03", signature("BioclimModel", "missing"),
          function(tasmax, tasmin, ...) bioclim_model_bio03(tasmax@pntr))

setMethod("bio07", signature("BioclimModel", "missing"),
          function(tasmax, tasmin, ...) bioclim_model_bio07(tasmax@pntr))

# ── Two-argument BioclimModel methods: (tas, pr) ──────────────────────────────

setMethod("bio08", signature("BioclimModel", "missing"),
          function(tas, pr, ...) bioclim_model_bio08(tas@pntr))

setMethod("bio09", signature("BioclimModel", "missing"),
          function(tas, pr, ...) bioclim_model_bio09(tas@pntr))

setMethod("bio18", signature("BioclimModel", "missing"),
          function(tas, pr, ...) bioclim_model_bio18(tas@pntr))

setMethod("bio19", signature("BioclimModel", "missing"),
          function(tas, pr, ...) bioclim_model_bio19(tas@pntr))

# ── bioclim() BioclimModel method ─────────────────────────────────────────────
#
# The BioclimModel method ignores the extra arguments (tasmax, tasmin, pr)
# because all data are already stored in the C++ object.

setMethod("bioclim", "BioclimModel",
          function(tas, tasmax, tasmin, pr, ...)
            bioclim_model_compute(tas@pntr))
