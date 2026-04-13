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

# ── Convert existing plain functions to S4 generics ──────────────────────────
#
# The original plain functions from bioclim.R are captured first, then each
# name is promoted to an S4 generic.  A catch-all "ANY" or signature-based
# default method preserves the existing behaviour (including input validation
# and error messages) for plain numeric vector inputs; a "BioclimModel" method
# delegates to the C++ object.
#
# Parameter names in each setGeneric() call intentionally match those in the
# original plain function so that existing roxygen2 documentation remains
# consistent with the generic's \usage entry.
#
# No @export or @describeIn tags are needed here: the existing @export tags on
# the plain functions in bioclim.R already arrange for NAMESPACE entries, and
# those entries remain valid once the names are promoted to S4 generics.

.bio01_fn   <- bio01
.bio02_fn   <- bio02
.bio03_fn   <- bio03
.bio04_fn   <- bio04
.bio05_fn   <- bio05
.bio06_fn   <- bio06
.bio07_fn   <- bio07
.bio08_fn   <- bio08
.bio09_fn   <- bio09
.bio10_fn   <- bio10
.bio11_fn   <- bio11
.bio12_fn   <- bio12
.bio13_fn   <- bio13
.bio14_fn   <- bio14
.bio15_fn   <- bio15
.bio16_fn   <- bio16
.bio17_fn   <- bio17
.bio18_fn   <- bio18
.bio19_fn   <- bio19
.bioclim_fn <- bioclim

# ── Method documentation ──────────────────────────────────────────────────────

#' S4 Methods for BioclimModel Objects and Plain Numeric Vectors
#'
#' S4 method implementations for all 19 bioclimatic variable functions and
#' \code{bioclim}. Two dispatch paths exist for each:
#' \itemize{
#'   \item A \code{\link{BioclimModel-class}} path that delegates the
#'     computation to the underlying C++ object via Rcpp.
#'   \item A default (\code{ANY}) path that calls the original plain-function
#'     implementation, preserving full backward compatibility.
#' }
#'
#' @aliases
#'   bio01,ANY-method              bio01,BioclimModel-method
#'   bio02,ANY,ANY-method          bio02,BioclimModel,missing-method
#'   bio03,ANY,ANY-method          bio03,BioclimModel,missing-method
#'   bio04,ANY-method              bio04,BioclimModel-method
#'   bio05,ANY-method              bio05,BioclimModel-method
#'   bio06,ANY-method              bio06,BioclimModel-method
#'   bio07,ANY,ANY-method          bio07,BioclimModel,missing-method
#'   bio08,ANY,ANY-method          bio08,BioclimModel,missing-method
#'   bio08,BioclimModel,NULL-method
#'   bio09,ANY,ANY-method          bio09,BioclimModel,missing-method
#'   bio09,BioclimModel,NULL-method
#'   bio10,ANY-method              bio10,BioclimModel-method
#'   bio11,ANY-method              bio11,BioclimModel-method
#'   bio12,ANY-method              bio12,BioclimModel-method
#'   bio13,ANY-method              bio13,BioclimModel-method
#'   bio14,ANY-method              bio14,BioclimModel-method
#'   bio15,ANY-method              bio15,BioclimModel-method
#'   bio16,ANY-method              bio16,BioclimModel-method
#'   bio17,ANY-method              bio17,BioclimModel-method
#'   bio18,ANY,ANY-method          bio18,BioclimModel,missing-method
#'   bio18,BioclimModel,NULL-method
#'   bio19,ANY,ANY-method          bio19,BioclimModel,missing-method
#'   bio19,BioclimModel,NULL-method
#'   bioclim,ANY-method            bioclim,BioclimModel-method
#' @name BioclimModel-methods
NULL

# ── Single-argument generics ─────────────────────────────────────────────────

setGeneric("bio01", function(tas)    standardGeneric("bio01"))
setMethod("bio01", "ANY",          function(tas)    .bio01_fn(tas))
setMethod("bio01", "BioclimModel", function(tas)    bioclim_model_bio01(tas@pntr))

setGeneric("bio04", function(tas)    standardGeneric("bio04"))
setMethod("bio04", "ANY",          function(tas)    .bio04_fn(tas))
setMethod("bio04", "BioclimModel", function(tas)    bioclim_model_bio04(tas@pntr))

setGeneric("bio05", function(tasmax) standardGeneric("bio05"))
setMethod("bio05", "ANY",          function(tasmax) .bio05_fn(tasmax))
setMethod("bio05", "BioclimModel", function(tasmax) bioclim_model_bio05(tasmax@pntr))

setGeneric("bio06", function(tasmin) standardGeneric("bio06"))
setMethod("bio06", "ANY",          function(tasmin) .bio06_fn(tasmin))
setMethod("bio06", "BioclimModel", function(tasmin) bioclim_model_bio06(tasmin@pntr))

setGeneric("bio10", function(tas)    standardGeneric("bio10"))
setMethod("bio10", "ANY",          function(tas)    .bio10_fn(tas))
setMethod("bio10", "BioclimModel", function(tas)    bioclim_model_bio10(tas@pntr))

setGeneric("bio11", function(tas)    standardGeneric("bio11"))
setMethod("bio11", "ANY",          function(tas)    .bio11_fn(tas))
setMethod("bio11", "BioclimModel", function(tas)    bioclim_model_bio11(tas@pntr))

setGeneric("bio12", function(pr)     standardGeneric("bio12"))
setMethod("bio12", "ANY",          function(pr)     .bio12_fn(pr))
setMethod("bio12", "BioclimModel", function(pr)     bioclim_model_bio12(pr@pntr))

setGeneric("bio13", function(pr)     standardGeneric("bio13"))
setMethod("bio13", "ANY",          function(pr)     .bio13_fn(pr))
setMethod("bio13", "BioclimModel", function(pr)     bioclim_model_bio13(pr@pntr))

setGeneric("bio14", function(pr)     standardGeneric("bio14"))
setMethod("bio14", "ANY",          function(pr)     .bio14_fn(pr))
setMethod("bio14", "BioclimModel", function(pr)     bioclim_model_bio14(pr@pntr))

setGeneric("bio15", function(pr)     standardGeneric("bio15"))
setMethod("bio15", "ANY",          function(pr)     .bio15_fn(pr))
setMethod("bio15", "BioclimModel", function(pr)     bioclim_model_bio15(pr@pntr))

setGeneric("bio16", function(pr)     standardGeneric("bio16"))
setMethod("bio16", "ANY",          function(pr)     .bio16_fn(pr))
setMethod("bio16", "BioclimModel", function(pr)     bioclim_model_bio16(pr@pntr))

setGeneric("bio17", function(pr)     standardGeneric("bio17"))
setMethod("bio17", "ANY",          function(pr)     .bio17_fn(pr))
setMethod("bio17", "BioclimModel", function(pr)     bioclim_model_bio17(pr@pntr))

# ── Two-argument generics: (tasmax, tasmin) ───────────────────────────────────

setGeneric("bio02", function(tasmax, tasmin) standardGeneric("bio02"))
setMethod("bio02", signature("ANY", "ANY"),
          function(tasmax, tasmin) .bio02_fn(tasmax, tasmin))
setMethod("bio02", signature("BioclimModel", "missing"),
          function(tasmax, tasmin) bioclim_model_bio02(tasmax@pntr))

setGeneric("bio03", function(tasmax, tasmin) standardGeneric("bio03"))
setMethod("bio03", signature("ANY", "ANY"),
          function(tasmax, tasmin) .bio03_fn(tasmax, tasmin))
setMethod("bio03", signature("BioclimModel", "missing"),
          function(tasmax, tasmin) bioclim_model_bio03(tasmax@pntr))

setGeneric("bio07", function(tasmax, tasmin) standardGeneric("bio07"))
setMethod("bio07", signature("ANY", "ANY"),
          function(tasmax, tasmin) .bio07_fn(tasmax, tasmin))
setMethod("bio07", signature("BioclimModel", "missing"),
          function(tasmax, tasmin) bioclim_model_bio07(tasmax@pntr))

# ── Two-argument generics: (tas, pr) ─────────────────────────────────────────

setGeneric("bio08", function(tas, pr = NULL) standardGeneric("bio08"))
setMethod("bio08", signature("ANY", "ANY"),
          function(tas, pr = NULL) .bio08_fn(tas, pr))
setMethod("bio08", signature("BioclimModel", "missing"),
          function(tas, pr = NULL) bioclim_model_bio08(tas@pntr))
setMethod("bio08", signature("BioclimModel", "NULL"),
          function(tas, pr = NULL) bioclim_model_bio08(tas@pntr))

setGeneric("bio09", function(tas, pr = NULL) standardGeneric("bio09"))
setMethod("bio09", signature("ANY", "ANY"),
          function(tas, pr = NULL) .bio09_fn(tas, pr))
setMethod("bio09", signature("BioclimModel", "missing"),
          function(tas, pr = NULL) bioclim_model_bio09(tas@pntr))
setMethod("bio09", signature("BioclimModel", "NULL"),
          function(tas, pr = NULL) bioclim_model_bio09(tas@pntr))

setGeneric("bio18", function(tas, pr = NULL) standardGeneric("bio18"))
setMethod("bio18", signature("ANY", "ANY"),
          function(tas, pr = NULL) .bio18_fn(tas, pr))
setMethod("bio18", signature("BioclimModel", "missing"),
          function(tas, pr = NULL) bioclim_model_bio18(tas@pntr))
setMethod("bio18", signature("BioclimModel", "NULL"),
          function(tas, pr = NULL) bioclim_model_bio18(tas@pntr))

setGeneric("bio19", function(tas, pr = NULL) standardGeneric("bio19"))
setMethod("bio19", signature("ANY", "ANY"),
          function(tas, pr = NULL) .bio19_fn(tas, pr))
setMethod("bio19", signature("BioclimModel", "missing"),
          function(tas, pr = NULL) bioclim_model_bio19(tas@pntr))
setMethod("bio19", signature("BioclimModel", "NULL"),
          function(tas, pr = NULL) bioclim_model_bio19(tas@pntr))

# ── bioclim() generic ─────────────────────────────────────────────────────────
#
# The BioclimModel method ignores the three extra arguments (tasmax, tasmin,
# pr) because all data are already stored in the C++ object.  Default NULL
# values allow `bioclim(model)` to be called without those arguments.

setGeneric("bioclim",
           function(tas, tasmax = NULL, tasmin = NULL, pr = NULL)
             standardGeneric("bioclim"))
setMethod("bioclim", "ANY",
          function(tas, tasmax = NULL, tasmin = NULL, pr = NULL)
            .bioclim_fn(tas, tasmax, tasmin, pr))
setMethod("bioclim", "BioclimModel",
          function(tas, tasmax = NULL, tasmin = NULL, pr = NULL)
            bioclim_model_compute(tas@pntr))

# ── Re-register BioclimData methods ──────────────────────────────────────────
#
# setGeneric() calls above with NULL defaults change the generic signatures
# so that arguments like `pr` receive NULL (not missing) when omitted.
# The signature("BioclimData", "missing") methods from BioclimData.R no longer
# fire in that case. Register single-dispatch BioclimData methods here so they
# take precedence over ANY for all pr=NULL call patterns.

setMethod("bio03", signature("BioclimData", "missing"),
          function(tasmax, tasmin) bio03_cpp(tasmax@tasmax, tasmax@tasmin))

setMethod("bio08", "BioclimData",
          function(tas, pr = NULL) bio08_cpp(tas@tas, tas@pr))

setMethod("bio09", "BioclimData",
          function(tas, pr = NULL) bio09_cpp(tas@tas, tas@pr))

setMethod("bio18", "BioclimData",
          function(tas, pr = NULL) bio18_cpp(tas@tas, tas@pr))

setMethod("bio19", "BioclimData",
          function(tas, pr = NULL) bio19_cpp(tas@tas, tas@pr))

setMethod("bioclim", "BioclimData",
          function(tas, tasmax = NULL, tasmin = NULL, pr = NULL)
            bioclim_cpp(tas@tas, tas@tasmax, tas@tasmin, tas@pr))
