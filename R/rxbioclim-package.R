#' rxbioclim: Bioclimatic Variables from Monthly Climate Data
#'
#' Computes the 19 standard bioclimatic variables (BIO01-BIO19) from monthly
#' climate data following the WorldClim specification. This is an R
#' implementation of the xbioclim C++ library, with a compiled C++ back-end
#' exposed through Rcpp Modules.
#'
#' @docType package
#' @name rxbioclim-package
#' @useDynLib rxbioclim
#' @importFrom Rcpp evalCpp
#' @importFrom methods new
"_PACKAGE"
