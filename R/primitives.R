#' Primitive Helper Functions for Bioclimatic Variable Computation
#'
#' Internal helper functions used by the bioclimatic variable functions.
#' These mirror the primitives in the xbioclim C++ library.
#'
#' @name primitives
#' @keywords internal
NULL

#' Population standard deviation
#'
#' Computes the population standard deviation (denominator N, not N-1)
#' matching the xbioclim convention.
#'
#' @param x A numeric vector.
#' @param na.rm Logical. If `TRUE`, missing values are skipped.
#' @return A single numeric value.
#' @keywords internal
sd_pop <- function(x, na.rm = FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  n <- length(x)
  if (n == 0L) return(NaN)
  sqrt(sum((x - mean(x))^2) / n)
}

#' Compute rolling quarter sums with circular wrapping
#'
#' For each starting month (1-12), computes the sum of 3 consecutive months
#' with circular wrapping (month 13 = month 1, month 14 = month 2).
#'
#' @param x A numeric vector of length 12 (monthly values).
#' @param na.rm Logical. If `TRUE`, missing values are skipped.
#' @return A numeric vector of length 12 with rolling quarter sums.
#' @keywords internal
rolling_quarter_sum <- function(x, na.rm = FALSE) {
  x_ext <- c(x, x[1:2])
  vapply(seq_len(12), function(i) sum(x_ext[i:(i + 2L)], na.rm = na.rm), numeric(1))
}

#' Compute rolling quarter means with circular wrapping
#'
#' For each starting month (1-12), computes the mean of 3 consecutive months
#' with circular wrapping (month 13 = month 1, month 14 = month 2).
#'
#' @param x A numeric vector of length 12 (monthly values).
#' @param na.rm Logical. If `TRUE`, missing values are skipped.
#' @return A numeric vector of length 12 with rolling quarter means.
#' @keywords internal
rolling_quarter_mean <- function(x, na.rm = FALSE) {
  x_ext <- c(x, x[1:2])
  vapply(seq_len(12), function(i) mean(x_ext[i:(i + 2L)], na.rm = na.rm), numeric(1))
}

#' Find the starting month of the quarter with the maximum sum
#'
#' @param x A numeric vector of length 12 (monthly values).
#' @param na.rm Logical. If `TRUE`, missing values are skipped; a quarter needs
#'   at least one non-missing value to be considered.
#' @return An integer (1-12) indicating the starting month.
#' @keywords internal
quarter_argmax <- function(x, na.rm = FALSE) {
  qsums <- rolling_quarter_sum(x, na.rm = na.rm)
  if (na.rm) {
    x_ext <- c(x, x[1:2])
    n_valid <- vapply(seq_len(12), function(i) sum(!is.na(x_ext[i:(i + 2L)])), integer(1))
    qsums[n_valid == 0L] <- -Inf
  }
  which.max(qsums)
}

#' Find the starting month of the quarter with the minimum sum
#'
#' @param x A numeric vector of length 12 (monthly values).
#' @param na.rm Logical. If `TRUE`, missing values are skipped; a quarter needs
#'   at least one non-missing value to be considered.
#' @return An integer (1-12) indicating the starting month.
#' @keywords internal
quarter_argmin <- function(x, na.rm = FALSE) {
  qsums <- rolling_quarter_sum(x, na.rm = na.rm)
  if (na.rm) {
    x_ext <- c(x, x[1:2])
    n_valid <- vapply(seq_len(12), function(i) sum(!is.na(x_ext[i:(i + 2L)])), integer(1))
    qsums[n_valid == 0L] <- Inf
  }
  which.min(qsums)
}

#' Get the 3-month values for a quarter starting at a given month
#'
#' @param x A numeric vector of length 12 (monthly values).
#' @param start The starting month (1-12).
#' @return A numeric vector of length 3.
#' @keywords internal
quarter_values <- function(x, start) {
  idx <- ((start - 1L):(start + 1L)) %% 12L + 1L
  x[idx]
}

#' Validate monthly climate input
#'
#' Checks that input is a numeric vector of length 12.
#'
#' @param x The input to validate.
#' @param name Name of the variable for error messages.
#' @return Invisible NULL. Throws an error if validation fails.
#' @keywords internal
validate_monthly <- function(x, name = "input") {
  if (!is.numeric(x)) {
    stop(sprintf("'%s' must be numeric", name), call. = FALSE)
  }
  if (length(x) != 12L) {
    stop(sprintf("'%s' must have length 12 (one value per month), got %d",
                 name, length(x)), call. = FALSE)
  }
  invisible(NULL)
}
