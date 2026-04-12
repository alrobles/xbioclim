#' Error and Warning Message Store
#'
#' @description
#' rxbioclim mirrors the `SpatMessages` pattern from the terra package to
#' provide a clean mechanism for propagating error and warning messages across
#' the C++/R boundary.
#'
#' **How the cross-boundary workflow works:**
#'
#' 1. A C++ routine performs its computation and, instead of throwing directly,
#'    records any error or warning strings into a session-level message store.
#' 2. After every `.Call()` invocation, the R-side wrapper calls
#'    [check_messages()] to inspect the store and re-raise any messages as
#'    native R conditions (`stop()` for errors, `warning()` for warnings).
#' 3. Users can also inspect the store directly with [bioclim_errors()] and
#'    [bioclim_warnings()] before R conditions are raised, or clear it with
#'    [clear_messages()].
#'
#' **Functions available to users:**
#' * [bioclim_errors()] – retrieve stored error messages.
#' * [bioclim_warnings()] – retrieve stored warning messages.
#' * [has_error()] – `TRUE` if at least one error is stored.
#' * [has_warning()] – `TRUE` if at least one warning is stored.
#' * [clear_messages()] – discard all stored messages.
#'
#' **Functions used internally (and by future C++ glue code):**
#' * `push_error()` – store an error message.
#' * `push_warning()` – store a warning message.
#' * `check_messages()` – raise stored messages as R conditions and clear them.
#'
#' @name rxbioclim-messages
NULL

# ---------------------------------------------------------------------------
# Internal message store
# ---------------------------------------------------------------------------

# A dedicated environment so messages survive calls without polluting the
# global namespace.  Analogous to a C++ SpatMessages struct kept alive for
# the duration of an R session.
.rxbioclim_env <- new.env(parent = emptyenv())
.rxbioclim_env$errors   <- character(0)
.rxbioclim_env$warnings <- character(0)

# ---------------------------------------------------------------------------
# Internal push functions  (called by C++ glue code via .Call, or by R
# helpers that need to record a deferred message)
# ---------------------------------------------------------------------------

#' Store an error message in the rxbioclim message store
#'
#' This function is called by C++ routines (via `.Call`) or internal R helpers
#' to record an error without immediately throwing an R condition.  Call
#' [check_messages()] afterwards to convert the stored message into a proper
#' R error.
#'
#' @param msg A single character string describing the error.
#' @return Invisible `NULL`.
#' @keywords internal
push_error <- function(msg) {
  stopifnot(is.character(msg), length(msg) == 1L)
  .rxbioclim_env$errors <- c(.rxbioclim_env$errors, msg)
  invisible(NULL)
}

#' Store a warning message in the rxbioclim message store
#'
#' This function is called by C++ routines (via `.Call`) or internal R helpers
#' to record a warning without immediately issuing an R condition.  Call
#' [check_messages()] afterwards to convert the stored message into a proper
#' R warning.
#'
#' @param msg A single character string describing the warning.
#' @return Invisible `NULL`.
#' @keywords internal
push_warning <- function(msg) {
  stopifnot(is.character(msg), length(msg) == 1L)
  .rxbioclim_env$warnings <- c(.rxbioclim_env$warnings, msg)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Internal check function (called by every R wrapper after a C++ call)
# ---------------------------------------------------------------------------

#' Propagate stored messages as R conditions
#'
#' Inspects the internal message store, issues any recorded warnings via
#' [base::warning()], clears them, then raises any recorded error via
#' [base::stop()] and clears it.  This function is called automatically by
#' every R wrapper function immediately after invoking a C++ routine, matching
#' the pattern used by the terra package.
#'
#' Calling this function when the store is empty is a no-op.
#'
#' @return Invisible `NULL` (unless an error is stored, in which case it
#'   throws).
#' @keywords internal
check_messages <- function() {
  if (has_warning()) {
    msgs <- .rxbioclim_env$warnings
    .rxbioclim_env$warnings <- character(0)
    for (msg in msgs) warning(msg, call. = FALSE)
  }
  if (has_error()) {
    msg <- paste(.rxbioclim_env$errors, collapse = "\n")
    .rxbioclim_env$errors <- character(0)
    stop(msg, call. = FALSE)
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Exported user-facing functions
# ---------------------------------------------------------------------------

#' Retrieve stored error messages
#'
#' Returns the character vector of error messages currently held in the
#' rxbioclim message store.  Under normal usage the store is automatically
#' flushed by [check_messages()] after every C++ call, but you can inspect it
#' manually before that point if needed.
#'
#' @return A character vector (possibly empty).
#' @export
#' @seealso [bioclim_warnings()], [has_error()], [clear_messages()]
#' @examples
#' clear_messages()
#' bioclim_errors()   # character(0)
bioclim_errors <- function() {
  .rxbioclim_env$errors
}

#' Retrieve stored warning messages
#'
#' Returns the character vector of warning messages currently held in the
#' rxbioclim message store.  Under normal usage the store is automatically
#' flushed by [check_messages()] after every C++ call, but you can inspect it
#' manually before that point if needed.
#'
#' @return A character vector (possibly empty).
#' @export
#' @seealso [bioclim_errors()], [has_warning()], [clear_messages()]
#' @examples
#' clear_messages()
#' bioclim_warnings()   # character(0)
bioclim_warnings <- function() {
  .rxbioclim_env$warnings
}

#' Check whether any errors are stored
#'
#' @return `TRUE` if the message store contains at least one error, `FALSE`
#'   otherwise.
#' @export
#' @seealso [bioclim_errors()], [has_warning()]
#' @examples
#' clear_messages()
#' has_error()   # FALSE
has_error <- function() {
  length(.rxbioclim_env$errors) > 0L
}

#' Check whether any warnings are stored
#'
#' @return `TRUE` if the message store contains at least one warning, `FALSE`
#'   otherwise.
#' @export
#' @seealso [bioclim_warnings()], [has_error()]
#' @examples
#' clear_messages()
#' has_warning()   # FALSE
has_warning <- function() {
  length(.rxbioclim_env$warnings) > 0L
}

#' Clear all stored messages
#'
#' Discards all error and warning messages currently held in the rxbioclim
#' message store.  This is called automatically by [check_messages()] after
#' propagating messages to R conditions, but you can call it manually to reset
#' state between operations.
#'
#' @return Invisible `NULL`.
#' @export
#' @seealso [bioclim_errors()], [bioclim_warnings()]
#' @examples
#' clear_messages()
#' has_error()    # FALSE
#' has_warning()  # FALSE
clear_messages <- function() {
  .rxbioclim_env$errors   <- character(0)
  .rxbioclim_env$warnings <- character(0)
  invisible(NULL)
}
