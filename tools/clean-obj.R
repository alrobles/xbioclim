#!/usr/bin/env Rscript
# tools/clean-obj.R
#
# Remove build artifacts that roxygen2::roxygenise() may leave behind
# when it compiles the package with debug flags (-g -O0 -UNDEBUG).
# Run this from the package root before a release R CMD INSTALL.

root <- getwd()

stopifnot(
  "Run this script from the package root (it needs DESCRIPTION and src/)" =
    file.exists(file.path(root, "DESCRIPTION")) &&
    file.exists(file.path(root, "src"))
)

cat("Cleaning debug build artifacts from", root, "\n")

rm_pattern <- function(pattern, recursive = FALSE) {
  paths <- Sys.glob(file.path(root, pattern))
  if (length(paths) == 0L) {
    return(invisible(NULL))
  }
  for (p in paths) {
    unlink(p, recursive = recursive, force = TRUE)
  }
  invisible(NULL)
}

rm_pattern("src/*.o")
rm_pattern("src/*.so")
rm_pattern("src/*.dll")
rm_pattern("src/symbols.rds")
rm_pattern("src/Makevars")
rm_pattern("config.log")
rm_pattern("config.status")
rm_pattern("src/config.log")
rm_pattern("src/config.status")
rm_pattern("autom4te.cache", recursive = TRUE)
rm_pattern("tools/autom4te.cache", recursive = TRUE)
rm_pattern("src/autom4te.cache", recursive = TRUE)

invisible(NULL)
