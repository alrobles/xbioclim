# Developer convenience makefile.  Not used by R CMD build/INSTALL.
# This file is excluded from CRAN source tarballs via .Rbuildignore.

.PHONY: all clean

all:
	@:

clean:
	Rscript tools/clean-obj.R
