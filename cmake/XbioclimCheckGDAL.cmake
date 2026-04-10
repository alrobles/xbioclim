# XbioclimCheckGDAL.cmake — verify GDAL is available with a helpful message
#
# Include this module after find_package(GDAL ...) to get a clear error
# message pointing users at the right install command.

if(NOT GDAL_FOUND)
    message(FATAL_ERROR
        "GDAL (>= 3.0) not found.\n"
        "Install it via your package manager:\n"
        "  Ubuntu/Debian: sudo apt install libgdal-dev\n"
        "  macOS (Homebrew): brew install gdal\n"
        "Or set GDAL_DIR / CMAKE_PREFIX_PATH to a custom installation.")
endif()
