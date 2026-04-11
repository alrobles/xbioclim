# FindGDAL.cmake — Thin wrapper around CMake's built-in FindGDAL.
#
# Re-export the search so that `find_package(GDAL ...)` works even when
# CMAKE_MODULE_PATH includes this directory first.  We delegate to the
# standard module shipped with CMake.

# Avoid infinite recursion: temporarily remove our directory from the
# module path so the built-in FindGDAL is found instead.
set(_xbioclim_gdal_args "")
if(GDAL_FIND_REQUIRED)
    list(APPEND _xbioclim_gdal_args REQUIRED)
endif()

list(REMOVE_ITEM CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")
find_package(GDAL ${GDAL_FIND_VERSION} ${_xbioclim_gdal_args} QUIET)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

unset(_xbioclim_gdal_args)
