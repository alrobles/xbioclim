# FindGDAL.cmake — Thin wrapper around CMake's built-in FindGDAL.
#
# Re-export the search so that `find_package(GDAL ...)` works even when
# CMAKE_MODULE_PATH includes this directory first.  We delegate to the
# standard module shipped with CMake.
#
# After the built-in module runs, GDAL_FOUND, GDAL_INCLUDE_DIRS and
# GDAL_LIBRARIES will be set.

# Avoid infinite recursion: temporarily remove our directory from the
# module path so the built-in FindGDAL is found instead.
list(REMOVE_ITEM CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")
find_package(GDAL ${GDAL_FIND_VERSION} ${GDAL_FIND_REQUIRED_ARG} QUIET)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")
