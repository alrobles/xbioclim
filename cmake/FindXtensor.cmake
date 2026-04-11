# FindXtensor.cmake — Thin wrapper so `find_package(xtensor ...)` resolves
# through CMAKE_MODULE_PATH even when this directory is listed.
#
# Delegates to the xtensorConfig.cmake shipped by the xtensor installation.

# Avoid infinite recursion: temporarily remove our directory from the
# module path so the upstream config package is found instead.
set(_xbioclim_xt_args "")
if(xtensor_FIND_REQUIRED)
    list(APPEND _xbioclim_xt_args REQUIRED)
endif()

list(REMOVE_ITEM CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")
find_package(xtensor ${xtensor_FIND_VERSION} ${_xbioclim_xt_args} QUIET CONFIG)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

unset(_xbioclim_xt_args)
