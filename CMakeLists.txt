cmake_minimum_required(VERSION 3.15)
project(xbioclim VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Dependencies ---
find_package(GDAL 3.0 REQUIRED)
find_package(xtensor REQUIRED)

option(XBIOCLIM_USE_OPENMP "Enable OpenMP CPU parallelism" ON)
if(XBIOCLIM_USE_OPENMP)
    find_package(OpenMP)
endif()

option(XBIOCLIM_USE_CUDA "Enable CUDA GPU kernels" OFF)
if(XBIOCLIM_USE_CUDA)
    enable_language(CUDA)
endif()

# --- Library target ---
add_library(xbioclim
    src/primitives.cpp
    src/bioclim.cpp
    src/gdal_io.cpp
)

target_include_directories(xbioclim PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(xbioclim PUBLIC
    GDAL::GDAL
    xtensor
    $<$<AND:$<BOOL:${XBIOCLIM_USE_OPENMP}>,$<TARGET_EXISTS:OpenMP::OpenMP_CXX>>:OpenMP::OpenMP_CXX>
)

if(XBIOCLIM_USE_OPENMP AND OpenMP_CXX_FOUND)
    target_compile_definitions(xbioclim PUBLIC XBIOCLIM_USE_OPENMP)
endif()

# --- Export target ---
add_library(xbioclim::xbioclim ALIAS xbioclim)

# --- Tests ---
option(XBIOCLIM_BUILD_TESTS "Build the test suite" ON)
if(XBIOCLIM_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
