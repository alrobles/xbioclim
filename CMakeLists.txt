cmake_minimum_required(VERSION 3.15)
project(xbioclim VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Custom CMake module path ---
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# --- Dependencies ---
find_package(GDAL 3.0 REQUIRED)
include(XbioclimCheckGDAL)
find_package(xtensor REQUIRED)
include(XbioclimCheckXtensor)

option(XBIOCLIM_USE_OPENMP "Enable OpenMP CPU parallelism" ON)
if(XBIOCLIM_USE_OPENMP)
    find_package(OpenMP)
endif()

option(XBIOCLIM_USE_OPENMP_OFFLOAD "Enable OpenMP 4.5 GPU offload" OFF)
if(XBIOCLIM_USE_OPENMP_OFFLOAD)
    if(NOT OpenMP_CXX_FOUND)
        find_package(OpenMP REQUIRED)
    endif()
endif()

option(XBIOCLIM_USE_CUDA "Enable CUDA GPU kernels" OFF)
if(XBIOCLIM_USE_CUDA)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

# --- Library target ---
add_library(xbioclim
    src/primitives.cpp
    src/bioclim.cpp
    src/gdal_io.cpp
)

if(XBIOCLIM_USE_CUDA)
    target_sources(xbioclim PRIVATE src/bioclim_cuda.cu)
    set_target_properties(xbioclim PROPERTIES CUDA_ARCHITECTURES "80;86;90")
    target_compile_definitions(xbioclim PUBLIC XBIOCLIM_USE_CUDA)
endif()

target_include_directories(xbioclim PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(xbioclim PUBLIC
    GDAL::GDAL
    xtensor
    $<$<AND:$<BOOL:${XBIOCLIM_USE_OPENMP}>,$<TARGET_EXISTS:OpenMP::OpenMP_CXX>>:OpenMP::OpenMP_CXX>
    $<$<AND:$<BOOL:${XBIOCLIM_USE_OPENMP_OFFLOAD}>,$<NOT:$<BOOL:${XBIOCLIM_USE_OPENMP}>>,$<TARGET_EXISTS:OpenMP::OpenMP_CXX>>:OpenMP::OpenMP_CXX>
)

if(XBIOCLIM_USE_OPENMP AND OpenMP_CXX_FOUND)
    target_compile_definitions(xbioclim PUBLIC XBIOCLIM_USE_OPENMP)
endif()

if(XBIOCLIM_USE_OPENMP_OFFLOAD AND OpenMP_CXX_FOUND)
    target_compile_definitions(xbioclim PUBLIC XBIOCLIM_USE_OPENMP_OFFLOAD)
    # Allow users to specify GPU offload targets via -DXBIOCLIM_OPENMP_OFFLOAD_FLAGS=...
    # Examples:
    #   GCC  (NVIDIA): -foffload=nvptx-none -foffload-options=nvptx-none=-misa=sm_80
    #   Clang (NVIDIA): -fopenmp-targets=nvptx64-nvidia-cuda
    #   Clang (AMD):    -fopenmp-targets=amdgcn-amd-amdhsa
    # When empty (default), OpenMP target regions fall back to host CPU execution.
    if(XBIOCLIM_OPENMP_OFFLOAD_FLAGS)
        target_compile_options(xbioclim PRIVATE ${XBIOCLIM_OPENMP_OFFLOAD_FLAGS})
        target_link_options(xbioclim PRIVATE ${XBIOCLIM_OPENMP_OFFLOAD_FLAGS})
    endif()
endif()

# --- Export target ---
add_library(xbioclim::xbioclim ALIAS xbioclim)

# --- CLI executable ---
add_executable(xbioclim_cli src/main.cpp)
target_link_libraries(xbioclim_cli PRIVATE xbioclim)

# --- Tests ---
option(XBIOCLIM_BUILD_TESTS "Build the test suite" ON)
if(XBIOCLIM_BUILD_TESTS)
    # --- Test data generator (test-only) ---
    add_executable(xbioclim_generate_test_data tools/generate_mock_tiff.cpp)
    target_link_libraries(xbioclim_generate_test_data PRIVATE GDAL::GDAL)

    enable_testing()
    add_subdirectory(tests)
endif()

# --- Install targets ---
include(GNUInstallDirs)

install(TARGETS xbioclim
    EXPORT xbioclimTargets
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT xbioclimTargets
    FILE xbioclimTargets.cmake
    NAMESPACE xbioclim::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/xbioclim
)
