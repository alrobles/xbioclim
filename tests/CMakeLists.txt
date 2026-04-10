include(FetchContent)

find_package(Catch2 3 QUIET COMPONENTS Catch2 Catch2WithMain)
if(NOT Catch2_FOUND)
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.5.2
        GIT_SHALLOW    ON
    )
    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/extras")
endif()

add_executable(test_xbioclim
    test_primitives.cpp
    test_bioclim.cpp
    test_gdal_io.cpp
)

target_link_libraries(test_xbioclim PRIVATE
    xbioclim
    Catch2::Catch2WithMain
)

include(Catch)
# Catch2 v3 uses exit code 4 to signal skipped tests; tell CTest to treat it
# as a pass so that tests using SKIP() don't show as failures.
catch_discover_tests(test_xbioclim
    PROPERTIES SKIP_RETURN_CODE 4
)
