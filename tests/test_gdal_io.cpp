#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "xbioclim/gdal_io.hpp"

#include <filesystem>
#include <string>
#include <vector>

using Catch::Matchers::WithinAbs;
using namespace xbioclim;

// ---------------------------------------------------------------------------
// Helper: return the path to the tests/data directory relative to the build.
// CI generates test fixtures using tools/generate_test_data.py before running
// ctest, so the directory is expected to exist if the CI workflow is followed.
// ---------------------------------------------------------------------------
static std::filesystem::path data_dir() {
    // Allow override via environment variable for flexibility
    const char* env = std::getenv("XBIOCLIM_TEST_DATA");
    if (env) return std::filesystem::path(env);
    // Default: tests/data relative to the source root (set via CMake or cwd)
    return std::filesystem::path("tests/data");
}

static bool test_data_available() {
    auto p = data_dir() / "tas_01.tif";
    return std::filesystem::exists(p);
}

static std::vector<std::string> make_file_list(const std::string& var) {
    std::vector<std::string> files;
    auto base = data_dir();
    for (int m = 1; m <= 12; ++m) {
        char buf[8];
        std::snprintf(buf, sizeof(buf), "%02d", m);
        files.push_back((base / (var + "_" + buf + ".tif")).string());
    }
    return files;
}

// ---------------------------------------------------------------------------
// Tests that require the synthetic GeoTIFF fixtures to be present.
// They are tagged [gdal_io] and will be skipped gracefully if data is absent.
// ---------------------------------------------------------------------------

TEST_CASE("GdalReader opens valid datasets", "[gdal_io]") {
    if (!test_data_available()) {
        SKIP("Test data not found; run tools/generate_test_data.py first");
    }

    REQUIRE_NOTHROW(GdalReader(
        make_file_list("tas"),
        make_file_list("tasmax"),
        make_file_list("tasmin"),
        make_file_list("pr")));
}

TEST_CASE("GdalReader reports correct raster dimensions", "[gdal_io]") {
    if (!test_data_available()) {
        SKIP("Test data not found; run tools/generate_test_data.py first");
    }

    GdalReader reader(
        make_file_list("tas"),
        make_file_list("tasmax"),
        make_file_list("tasmin"),
        make_file_list("pr"));

    // The synthetic fixtures are 10×10
    CHECK(reader.raster_x_size() == 10);
    CHECK(reader.raster_y_size() == 10);
}

TEST_CASE("GdalReader read_block returns correct shape", "[gdal_io]") {
    if (!test_data_available()) {
        SKIP("Test data not found; run tools/generate_test_data.py first");
    }

    GdalReader reader(
        make_file_list("tas"),
        make_file_list("tasmax"),
        make_file_list("tasmin"),
        make_file_list("pr"));

    ClimateBlock block = reader.read_block(0, 0, 10, 10);

    CHECK(block.tas.shape(0) == 100u);    // 10*10 pixels
    CHECK(block.tas.shape(1) == 12u);     // 12 months
    CHECK(block.tasmax.shape(0) == 100u);
    CHECK(block.tasmin.shape(0) == 100u);
    CHECK(block.pr.shape(0) == 100u);
}

TEST_CASE("GdalReader read_block returns expected pixel values", "[gdal_io]") {
    if (!test_data_available()) {
        SKIP("Test data not found; run tools/generate_test_data.py first");
    }

    GdalReader reader(
        make_file_list("tas"),
        make_file_list("tasmax"),
        make_file_list("tasmin"),
        make_file_list("pr"));

    ClimateBlock block = reader.read_block(0, 0, 10, 10);

    // All pixels are uniform; month 0 (January) should be:
    //   tas[0]    = 1.0
    //   tasmax[0] = 2.0
    //   tasmin[0] = 0.0
    //   pr[0]     = 1.0
    CHECK_THAT(block.tas   (0, 0), WithinAbs(1.0f, 1e-3f));
    CHECK_THAT(block.tasmax(0, 0), WithinAbs(2.0f, 1e-3f));
    CHECK_THAT(block.tasmin(0, 0), WithinAbs(0.0f, 1e-3f));
    CHECK_THAT(block.pr    (0, 0), WithinAbs(1.0f, 1e-3f));

    // Month 11 (December):
    //   tas[11] = 12.0
    CHECK_THAT(block.tas(0, 11), WithinAbs(12.0f, 1e-3f));
}
