#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "xbioclim/gdal_io.hpp"

#include <gdal_priv.h>
#include <cpl_conv.h>

#include <filesystem>
#include <string>
#include <vector>

using Catch::Matchers::WithinAbs;
using namespace xbioclim;

// ---------------------------------------------------------------------------
// Helper: return the path to the tests/data directory relative to the build.
// CI generates test fixtures using xbioclim_generate_test_data before running
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
        SKIP("Test data not found; run xbioclim_generate_test_data --outdir tests/data first");
    }

    REQUIRE_NOTHROW(GdalReader(
        make_file_list("tas"),
        make_file_list("tasmax"),
        make_file_list("tasmin"),
        make_file_list("pr")));
}

TEST_CASE("GdalReader reports correct raster dimensions", "[gdal_io]") {
    if (!test_data_available()) {
        SKIP("Test data not found; run xbioclim_generate_test_data --outdir tests/data first");
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
        SKIP("Test data not found; run xbioclim_generate_test_data --outdir tests/data first");
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

TEST_CASE("GdalWriter roundtrip writes and reads back correct values",
          "[gdal_io]") {
    // Create a temporary directory for output
    auto tmpdir = std::filesystem::temp_directory_path() / "xbioclim_test_write";
    std::filesystem::create_directories(tmpdir);
    REQUIRE(std::filesystem::exists(tmpdir));

    // Build a known BioBlock with 4 pixels (2x2 tile)
    const std::size_t N = 4;
    BioBlock bio;
    auto init = [&](Array1D& arr, float val) {
        arr = Array1D::from_shape({N});
        arr.fill(val);
    };
    init(bio.bio01, 6.5f);
    init(bio.bio02, 2.0f);
    init(bio.bio03, 15.38f);
    init(bio.bio04, 345.0f);
    init(bio.bio05, 13.0f);
    init(bio.bio06, 0.0f);
    init(bio.bio07, 13.0f);
    init(bio.bio08, 11.0f);
    init(bio.bio09, 2.0f);
    init(bio.bio10, 11.0f);
    init(bio.bio11, 2.0f);
    init(bio.bio12, 78.0f);
    init(bio.bio13, 12.0f);
    init(bio.bio14, 1.0f);
    init(bio.bio15, 53.1f);
    init(bio.bio16, 11.0f);
    init(bio.bio17, 2.0f);
    init(bio.bio18, 11.0f);
    init(bio.bio19, 2.0f);

    // Write
    std::array<double, 6> gt = {0.0, 1.0, 0.0, 0.0, 0.0, -1.0};
    {
        GdalWriter writer(tmpdir.string(), "bio", 2, 2, "", gt);
        writer.write_block(0, 0, 2, 2, bio);
        writer.finalize();
    }

    // Read back BIO01 and BIO12 with GDAL directly
    GDALAllRegister();
    {
        auto path1 = (tmpdir / "bio01.tif").string();
        INFO("Checking BIO01 file: " << path1);
        REQUIRE(std::filesystem::exists(path1));
        GDALDataset* ds = static_cast<GDALDataset*>(
            GDALOpen(path1.c_str(), GA_ReadOnly));
        REQUIRE(ds != nullptr);
        CHECK(ds->GetRasterXSize() == 2);
        CHECK(ds->GetRasterYSize() == 2);
        float buf[4];
        CPLErr err1 = ds->GetRasterBand(1)->RasterIO(
            GF_Read, 0, 0, 2, 2, buf, 2, 2, GDT_Float32, 0, 0);
        REQUIRE(err1 == CE_None);
        for (int i = 0; i < 4; ++i)
            CHECK_THAT(buf[i], WithinAbs(6.5f, 1e-3f));
        GDALClose(ds);
    }
    {
        auto path12 = (tmpdir / "bio12.tif").string();
        INFO("Checking BIO12 file: " << path12);
        REQUIRE(std::filesystem::exists(path12));
        GDALDataset* ds = static_cast<GDALDataset*>(
            GDALOpen(path12.c_str(), GA_ReadOnly));
        REQUIRE(ds != nullptr);
        float buf[4];
        CPLErr err2 = ds->GetRasterBand(1)->RasterIO(
            GF_Read, 0, 0, 2, 2, buf, 2, 2, GDT_Float32, 0, 0);
        REQUIRE(err2 == CE_None);
        for (int i = 0; i < 4; ++i)
            CHECK_THAT(buf[i], WithinAbs(78.0f, 1e-3f));
        GDALClose(ds);
    }

    // Cleanup
    std::filesystem::remove_all(tmpdir);
}

TEST_CASE("GdalReader read_block returns expected pixel values", "[gdal_io]") {
    if (!test_data_available()) {
        SKIP("Test data not found; run xbioclim_generate_test_data --outdir tests/data first");
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
