/**
 * generate_mock_tiff.cpp — Create synthetic 10×10 GeoTIFF test fixtures.
 *
 * Produces 48 compressed Int16 GeoTIFFs (4 variables × 12 months) that
 * match the data layout described in docs/ROADMAP.md §1.1:
 *
 *   tas[m]    = m          (m = 1..12)
 *   tasmax[m] = m + 1
 *   tasmin[m] = m − 1
 *   pr[m]     = m
 *
 * Output filename convention: <var>_<MM>.tif  (e.g. tas_01.tif, pr_12.tif)
 *
 * Usage:
 *   generate_mock_tiff [outdir]   (default: tests/data)
 */

#include <gdal_priv.h>
#include <cpl_conv.h>
#include <ogr_spatialref.h>

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

static constexpr int XSIZE = 10;
static constexpr int YSIZE = 10;
static constexpr int16_t NODATA = -9999;

// Geotransform for a 10-pixel × 10-pixel global grid (36° × 18° cells).
static const double GT[6] = {-180.0, 36.0, 0.0, 90.0, 0.0, -18.0};

static void write_tif(const std::string& path, int16_t value) {
    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("GTiff");
    if (!drv) {
        std::cerr << "Error: GTiff driver not available\n";
        std::exit(1);
    }

    char opt_compress[] = "COMPRESS=LZW";
    char* opts[] = {opt_compress, nullptr};
    GDALDataset* ds = drv->Create(path.c_str(), XSIZE, YSIZE, 1,
                                  GDT_Int16, opts);
    if (!ds) {
        std::cerr << "Error: cannot create " << path << "\n";
        std::exit(1);
    }

    OGRSpatialReference srs;
    srs.importFromEPSG(4326);
    srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    char* wkt = nullptr;
    srs.exportToWkt(&wkt);
    ds->SetProjection(wkt);
    CPLFree(wkt);

    double gt[6];
    for (int i = 0; i < 6; ++i) gt[i] = GT[i];
    ds->SetGeoTransform(gt);

    std::vector<int16_t> buf(static_cast<std::size_t>(XSIZE) * YSIZE, value);
    GDALRasterBand* band = ds->GetRasterBand(1);
    band->SetNoDataValue(NODATA);
    CPLErr err = band->RasterIO(GF_Write, 0, 0, XSIZE, YSIZE,
                                buf.data(), XSIZE, YSIZE, GDT_Int16, 0, 0);
    if (err != CE_None) {
        std::cerr << "Error: RasterIO write failed for " << path << "\n";
        GDALClose(ds);
        std::exit(1);
    }

    GDALClose(ds);
}

int main(int argc, char* argv[]) {
    std::string outdir = "tests/data";
    if (argc > 1) {
        outdir = argv[1];
    }

    std::filesystem::create_directories(outdir);
    GDALAllRegister();

    struct VarDef {
        const char* name;
        int (*value_fn)(int);
    };

    VarDef vars[] = {
        {"tas",    [](int m) -> int { return m; }},
        {"tasmax", [](int m) -> int { return m + 1; }},
        {"tasmin", [](int m) -> int { return m - 1; }},
        {"pr",     [](int m) -> int { return m; }},
    };

    int count = 0;
    for (const auto& var : vars) {
        for (int m = 1; m <= 12; ++m) {
            char fname[64];
            std::snprintf(fname, sizeof(fname), "%s_%02d.tif", var.name, m);
            std::string path = outdir + "/" + fname;
            write_tif(path, static_cast<int16_t>(var.value_fn(m)));
            ++count;
        }
    }

    std::cout << "Generated " << count << " GeoTIFFs in " << outdir << "/\n";
    return 0;
}
