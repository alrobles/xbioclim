#include "xbioclim/gdal_io.hpp"

#include <gdal_priv.h>
#include <cpl_conv.h>

#include <array>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace xbioclim {

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

static constexpr int   NUM_VARS   = 4;   // tas, tasmax, tasmin, pr
static constexpr int   NUM_MONTHS = 12;
static constexpr float NODATA_OUT = std::numeric_limits<float>::quiet_NaN();

// ---------------------------------------------------------------------------
// GdalReader::Impl
// ---------------------------------------------------------------------------

struct GdalReader::Impl {
    // datasets[var][month], var: 0=tas 1=tasmax 2=tasmin 3=pr
    GDALDataset* ds[NUM_VARS][NUM_MONTHS] = {};

    int x_size = 0;
    int y_size = 0;
    std::string proj_wkt;
    std::array<double, 6> gt{};

    void open(int var_idx, const std::vector<std::string>& files) {
        if (static_cast<int>(files.size()) != NUM_MONTHS) {
            throw std::runtime_error(
                "Expected exactly 12 files per variable, got " +
                std::to_string(files.size()));
        }
        for (int m = 0; m < NUM_MONTHS; ++m) {
            ds[var_idx][m] = static_cast<GDALDataset*>(
                GDALOpen(files[m].c_str(), GA_ReadOnly));
            if (!ds[var_idx][m]) {
                throw std::runtime_error("Cannot open: " + files[m]);
            }
        }
    }

    void validate_and_cache_metadata() {
        // Use the first dataset as reference
        GDALDataset* ref = ds[0][0];
        x_size = ref->GetRasterXSize();
        y_size = ref->GetRasterYSize();
        proj_wkt = ref->GetProjectionRef();
        ref->GetGeoTransform(gt.data());

        for (int v = 0; v < NUM_VARS; ++v) {
            for (int m = 0; m < NUM_MONTHS; ++m) {
                GDALDataset* d = ds[v][m];
                if (d->GetRasterXSize() != x_size ||
                    d->GetRasterYSize() != y_size) {
                    throw std::runtime_error(
                        "Raster dimension mismatch across input files");
                }
            }
        }
    }
};

// ---------------------------------------------------------------------------
// GdalReader
// ---------------------------------------------------------------------------

GdalReader::GdalReader(const std::vector<std::string>& tas_files,
                       const std::vector<std::string>& tasmax_files,
                       const std::vector<std::string>& tasmin_files,
                       const std::vector<std::string>& pr_files)
    : impl_(new Impl())
{
    GDALAllRegister();
    impl_->open(0, tas_files);
    impl_->open(1, tasmax_files);
    impl_->open(2, tasmin_files);
    impl_->open(3, pr_files);
    impl_->validate_and_cache_metadata();
}

GdalReader::~GdalReader() {
    if (impl_) {
        for (int v = 0; v < NUM_VARS; ++v)
            for (int m = 0; m < NUM_MONTHS; ++m)
                if (impl_->ds[v][m])
                    GDALClose(impl_->ds[v][m]);
        delete impl_;
    }
}

ClimateBlock GdalReader::read_block(int x_off, int y_off,
                                    int block_w, int block_h) const {
    const std::size_t n_pixels = static_cast<std::size_t>(block_w) * block_h;

    ClimateBlock block;
    block.tas    = Array2D::from_shape({n_pixels, 12});
    block.tasmax = Array2D::from_shape({n_pixels, 12});
    block.tasmin = Array2D::from_shape({n_pixels, 12});
    block.pr     = Array2D::from_shape({n_pixels, 12});

    Array2D* arrays[NUM_VARS] = {
        &block.tas, &block.tasmax, &block.tasmin, &block.pr
    };

    std::vector<float> buf(n_pixels);

    for (int v = 0; v < NUM_VARS; ++v) {
        for (int m = 0; m < NUM_MONTHS; ++m) {
            GDALRasterBand* band = impl_->ds[v][m]->GetRasterBand(1);
            CPLErr err = band->RasterIO(
                GF_Read,
                x_off, y_off,
                block_w, block_h,
                buf.data(),
                block_w, block_h,
                GDT_Float32,
                0, 0);
            if (err != CE_None) {
                throw std::runtime_error("RasterIO read failed");
            }
            // Copy buffer column (month m) into the array
            for (std::size_t p = 0; p < n_pixels; ++p) {
                (*arrays[v])(p, static_cast<std::size_t>(m)) = buf[p];
            }
        }
    }
    return block;
}

int GdalReader::raster_x_size() const { return impl_->x_size; }
int GdalReader::raster_y_size() const { return impl_->y_size; }
std::string GdalReader::projection_wkt() const { return impl_->proj_wkt; }
std::array<double, 6> GdalReader::geotransform() const { return impl_->gt; }

// ---------------------------------------------------------------------------
// GdalWriter::Impl
// ---------------------------------------------------------------------------

static constexpr int NUM_BIO = 19;

struct GdalWriter::Impl {
    GDALDataset* ds[NUM_BIO] = {};
};

// ---------------------------------------------------------------------------
// GdalWriter
// ---------------------------------------------------------------------------

GdalWriter::GdalWriter(const std::string& output_dir,
                       const std::string& prefix,
                       int x_size, int y_size,
                       const std::string& projection_wkt,
                       const std::array<double, 6>& gt)
    : impl_(new Impl())
{
    GDALAllRegister();
    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("GTiff");
    if (!drv) {
        throw std::runtime_error("GTiff driver not available");
    }

    char** opts = nullptr;
    opts = CSLSetNameValue(opts, "COMPRESS", "LZW");
    opts = CSLSetNameValue(opts, "TILED", "YES");
    opts = CSLSetNameValue(opts, "BLOCKXSIZE", "512");
    opts = CSLSetNameValue(opts, "BLOCKYSIZE", "512");
    opts = CSLSetNameValue(opts, "BIGTIFF", "IF_SAFER");

    for (int b = 0; b < NUM_BIO; ++b) {
        std::string path = output_dir + "/" + prefix +
                           std::to_string(b + 1) + ".tif";
        GDALDataset* out_ds = drv->Create(
            path.c_str(), x_size, y_size, 1, GDT_Float32, opts);
        if (!out_ds) {
            CSLDestroy(opts);
            throw std::runtime_error("Cannot create output: " + path);
        }
        out_ds->SetProjection(projection_wkt.c_str());
        double gt_arr[6];
        for (int i = 0; i < 6; ++i) gt_arr[i] = gt[i];
        out_ds->SetGeoTransform(gt_arr);
        out_ds->GetRasterBand(1)->SetNoDataValue(
            static_cast<double>(NODATA_OUT));
        impl_->ds[b] = out_ds;
    }

    CSLDestroy(opts);
}

GdalWriter::~GdalWriter() {
    if (impl_) {
        finalize();
        delete impl_;
    }
}

void GdalWriter::write_block(int x_off, int y_off, const BioBlock& bio) {
    const Array1D* arrays[NUM_BIO] = {
        &bio.bio01, &bio.bio02, &bio.bio03, &bio.bio04, &bio.bio05,
        &bio.bio06, &bio.bio07, &bio.bio08, &bio.bio09, &bio.bio10,
        &bio.bio11, &bio.bio12, &bio.bio13, &bio.bio14, &bio.bio15,
        &bio.bio16, &bio.bio17, &bio.bio18, &bio.bio19
    };

    // Determine block dimensions from the first BIO array size
    const std::size_t n_pixels = arrays[0]->size();
    // Infer block_w × block_h: we assume the caller knows the geometry, but
    // we need block_w and block_h. We derive them from the dataset dimensions
    // and the offsets (callers must supply consistent geometry).
    const int ds_x = impl_->ds[0]->GetRasterXSize();
    const int block_w = std::min(ds_x - x_off,
                                 static_cast<int>(n_pixels));
    const int block_h = static_cast<int>(n_pixels) / block_w;

    for (int b = 0; b < NUM_BIO; ++b) {
        const Array1D& arr = *arrays[b];
        const float* data_ptr = arr.data();
        GDALRasterBand* band = impl_->ds[b]->GetRasterBand(1);
        CPLErr err = band->RasterIO(
            GF_Write,
            x_off, y_off,
            block_w, block_h,
            const_cast<float*>(data_ptr),
            block_w, block_h,
            GDT_Float32,
            0, 0);
        if (err != CE_None) {
            throw std::runtime_error("RasterIO write failed for BIO" +
                                     std::to_string(b + 1));
        }
    }
}

void GdalWriter::finalize() {
    for (int b = 0; b < NUM_BIO; ++b) {
        if (impl_->ds[b]) {
            GDALClose(impl_->ds[b]);
            impl_->ds[b] = nullptr;
        }
    }
}

} // namespace xbioclim
