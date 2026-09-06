// gdal_io.cpp — implementation of GdalReader and GdalWriter.
//
// All GDAL-dependent code is inside #ifdef HAVE_GDAL blocks.  When GDAL is
// absent the constructors throw std::runtime_error with a clear message.

#include "gdal_io.hpp"

#include <sstream>
#include <stdexcept>

#ifdef HAVE_GDAL
#include <gdal_priv.h>
#include <cpl_conv.h>
#include <ogr_spatialref.h>
#endif

namespace xbioclim {

// ============================================================================
// GdalReader implementation
// ============================================================================

GdalReader::GdalReader(const std::string& path)
    : path_(path)
{
#ifdef HAVE_GDAL
    GDALAllRegister();
    ds_ = static_cast<GDALDataset*>(GDALOpen(path.c_str(), GA_ReadOnly));
    if (ds_ == nullptr) {
        std::ostringstream oss;
        oss << "GdalReader: cannot open '" << path << "': "
            << CPLGetLastErrorMsg();
        throw std::runtime_error(oss.str());
    }
#else
    (void)path;
    throw std::runtime_error(
        "GdalReader: xbioclim was built without GDAL support. "
        "Install GDAL >= 2.0.1 and reinstall the package."
    );
#endif
}

GdalReader::~GdalReader() {
#ifdef HAVE_GDAL
    if (ds_ != nullptr) {
        GDALClose(ds_);
        ds_ = nullptr;
    }
#endif
}

int GdalReader::nrows() const {
#ifdef HAVE_GDAL
    return ds_->GetRasterYSize();
#else
    return 0;
#endif
}

int GdalReader::ncols() const {
#ifdef HAVE_GDAL
    return ds_->GetRasterXSize();
#else
    return 0;
#endif
}

int GdalReader::nbands() const {
#ifdef HAVE_GDAL
    return ds_->GetRasterCount();
#else
    return 0;
#endif
}

std::vector<double> GdalReader::geotransform() const {
#ifdef HAVE_GDAL
    std::vector<double> gt(6, 0.0);
    ds_->GetGeoTransform(gt.data());
    return gt;
#else
    return std::vector<double>(6, 0.0);
#endif
}

std::string GdalReader::crs() const {
#ifdef HAVE_GDAL
    const char* wkt = ds_->GetProjectionRef();
    if (wkt == nullptr) return "";
    return std::string(wkt);
#else
    return "";
#endif
}

double GdalReader::scale(int band) const {
#ifdef HAVE_GDAL
    GDALRasterBand* b = ds_->GetRasterBand(band);
    if (b == nullptr) return 1.0;
    int has_scale = 0;
    double s = b->GetScale(&has_scale);
    return has_scale ? s : 1.0;
#else
    (void)band;
    return 1.0;
#endif
}

double GdalReader::offset(int band) const {
#ifdef HAVE_GDAL
    GDALRasterBand* b = ds_->GetRasterBand(band);
    if (b == nullptr) return 0.0;
    int has_offset = 0;
    double o = b->GetOffset(&has_offset);
    return has_offset ? o : 0.0;
#else
    (void)band;
    return 0.0;
#endif
}

void GdalReader::read_window(int xoff, int yoff,
                             int xsize, int ysize,
                             int band,
                             std::vector<double>& buf) const {
#ifdef HAVE_GDAL
    GDALRasterBand* b = ds_->GetRasterBand(band);
    if (b == nullptr) {
        std::ostringstream oss;
        oss << "GdalReader::read_window: invalid band " << band;
        throw std::runtime_error(oss.str());
    }

    const int n = xsize * ysize;
    buf.resize(static_cast<std::size_t>(n));

    CPLErr err = b->RasterIO(GF_Read,
                             xoff, yoff, xsize, ysize,
                             buf.data(),
                             xsize, ysize,
                             GDT_Float64,
                             0, 0);
    if (err != CE_None) {
        std::ostringstream oss;
        oss << "GdalReader::read_window: RasterIO failed for band " << band
            << ": " << CPLGetLastErrorMsg();
        throw std::runtime_error(oss.str());
    }

    // Apply scale / offset if present (physical = raw * scale + offset).
    int has_scale  = 0;
    int has_offset = 0;
    double sc = b->GetScale(&has_scale);
    double of = b->GetOffset(&has_offset);

    if (has_scale || has_offset) {
        if (!has_scale)  sc = 1.0;
        if (!has_offset) of = 0.0;
        for (int i = 0; i < n; ++i) {
            buf[static_cast<std::size_t>(i)] =
                buf[static_cast<std::size_t>(i)] * sc + of;
        }
    }
#else
    (void)xoff; (void)yoff; (void)xsize; (void)ysize; (void)band;
    buf.clear();
    throw std::runtime_error(
        "GdalReader::read_window: xbioclim was built without GDAL support."
    );
#endif
}

// ============================================================================
// GdalWriter implementation
// ============================================================================

GdalWriter::GdalWriter(const std::string& path,
                       int nrows, int ncols, int nbands,
                       const std::vector<double>& geotransform,
                       const std::string& crs,
                       bool cog_compatible)
    : closed_(false)
{
#ifdef HAVE_GDAL
    GDALAllRegister();

    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GTiff");
    if (driver == nullptr) {
        throw std::runtime_error("GdalWriter: GTiff driver not available.");
    }

    // Build creation options.
    char** opts = nullptr;
    if (cog_compatible) {
        opts = CSLSetNameValue(opts, "TILED",    "YES");
        opts = CSLSetNameValue(opts, "BLOCKXSIZE", "256");
        opts = CSLSetNameValue(opts, "BLOCKYSIZE", "256");
        opts = CSLSetNameValue(opts, "COMPRESS", "LZW");
        opts = CSLSetNameValue(opts, "BIGTIFF",  "IF_SAFER");
    }

    ds_ = driver->Create(path.c_str(), ncols, nrows, nbands,
                         GDT_Float64, opts);
    CSLDestroy(opts);

    if (ds_ == nullptr) {
        std::ostringstream oss;
        oss << "GdalWriter: cannot create '" << path << "': "
            << CPLGetLastErrorMsg();
        throw std::runtime_error(oss.str());
    }

    // Set geotransform (6 elements required).
    if (geotransform.size() == 6) {
        std::vector<double> gt_copy(geotransform);
        ds_->SetGeoTransform(gt_copy.data());
    }

    // Set CRS.
    if (!crs.empty()) {
        ds_->SetProjection(crs.c_str());
    }
#else
    (void)path; (void)nrows; (void)ncols; (void)nbands;
    (void)geotransform; (void)crs; (void)cog_compatible;
    throw std::runtime_error(
        "GdalWriter: xbioclim was built without GDAL support. "
        "Install GDAL >= 2.0.1 and reinstall the package."
    );
#endif
}

GdalWriter::~GdalWriter() {
    try { close(); } catch (...) {}
}

void GdalWriter::write_window(int xoff, int yoff,
                              int xsize, int ysize,
                              int band,
                              const std::vector<double>& buf) {
#ifdef HAVE_GDAL
    if (closed_) {
        throw std::runtime_error(
            "GdalWriter::write_window: dataset is already closed.");
    }

    GDALRasterBand* b = ds_->GetRasterBand(band);
    if (b == nullptr) {
        std::ostringstream oss;
        oss << "GdalWriter::write_window: invalid band " << band;
        throw std::runtime_error(oss.str());
    }

    const int n = xsize * ysize;
    if (static_cast<int>(buf.size()) < n) {
        std::ostringstream oss;
        oss << "GdalWriter::write_window: buffer too small ("
            << buf.size() << " < " << n << ")";
        throw std::runtime_error(oss.str());
    }

    CPLErr err = b->RasterIO(GF_Write,
                             xoff, yoff, xsize, ysize,
                             const_cast<double*>(buf.data()),
                             xsize, ysize,
                             GDT_Float64,
                             0, 0);
    if (err != CE_None) {
        std::ostringstream oss;
        oss << "GdalWriter::write_window: RasterIO failed for band " << band
            << ": " << CPLGetLastErrorMsg();
        throw std::runtime_error(oss.str());
    }
#else
    (void)xoff; (void)yoff; (void)xsize; (void)ysize; (void)band; (void)buf;
    throw std::runtime_error(
        "GdalWriter::write_window: xbioclim was built without GDAL support."
    );
#endif
}

void GdalWriter::close() {
#ifdef HAVE_GDAL
    if (!closed_ && ds_ != nullptr) {
        GDALClose(ds_);
        ds_     = nullptr;
        closed_ = true;
    }
#else
    closed_ = true;
#endif
}

} // namespace xbioclim
