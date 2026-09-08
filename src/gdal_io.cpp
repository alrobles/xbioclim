// gdal_io.cpp — implementation of GdalReader and GdalWriter.
//
// All GDAL-dependent code is inside #ifdef HAVE_GDAL blocks.  When GDAL is
// absent the constructors throw std::runtime_error with a clear message.

#include "gdal_io.hpp"

#include <algorithm>
#include <sstream>
#include <stdexcept>

#ifdef HAVE_GDAL
#include <gdal_priv.h>
#include <cpl_conv.h>
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
    const std::size_t n = static_cast<std::size_t>(xsize) * ysize;
    buf.resize(n);
    read_window(xoff, yoff, xsize, ysize, band, buf.data());
}

void GdalReader::read_window(int xoff, int yoff,
                             int xsize, int ysize,
                             int band,
                             double* out) const {
#ifdef HAVE_GDAL
    GDALRasterBand* b = ds_->GetRasterBand(band);
    if (b == nullptr) {
        std::ostringstream oss;
        oss << "GdalReader::read_window: invalid band " << band;
        throw std::runtime_error(oss.str());
    }

    const int n = xsize * ysize;

    CPLErr err = b->RasterIO(GF_Read,
                             xoff, yoff, xsize, ysize,
                             out,
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
            out[i] = out[i] * sc + of;
        }
    }
#else
    (void)xoff; (void)yoff; (void)xsize; (void)ysize; (void)band; (void)out;
    throw std::runtime_error(
        "GdalReader::read_window: xbioclim was built without GDAL support."
    );
#endif
}

void GdalReader::read_bands_window(int xoff, int yoff,
                                   int xsize, int ysize,
                                   const std::vector<int>& bands,
                                   std::vector<double>& buf) const {
    const std::size_t n = static_cast<std::size_t>(xsize) * ysize;
    buf.resize(n * bands.size());
    read_bands_window(xoff, yoff, xsize, ysize, bands, buf.data());
}

void GdalReader::read_bands_window(int xoff, int yoff,
                                   int xsize, int ysize,
                                   const std::vector<int>& bands,
                                   double* out) const {
#ifdef HAVE_GDAL
    const std::size_t n = static_cast<std::size_t>(xsize) * ysize;

    CPLErr err = ds_->RasterIO(
        GF_Read,
        xoff, yoff, xsize, ysize,
        out,
        xsize, ysize,
        GDT_Float64,
        static_cast<int>(bands.size()),
        const_cast<int*>(bands.data()),
        static_cast<GIntBig>(sizeof(double)),
        0,
        static_cast<GIntBig>(n * sizeof(double)),
        nullptr
    );

    if (err != CE_None) {
        std::ostringstream oss;
        oss << "GdalReader::read_bands_window: RasterIO failed: "
            << CPLGetLastErrorMsg();
        throw std::runtime_error(oss.str());
    }

    // Apply scale / offset per band.
    for (std::size_t bi = 0; bi < bands.size(); ++bi) {
        GDALRasterBand* b = ds_->GetRasterBand(bands[bi]);
        if (b == nullptr) {
            std::ostringstream oss;
            oss << "GdalReader::read_bands_window: invalid band " << bands[bi];
            throw std::runtime_error(oss.str());
        }

        int has_scale  = 0;
        int has_offset = 0;
        double sc = b->GetScale(&has_scale);
        double of = b->GetOffset(&has_offset);

        if (has_scale || has_offset) {
            if (!has_scale)  sc = 1.0;
            if (!has_offset) of = 0.0;
            double* band_ptr = out + bi * n;
            for (std::size_t i = 0; i < n; ++i) {
                band_ptr[i] = band_ptr[i] * sc + of;
            }
        }
    }
#else
    (void)xoff; (void)yoff; (void)xsize; (void)ysize;
    (void)bands; (void)out;
    throw std::runtime_error(
        "GdalReader::read_bands_window: xbioclim was built without GDAL support."
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
                       bool cog_compatible,
                       GDALDataType dtype)
    : dtype_(dtype), closed_(false)
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
                         dtype, opts);
    CSLDestroy(opts);

    if (ds_ == nullptr) {
        std::ostringstream oss;
        oss << "GdalWriter: cannot create '" << path << "': "
            << CPLGetLastErrorMsg();
        throw std::runtime_error(oss.str());
    }

    // Set geotransform only if it is a valid, non-zero transform.
    if (geotransform.size() == 6 &&
        !std::all_of(geotransform.begin(), geotransform.end(),
                     [](double v) { return v == 0.0; })) {
        std::vector<double> gt_copy(geotransform);
        ds_->SetGeoTransform(gt_copy.data());
    }

    // Set CRS.
    if (!crs.empty()) {
        ds_->SetProjection(crs.c_str());
    }
#else
    (void)path; (void)nrows; (void)ncols; (void)nbands;
    (void)geotransform; (void)crs; (void)cog_compatible; (void)dtype;
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
    write_window(xoff, yoff, xsize, ysize, band, buf, dtype_);
}

void GdalWriter::write_window(int xoff, int yoff,
                              int xsize, int ysize,
                              int band,
                              const std::vector<double>& buf,
                              GDALDataType dtype) {
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

    CPLErr err;
    if (dtype == GDT_Float32) {
        std::vector<float> fbuf(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) {
            fbuf[static_cast<std::size_t>(i)] = static_cast<float>(buf[i]);
        }
        err = b->RasterIO(GF_Write,
                          xoff, yoff, xsize, ysize,
                          fbuf.data(),
                          xsize, ysize,
                          GDT_Float32,
                          0, 0);
    } else {
        err = b->RasterIO(GF_Write,
                          xoff, yoff, xsize, ysize,
                          const_cast<double*>(buf.data()),
                          xsize, ysize,
                          GDT_Float64,
                          0, 0);
    }

    if (err != CE_None) {
        std::ostringstream oss;
        oss << "GdalWriter::write_window: RasterIO failed for band " << band
            << ": " << CPLGetLastErrorMsg();
        throw std::runtime_error(oss.str());
    }
#else
    (void)xoff; (void)yoff; (void)xsize; (void)ysize; (void)band;
    (void)buf; (void)dtype;
    throw std::runtime_error(
        "GdalWriter::write_window: xbioclim was built without GDAL support."
    );
#endif
}

void GdalWriter::write_bands_window(int xoff, int yoff,
                                    int xsize, int ysize,
                                    const std::vector<int>& bands,
                                    const std::vector<double>& buf) {
    write_bands_window(xoff, yoff, xsize, ysize, bands, buf, dtype_);
}

void GdalWriter::write_bands_window(int xoff, int yoff,
                                    int xsize, int ysize,
                                    const std::vector<int>& bands,
                                    const std::vector<double>& buf,
                                    GDALDataType dtype) {
#ifdef HAVE_GDAL
    if (closed_) {
        throw std::runtime_error(
            "GdalWriter::write_bands_window: dataset is already closed.");
    }

    const std::size_t n = static_cast<std::size_t>(xsize) * ysize;
    const std::size_t needed = n * bands.size();
    if (buf.size() < needed) {
        std::ostringstream oss;
        oss << "GdalWriter::write_bands_window: buffer too small ("
            << buf.size() << " < " << needed << ")";
        throw std::runtime_error(oss.str());
    }

    CPLErr err;
    if (dtype == GDT_Float32) {
        std::vector<float> fbuf(buf.begin(), buf.begin() + needed);
        err = ds_->RasterIO(
            GF_Write,
            xoff, yoff, xsize, ysize,
            fbuf.data(),
            xsize, ysize,
            GDT_Float32,
            static_cast<int>(bands.size()),
            const_cast<int*>(bands.data()),
            static_cast<GIntBig>(sizeof(float)),
            0,
            static_cast<GIntBig>(n * sizeof(float)),
            nullptr
        );
    } else {
        err = ds_->RasterIO(
            GF_Write,
            xoff, yoff, xsize, ysize,
            const_cast<double*>(buf.data()),
            xsize, ysize,
            GDT_Float64,
            static_cast<int>(bands.size()),
            const_cast<int*>(bands.data()),
            static_cast<GIntBig>(sizeof(double)),
            0,
            static_cast<GIntBig>(n * sizeof(double)),
            nullptr
        );
    }

    if (err != CE_None) {
        std::ostringstream oss;
        oss << "GdalWriter::write_bands_window: RasterIO failed: "
            << CPLGetLastErrorMsg();
        throw std::runtime_error(oss.str());
    }
#else
    (void)xoff; (void)yoff; (void)xsize; (void)ysize;
    (void)bands; (void)buf; (void)dtype;
    throw std::runtime_error(
        "GdalWriter::write_bands_window: xbioclim was built without GDAL support."
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
