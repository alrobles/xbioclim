// gdal_mask.cpp — implementation of polygon masking for xbioclim.
//
// All GDAL-dependent code is inside #ifdef HAVE_GDAL blocks.  When GDAL is
// absent the functions throw std::runtime_error with a clear message.

#include "gdal_mask.hpp"
#include "gdal_io.hpp"

#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef HAVE_GDAL
#include <gdal_priv.h>
#include <gdal_alg.h>   // GDALRasterizeLayers
#include <ogr_api.h>
#include <ogrsf_frmts.h>
#include <cpl_conv.h>
#endif

namespace xbioclim {

// ============================================================================
// rasterize_mask
// ============================================================================

void rasterize_mask(const std::string& vector_path,
                    const std::string& ref_raster,
                    const std::string& output_path)
{
#ifdef HAVE_GDAL
    GDALAllRegister();

    // ── 1. Open reference raster to obtain spatial metadata ─────────────────
    GdalReader ref(ref_raster);
    const int nrows = ref.nrows();
    const int ncols = ref.ncols();
    std::vector<double> gt = ref.geotransform();
    const std::string   wkt = ref.crs();

    // ── 2. Open vector source ────────────────────────────────────────────────
    GDALDataset* vec_ds = static_cast<GDALDataset*>(
        GDALOpenEx(vector_path.c_str(),
                   GDAL_OF_VECTOR | GDAL_OF_READONLY,
                   nullptr, nullptr, nullptr));
    if (vec_ds == nullptr) {
        std::ostringstream oss;
        oss << "rasterize_mask: cannot open vector source '"
            << vector_path << "': " << CPLGetLastErrorMsg();
        throw std::runtime_error(oss.str());
    }

    // ── 3. Create output mask raster (GDT_Byte, single band, zeros) ─────────
    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("GTiff");
    if (drv == nullptr) {
        GDALClose(vec_ds);
        throw std::runtime_error(
            "rasterize_mask: GTiff driver not available.");
    }

    GDALDataset* out_ds = drv->Create(output_path.c_str(),
                                      ncols, nrows, 1,
                                      GDT_Byte, nullptr);
    if (out_ds == nullptr) {
        GDALClose(vec_ds);
        std::ostringstream oss;
        oss << "rasterize_mask: cannot create output raster '"
            << output_path << "': " << CPLGetLastErrorMsg();
        throw std::runtime_error(oss.str());
    }

    if (gt.size() == 6) {
        out_ds->SetGeoTransform(gt.data());
    }
    if (!wkt.empty()) {
        out_ds->SetProjection(wkt.c_str());
    }

    // Initialise to 0 (outside all polygons).
    GDALRasterBand* band = out_ds->GetRasterBand(1);
    band->Fill(0.0);

    // ── 4. Collect layers and burn value = 1 ────────────────────────────────
    const int n_layers = vec_ds->GetLayerCount();
    std::vector<OGRLayerH> layers;
    layers.reserve(static_cast<std::size_t>(n_layers));
    for (int i = 0; i < n_layers; ++i) {
        OGRLayer* lyr = vec_ds->GetLayer(i);
        if (lyr != nullptr) {
            layers.push_back(static_cast<OGRLayerH>(lyr));
        }
    }

    if (layers.empty()) {
        GDALClose(out_ds);
        GDALClose(vec_ds);
        throw std::runtime_error(
            "rasterize_mask: vector source contains no layers.");
    }

    // Burn value 1 for all polygon features in every layer.
    std::vector<double> burn_values(layers.size(), 1.0);

    int bands_idx[1] = { 1 };  // 1-based band index for GDALRasterizeLayers

    CPLErr err = GDALRasterizeLayers(
        out_ds,
        1,                              // nBands
        bands_idx,                      // panBandList (1-based band indices)
        static_cast<int>(layers.size()),
        layers.data(),
        nullptr,                        // pfnTransformer (identity)
        nullptr,                        // pTransformArg
        burn_values.data(),             // padfLayerBurnValues
        nullptr,                        // papszOptions
        nullptr,                        // pfnProgress
        nullptr                         // pProgressArg
    );

    GDALClose(out_ds);
    GDALClose(vec_ds);

    if (err != CE_None) {
        std::ostringstream oss;
        oss << "rasterize_mask: GDALRasterizeLayers failed: "
            << CPLGetLastErrorMsg();
        throw std::runtime_error(oss.str());
    }
#else
    (void)vector_path; (void)ref_raster; (void)output_path;
    throw std::runtime_error(
        "rasterize_mask: xbioclim was built without GDAL support. "
        "Install GDAL >= 2.0.1 and reinstall the package.");
#endif
}

// ============================================================================
// apply_mask
// ============================================================================

void apply_mask(const std::string& input_path,
                const std::string& mask_path,
                const std::string& output_path)
{
#ifdef HAVE_GDAL
    GDALAllRegister();

    GdalReader  inp(input_path);
    GdalReader  msk(mask_path);

    const int nrows  = inp.nrows();
    const int ncols  = inp.ncols();
    const int nbands = inp.nbands();

    // Validate mask dimensions match input.
    if (msk.nrows() != nrows || msk.ncols() != ncols) {
        std::ostringstream oss;
        oss << "apply_mask: mask dimensions ("
            << msk.nrows() << "x" << msk.ncols()
            << ") do not match input dimensions ("
            << nrows << "x" << ncols << ").";
        throw std::runtime_error(oss.str());
    }

    GdalWriter out(output_path,
                   nrows, ncols, nbands,
                   inp.geotransform(),
                   inp.crs());

    const double nan_val = std::numeric_limits<double>::quiet_NaN();

    // Process row by row to keep memory usage constant.
    std::vector<double> inp_buf;
    std::vector<double> msk_buf;

    for (int row = 0; row < nrows; ++row) {
        msk.read_window(0, row, ncols, 1, 1, msk_buf);

        for (int b = 1; b <= nbands; ++b) {
            inp.read_window(0, row, ncols, 1, b, inp_buf);

            for (int col = 0; col < ncols; ++col) {
                if (msk_buf[static_cast<std::size_t>(col)] == 0.0) {
                    inp_buf[static_cast<std::size_t>(col)] = nan_val;
                }
            }

            out.write_window(0, row, ncols, 1, b, inp_buf);
        }
    }

    out.close();
#else
    (void)input_path; (void)mask_path; (void)output_path;
    throw std::runtime_error(
        "apply_mask: xbioclim was built without GDAL support. "
        "Install GDAL >= 2.0.1 and reinstall the package.");
#endif
}

} // namespace xbioclim
