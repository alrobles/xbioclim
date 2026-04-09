#pragma once

#include "xbioclim/primitives.hpp"

#include <gdal_priv.h>

#include <array>
#include <string>
#include <vector>

namespace xbioclim {

/// Reads monthly climate rasters from 48 GeoTIFFs (4 variables × 12 months).
class GdalReader {
public:
    explicit GdalReader(const std::vector<std::string>& tas_files,
                        const std::vector<std::string>& tasmax_files,
                        const std::vector<std::string>& tasmin_files,
                        const std::vector<std::string>& pr_files);

    ~GdalReader();

    /// Read a tile into a ClimateBlock (pixel values as float32).
    ClimateBlock read_block(int x_off, int y_off,
                            int block_w, int block_h) const;

    int raster_x_size() const;
    int raster_y_size() const;
    std::string projection_wkt() const;
    std::array<double, 6> geotransform() const;

private:
    struct Impl;
    Impl* impl_;
};

/// Writes 19 BIO variable rasters to a set of output GeoTIFFs.
class GdalWriter {
public:
    GdalWriter(const std::string& output_dir,
               const std::string& prefix,
               int x_size, int y_size,
               const std::string& projection_wkt,
               const std::array<double, 6>& geotransform);

    ~GdalWriter();

    /// Write a BioBlock tile to the output GeoTIFFs.
    void write_block(int x_off, int y_off, const BioBlock& bio);

    /// Flush and close all open datasets.
    void finalize();

private:
    struct Impl;
    Impl* impl_;
};

} // namespace xbioclim
