# GDAL Tiled I/O Layer

This document describes the GDAL tiled I/O integration added in **Phase 2 /
Issue #22** of the rxbioclim roadmap.

---

## Overview

`rxbioclim` can optionally use GDAL for reading and writing raster files
without loading entire datasets into R memory.  When GDAL is detected at build
time, the package gains two diagnostic R functions and a C++ API for tiled
processing.

The package **always compiles** even when GDAL is absent (Windows CI, minimal
Linux environments, etc.).  GDAL-backed features fail at runtime with a clear
error message in that case.

---

## Build-time detection

`configure.ac` (added in PR #35) searches for `gdal-config` and runs a
link test.  If GDAL ≥ 2.0.1 is found, `src/Makevars` is generated with:

```
PKG_CPPFLAGS = -DHAVE_GDAL <gdal-config --cflags>
PKG_LIBS     = <gdal-config --libs>
```

All GDAL-dependent C++ code is guarded by `#ifdef HAVE_GDAL`.

---

## C++ API

### `rxbioclim::GdalReader` (`src/gdal_io.hpp` / `src/gdal_io.cpp`)

Opens a raster for **read-only** tiled access.

```cpp
rxbioclim::GdalReader reader("/path/to/raster.tif");

int rows  = reader.nrows();
int cols  = reader.ncols();
int bands = reader.nbands();

std::vector<double> gt = reader.geotransform();  // 6-element GDAL array
std::string wkt        = reader.crs();

// Per-band scale / offset (returns 1.0 / 0.0 if not set)
double sc = reader.scale(1);   // band 1 (1-based)
double of = reader.offset(1);

// Read a 256×256 tile from band 1 (values always float64;
// scale/offset applied automatically when present)
std::vector<double> buf;
reader.read_window(0, 0, 256, 256, 1, buf);
```

Only the pixels in the requested window are read into memory; the rest of
the raster stays on disk.

### `rxbioclim::GdalWriter` (`src/gdal_io.hpp` / `src/gdal_io.cpp`)

Creates a new GTiff for **write-only** tiled output.

```cpp
std::vector<double> gt = {-180.0, 0.5, 0.0, 90.0, 0.0, -0.5};
rxbioclim::GdalWriter writer(
    "/path/to/output.tif",
    nrows, ncols, nbands,
    gt,               // 6-element geotransform (may be empty)
    wkt_crs,          // WKT CRS string (may be empty)
    /*cog_compatible=*/true   // adds TILED + LZW + BIGTIFF=IF_SAFER
);

std::vector<double> tile(256 * 256, 0.0);
writer.write_window(0, 0, 256, 256, 1, tile);  // write first tile, band 1

writer.close();  // flush and finalize
```

Setting `cog_compatible = true` adds these GTiff creation options:

| Option         | Value      |
|----------------|------------|
| `TILED`        | `YES`      |
| `BLOCKXSIZE`   | `256`      |
| `BLOCKYSIZE`   | `256`      |
| `COMPRESS`     | `LZW`      |
| `BIGTIFF`      | `IF_SAFER` |

> **Note:** Full COG validation (e.g. `gdaladdo` overviews) is out of scope
> for this phase.  The options above produce a file that satisfies the tiling
> and compression requirements for COG-compatible workflows.

---

## R Interface

Two diagnostic functions are exported to R:

### `gdal_can_open(path)`

Returns `TRUE` if GDAL can open the file at `path`, `FALSE` otherwise.
Stops with `"built without GDAL"` if the package was compiled without GDAL.

```r
gdal_can_open(system.file("extdata", "tiny.tif", package = "rxbioclim"))
#> [1] TRUE
```

### `gdal_info(path)`

Returns a named list of metadata for the raster at `path`:

| Element        | Type      | Description                                    |
|----------------|-----------|------------------------------------------------|
| `path`         | character | Input path                                     |
| `nrows`        | integer   | Number of rows                                 |
| `ncols`        | integer   | Number of columns                              |
| `nbands`       | integer   | Number of bands                                |
| `geotransform` | numeric[6]| GDAL geotransform array                        |
| `crs`          | character | WKT CRS (empty string if undefined)            |
| `scale`        | numeric[] | Per-band scale factor (1.0 if not set)         |
| `offset`       | numeric[] | Per-band offset (0.0 if not set)               |

---

## Scale / Offset Handling

Many WorldClim rasters are stored as packed `Int16` values with a scale
factor (e.g. `0.1` for temperature in 0.1 °C units).  `GdalReader` detects
these via the GDAL band metadata and applies the transformation:

```
physical_value = stored_value × scale + offset
```

The returned buffer always contains `double` (`float64`) physical values.

---

## Test Raster

A minimal 3×3 Float32 GeoTIFF is provided at `inst/extdata/tiny.tif` for
smoke testing.  Pixel values are `1.0` to `9.0` (row-major), with no
geotransform or CRS set (defaults to zeros / empty string).

---

## Thread Safety

Each `GdalReader` / `GdalWriter` instance owns its `GDALDataset` exclusively.
Do **not** share a single object across threads; open one reader/writer per
thread instead.
