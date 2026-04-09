#!/usr/bin/env python3
"""Generate synthetic 10x10 GeoTIFF test fixtures for xbioclim.

Creates 48 compressed Int16 GeoTIFFs (4 variables x 12 months) in --outdir.
All pixels within a layer share the same uniform integer value as defined in
docs/ROADMAP.md §1.1:

  tas[m]    = m          (m = 1..12)
  tasmax[m] = m + 1
  tasmin[m] = m - 1
  pr[m]     = m

Output filename convention: <var>_<MM>.tif  (e.g. tas_01.tif, pr_12.tif)
"""

import argparse
import pathlib

import numpy as np
from osgeo import gdal, osr

# Geotransform for a 10-pixel x 10-pixel global grid (10-degree cells)
GT = (-180.0, 36.0, 0.0, 90.0, 0.0, -18.0)

VARS = {
    "tas":    lambda m: m,
    "tasmax": lambda m: m + 1,
    "tasmin": lambda m: m - 1,
    "pr":     lambda m: m,
}


def write_tif(path: pathlib.Path, value: int, nodata: int = -9999) -> None:
    """Write a 10x10 uniform Int16 GeoTIFF with LZW compression."""
    drv = gdal.GetDriverByName("GTiff")
    ds = drv.Create(
        str(path),
        10, 10, 1,
        gdal.GDT_Int16,
        ["COMPRESS=LZW"],
    )
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    ds.SetProjection(srs.ExportToWkt())
    ds.SetGeoTransform(GT)
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    band.WriteArray(np.full((10, 10), value, dtype=np.int16))
    ds.FlushCache()
    ds = None  # close


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Generate synthetic GeoTIFF fixtures for xbioclim tests."
    )
    ap.add_argument(
        "--outdir",
        default="tests/data",
        help="Output directory (default: tests/data)",
    )
    args = ap.parse_args()

    out = pathlib.Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    count = 0
    for var, fn in VARS.items():
        for m in range(1, 13):
            filename = out / f"{var}_{m:02d}.tif"
            write_tif(filename, fn(m))
            count += 1

    print(f"Generated {count} GeoTIFFs in {out}/")


if __name__ == "__main__":
    main()
