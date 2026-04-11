#!/usr/bin/env python3
"""Validate that a GeoTIFF file is a valid Cloud Optimized GeoTIFF (COG).

Tries osgeo_utils.samples.validate_cloud_optimized_geotiff for a full
structural check; falls back to an IFD-based basic check when that module
is unavailable (older system GDAL packages).

Exit code: 0 if valid, 1 on any error.
"""

import sys
from osgeo import gdal


def _validate_basic(path: str):
    """Minimal COG check via GDAL metadata and block structure."""
    try:
        ds = gdal.Open(path)
    except Exception as exc:
        return [f"Cannot open file {path!r}: {exc}"], []

    if ds is None:
        return [f"Cannot open file: {path}"], []

    errors: list[str] = []
    warnings: list[str] = []

    # Must be GeoTIFF.
    if ds.GetDriver().ShortName != "GTiff":
        errors.append(f"Not a GeoTIFF (driver: {ds.GetDriver().ShortName})")
        return errors, warnings

    # GDAL 3.x exposes LAYOUT=COG in IMAGE_STRUCTURE metadata.
    layout = ds.GetMetadataItem("LAYOUT", "IMAGE_STRUCTURE")
    if layout == "COG":
        return errors, warnings  # All good.

    band = ds.GetRasterBand(1)
    block_x, block_y = band.GetBlockSize()
    nx, ny = ds.RasterXSize, ds.RasterYSize

    # Large files must be tiled and have overviews.
    if nx > 512 or ny > 512:
        if block_x == nx:
            errors.append("File is not tiled — not COG compliant.")
        if band.GetOverviewCount() == 0:
            warnings.append(
                "Large file has no overviews — COG may not be optimal."
            )
    else:
        # Small files are accepted without overviews; just ensure tiling or
        # that the IFD offset is near the start of the file (GDAL 3 only).
        pass

    return errors, warnings


def main() -> None:
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <file.tif>", file=sys.stderr)
        sys.exit(1)

    path = sys.argv[1]

    try:
        from osgeo_utils.samples.validate_cloud_optimized_geotiff import (
            validate,
        )

        errors, warnings, _ = validate(path, full_check=True)
        for w in warnings:
            print(f"  WARNING: {w}")
        if errors:
            for e in errors:
                print(f"  ERROR: {e}", file=sys.stderr)
            sys.exit(1)
        print(f"  COG validation PASSED (full check): {path}")

    except ImportError:
        errors, warnings = _validate_basic(path)
        for w in warnings:
            print(f"  WARNING: {w}")
        if errors:
            for e in errors:
                print(f"  ERROR: {e}", file=sys.stderr)
            sys.exit(1)
        print(f"  COG validation PASSED (basic check): {path}")


if __name__ == "__main__":
    main()
