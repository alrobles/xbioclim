# xbioclim — Release Guidelines

This document describes the release process for xbioclim, including the
mandatory COG (Cloud Optimized GeoTIFF) output validation that every release
artifact must pass before the GitHub Release is published.

---

## Release Checklist

Before tagging a release:

- [ ] All 19 BIO variables pass analytical tests (ε < 1 × 10⁻³).
- [ ] CI is green on Ubuntu 22.04 for both GCC and Clang builds.
- [ ] `clang-tidy` and `cppcheck` pass with zero warnings.
- [ ] Doxygen API documentation builds without errors.
- [ ] **COG output validated** — every output `.tif` passes
      `gdalinfo --checksum` and the structural COG check (see below).
- [ ] Performance benchmarks are documented in `docs/BENCHMARKS.md`.
- [ ] The `CHANGELOG` / GitHub release notes are up to date.
- [ ] Semantic versioning tag (`vMAJOR.MINOR.PATCH`) is created on `main`.

---

## COG Output Validation

### What is a Cloud Optimized GeoTIFF?

A **Cloud Optimized GeoTIFF (COG)** is a regular GeoTIFF file with an
internal layout that supports efficient range-request access (e.g. from
object storage such as S3 or GCS). Key structural requirements:

1. Image File Directory (IFD) is near the beginning of the file.
2. Tiles are stored in a predictable top-left–to–bottom-right order.
3. Overview levels (reduced-resolution copies) are included for large files,
   and overview tiles are stored *before* full-resolution tiles.

xbioclim's `GdalWriter` uses the GDAL **COG driver** (available since
GDAL 3.1) to write output files that satisfy all requirements automatically.

### Running the Validation Locally

After running the CLI:

```bash
# Generate outputs (replace paths with your real input files)
./build/xbioclim_cli \
  --tas    tests/data/tas_??.tif \
  --tasmax tests/data/tasmax_??.tif \
  --tasmin tests/data/tasmin_??.tif \
  --pr     tests/data/pr_??.tif \
  --outdir output --prefix bio --tile-size 10

# Validate every output file
bash scripts/validate_cog_outputs.sh output
```

The script performs two checks per file:

| Check | Tool | Failure condition |
|-------|------|-------------------|
| Checksum integrity | `gdalinfo -checksum` | non-zero exit code |
| COG structure | `scripts/check_cog.py` | any structural error |

The script exits with code **1** if any file fails; CI jobs that call it will
therefore block the subsequent release step.

### Validating a Single File

```bash
# Checksum
gdalinfo -checksum output/bio1.tif

# COG structure
python3 scripts/check_cog.py output/bio1.tif
```

### What to Do When Validation Fails

1. **`gdalinfo --checksum` fails** — the file is likely corrupt or empty.
   Re-run the pipeline from scratch; if the failure persists, check for
   disk-space or permission issues.

2. **COG structure errors** (e.g. "IFD not at the start of the file") —
   the output was probably written by an older build that did not use the
   COG driver. Rebuild from the current `main` branch and re-run.

3. **Missing overviews warning for large files** — add overview levels with:
   ```bash
   gdaladdo -r average output/bio1.tif 2 4 8 16
   ```
   Then re-validate.

---

## CI Integration

The validation runs automatically in two places:

### `pipeline.yml` — `performance-benchmarks` job

After a full end-to-end CLI run on every push to `main`/`develop` and on
pull requests, the job calls:

```yaml
- name: Validate COG outputs (checksum + structure)
  run: bash scripts/validate_cog_outputs.sh output
```

Failure here blocks the downstream `release-integration` job.

### `pipeline.yml` — `release-integration` job (tags only)

Before publishing artifacts to the GitHub Release, the job re-validates
the downloaded performance results:

```yaml
- name: Validate release artifacts (COG + checksum)
  run: bash scripts/validate_cog_outputs.sh performance
```

### `release.yml` — Release workflow (tags only)

The end-to-end smoke test runs the CLI and then validates outputs with the
same script before any artifacts are uploaded or the GitHub Release is
created:

```yaml
- name: End-to-end CLI smoke test + COG validation
  run: |
    ...
    bash scripts/validate_cog_outputs.sh e2e_out
```

If this step fails the release is **blocked** — no artifacts are uploaded
and no GitHub Release is created.

---

## Tagging a Release

```bash
git tag -a v0.3.0 -m "Release v0.3.0"
git push origin v0.3.0
```

Pushing the tag triggers both `release.yml` and the tag path of
`pipeline.yml`. Both workflows must be green before the release is
considered valid.
