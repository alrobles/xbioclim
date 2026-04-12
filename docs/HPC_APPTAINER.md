# Running xbioclim on Slurm HPC Clusters with Apptainer

## Overview

This guide explains how to build and run **xbioclim** on High-Performance
Computing (HPC) clusters using [Apptainer](https://apptainer.org/) (formerly
Singularity) containers. Apptainer is the standard container runtime in HPC
environments because it:

- **Requires no root privileges at runtime** — jobs run as your regular user.
- **Produces portable, single-file images** (`.sif`) — copy one file to any
  cluster and run it.
- **Ensures reproducibility** — every dependency version is locked inside the
  image.

The container bundles all of xbioclim's dependencies (GDAL, xtensor, xtl,
OpenMP) so you do not need to install or compile anything on the cluster itself.

---

## Prerequisites

| Requirement | Details |
|-------------|---------|
| **Apptainer** ≥ 1.1 (or Singularity ≥ 3.8) | Installed on the build machine and available on the cluster |
| **Root or `--fakeroot`** | Needed **only** during image build, not at runtime |
| **Slurm cluster access** | For submitting batch jobs |
| **Input raster files** | Monthly GeoTIFF files in the layout `${INPUT_DIR}/${VAR}_${MM}.tif` (e.g. `tas_01.tif` … `tas_12.tif`) for variables: `tas`, `tasmax`, `tasmin`, `pr` |

---

## Building the Apptainer Image

### 1. Clone the repository

```bash
git clone https://github.com/alrobles/xbioclim.git
cd xbioclim
```

The definition file is located at `apptainer/xbioclim.def`.

### 2. Build the image

**With root access:**

```bash
sudo apptainer build xbioclim.sif apptainer/xbioclim.def
```

**With `--fakeroot` (no root required):**

```bash
apptainer build --fakeroot xbioclim.sif apptainer/xbioclim.def
```

> **Note:** The build downloads and compiles several packages. Expect it to take
> approximately 10–20 minutes depending on your internet speed and CPU.

### 3. Transfer the image to the cluster

```bash
scp xbioclim.sif youruser@cluster.example.edu:~/xbioclim.sif
# or
rsync -avP xbioclim.sif youruser@cluster.example.edu:~/
```

---

## GPU / CUDA Variant

To build an image with CUDA support for GPU-accelerated computation:

1. Open `apptainer/xbioclim.def`.
2. Change the header to use the NVIDIA CUDA base image:
   ```
   Bootstrap: docker
   From: nvidia/cuda:12.3.2-devel-ubuntu22.04
   ```
3. Uncomment the GPU variant `%post` section (or add `-DXBIOCLIM_USE_CUDA=ON`
   to the cmake configure line).
4. Build the image as described above.

**Hardware requirements:** NVIDIA GPU with Compute Capability ≥ 8.0
(A100, H100, RTX 3090, etc.)

When running the CUDA-enabled image, add the `--nv` flag to grant GPU access:

```bash
apptainer exec --nv --bind /data:/data xbioclim.sif xbioclim_cli ...
```

---

## Running Interactively

### Inspect the container

```bash
apptainer shell xbioclim.sif
```

### Test the CLI

```bash
apptainer run xbioclim.sif --help
```

### Process rasters with a bind mount

```bash
apptainer exec \
    --bind /data:/data \
    xbioclim.sif \
    xbioclim_cli \
        --tas /data/tas_{01..12}.tif \
        --tasmax /data/tasmax_{01..12}.tif \
        --tasmin /data/tasmin_{01..12}.tif \
        --pr /data/pr_{01..12}.tif \
        --outdir /data/output \
        --prefix bio \
        --tile-size 512
```

---

## Submitting Slurm Array Jobs

The repository includes a ready-to-use Slurm script at
`scripts/run_slurm_apptainer.sh`. It splits the global raster into **360
longitude strips** (one per array task, from −180° to +179°).

### 1. Set the required environment variables

```bash
export INPUT_DIR=/scratch/myproject/rasters
export OUTPUT_DIR=/scratch/myproject/output
export SIF=$HOME/xbioclim.sif
```

### 2. Create the logs directory and submit

```bash
mkdir -p logs
sbatch scripts/run_slurm_apptainer.sh
```

### 3. Monitor the jobs

```bash
# List running/pending jobs
squeue -u $USER

# Show detailed accounting after completion
sacct -j <JOBID> --format=JobID,Elapsed,MaxRSS,State
```

---

## Bind Mounts

Apptainer containers have their own isolated filesystem. To let the container
access host directories (input data, output paths), you must **bind-mount**
them with `--bind`:

```bash
apptainer exec --bind /scratch:/scratch xbioclim.sif xbioclim_cli ...
```

### Common patterns

| Scenario | Flag |
|----------|------|
| Bind scratch | `--bind /scratch:/scratch` |
| Bind project space | `--bind /projects:/projects` |
| Bind multiple paths | `--bind /scratch:/scratch --bind /data:/data` |

### Auto-binding

Some clusters automatically bind `/home`, `/scratch`, and `/tmp` into every
container. Check your cluster documentation or run:

```bash
apptainer exec xbioclim.sif ls /scratch
```

> **Warning:** Do not rely on auto-binding for portability. Always use explicit
> `--bind` flags in Slurm scripts so the job works on any cluster.

---

## Troubleshooting

| Problem | Cause / Solution |
|---------|------------------|
| `FATAL: container creation failed` | You need root or `--fakeroot` privileges to **build** the image. Building is not needed at runtime. |
| `xbioclim_cli: command not found` | The image may not have built correctly. Verify with: `apptainer exec xbioclim.sif which xbioclim_cli` |
| GDAL errors about missing drivers | Ensure input files are valid GeoTIFF. Check with: `apptainer exec xbioclim.sif gdalinfo input.tif` |
| Slow OpenMP performance | Set `OMP_NUM_THREADS` to match `--cpus-per-task`. E.g. add `export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}` in your script. |
| `SIF file not found` | Verify the `$SIF` path and that the `.sif` file was transferred to the cluster correctly. |
| Permission denied on output directory | Ensure `OUTPUT_DIR` exists and is writable before submitting the job. |

---

## Environment Variables Reference

| Variable | Default | Description |
|----------|---------|-------------|
| `INPUT_DIR` | *(required)* | Directory containing the monthly raster files (`tas_01.tif`, etc.) |
| `OUTPUT_DIR` | *(required)* | Root output directory; per-strip subdirectories are created automatically |
| `SIF` | `$HOME/xbioclim.sif` | Path to the Apptainer `.sif` image |
| `TILE_SIZE` | `512` | Processing tile size in pixels (controls memory usage) |
| `OMP_NUM_THREADS` | *(system default)* | Number of OpenMP threads; set to `$SLURM_CPUS_PER_TASK` for best results |

---

## Further Reading

- [Apptainer Documentation](https://apptainer.org/docs/)
- [Slurm Workload Manager](https://slurm.schedmd.com/documentation.html)
- [xbioclim README](../README.md)
- [GDAL Documentation](https://gdal.org/en/stable/index.html)
- [xtensor Documentation](https://xtensor.readthedocs.io/)
