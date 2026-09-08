# Phase C: Overlapped I/O / Compute / Write Pipeline

## 1. Overview and bottleneck

`BioclimEngine::compute()` (on `feat/phase-a-multiband`, commit `10ca0d1`) processes the output raster tile by tile in a single serial loop:

```
for each tile:
    READ   → multi-band RasterIO of the 4 climate variables
    COMPUTE→ OpenMP per-pixel BIOCLIM loop
    WRITE  → multi-band RasterIO of the 19 output bands
```

The wall time of every tile is therefore approximately

```
tile_time = read_time + compute_time + write_time
```

Profiling shows that the read and compute stages can each be a large fraction of the per-tile time, especially with multi-band input files (12 bands per variable) and large tile sizes. The three stages are essentially independent for a given tile (the inputs, intermediate buffers, and output window are distinct), so they can be overlapped. The ideal limit is

```
tile_time ≈ max(read_time, compute_time, write_time)
```

provided the pipeline can keep one read, one compute, and one write in flight at the same time.

## 2. Proposed pipeline stages

The pipeline keeps the same logical stages as the current engine, but runs them concurrently across tiles:

- **READ stage**: For one `(xoff, yoff, xsize, ysize)` window, read the 12 monthly bands of `tas`, `tasmax`, `tasmin`, and `pr` into the tile buffers, plus the optional mask band. Inputs are already multi-band, so each variable is read with a single `GDALDataset::RasterIO` call when it is a multi-band file, or 12 single-band calls when it is a stack of 12 files.
- **COMPUTE stage**: Run the per-pixel OpenMP `compute_pixel` kernel on the filled tile buffers. Output is the 19-band `bio` tile, still in band-major layout `bio[bio_index * n_pix + pixel]` to match the existing `write_bands_window` call.
- **WRITE stage**: Call `GdalWriter::write_bands_window` for the 19 output bands of that tile.

The ordering of tiles is fixed: the engine processes the raster in row-major order (`yoff` outer, `xoff` inner), and writes must occur in that order if any driver or downstream consumer is sensitive to sequential access. The queue-based design below naturally preserves this order.

## 3. Buffer and ownership design: three-slot ring (generalisation of double-buffering)

A classic *front/back* double buffer would let the reader fill one buffer while the compute thread consumes the other. To also overlap writes, we extend the scheme to a bounded ring of **three reusable tile buffers**:

| Slot state | Owner        | Activity                          |
|------------|--------------|-----------------------------------|
| `free`     | free queue   | waiting to be filled              |
| `read`     | reader thread| filling with GDAL input           |
| `compute`  | compute thread| running `compute_pixel`           |
| `write`    | writer thread| writing via `GdalWriter`          |

At peak, all three slots can be live simultaneously: slot *N-1* being written, slot *N* being computed, and slot *N+1* being read.

Data structures:

```cpp
struct TileSlot {
    int xoff = 0, yoff = 0, xsize = 0, ysize = 0, n_pix = 0;
    std::vector<double> tas;      // 12 * n_pix
    std::vector<double> tasmax;   // 12 * n_pix
    std::vector<double> tasmin;   // 12 * n_pix
    std::vector<double> pr;       // 12 * n_pix
    std::vector<double> mask;     // n_pix, or empty
    std::vector<double> bio;      // 19 * n_pix
};
```

The buffers are resized per tile to `n_pix = xsize * ysize`; all slots are pre-allocated once to `tile_size_ * tile_size_` so partial edge tiles do not cause repeated heap allocations.

Three `std::condition_variable` / `std::mutex` queues implement ownership hand-off:

1. `free_slots`  → reader pops, fills, pushes to `ready_for_compute`.
2. `ready_for_compute` → compute pops, processes, pushes to `ready_for_write`.
3. `ready_for_write` → writer pops, writes, pushes back to `free_slots`.

A `nullptr` sentinel is pushed at the end of `ready_for_compute` to signal end-of-stream; the compute thread forwards the sentinel to `ready_for_write`, and the writer exits after seeing it.

## 4. Threading model

**Choice: three dedicated threads (reader, compute, writer).**

Alternatives considered:

- **One I/O thread + one compute thread, with compute doing the write**: Simpler, but the write is not overlapped with the next read. GDAL writes can be non-trivial (e.g. BigTIFF, compressed output, or network-backed VSI), so leaving write in the compute stage leaves speedup on the table.
- **Thread-pool / task graph (OpenMP tasks or a custom pool)**: More general, but adds scheduling overhead and makes GDAL dataset ownership harder to reason about. A tile is large enough that the cost of three `std::thread`s is negligible for the whole raster, and the strict one-thread-per-dataset rule is trivial to enforce.

The reader and writer are I/O-bound; the compute thread is CPU-bound and may itself fan out with OpenMP. Using dedicated threads means:

- Exactly one thread touches each `GdalReader`.
- Exactly one thread touches the `GdalWriter`.
- The main thread waits for the pipeline to finish and then closes the writer.

This is the safest model for a GDAL prototype, because most GDAL drivers do **not** support concurrent access to the same `GDALDataset*`.

## 5. GDAL concurrency safety

The serial code already opens a separate `GDALDataset` for each input variable and for the output. The proposed pipeline preserves the same object ownership:

- All `GdalReader` objects live only in the reader thread.
- The `GdalWriter` object lives only in the writer thread.
- The three threads never share a `GDALDataset*`, with the possible exception of the 12 single-band readers belonging to the same variable (they are still one `GdalReader`/dataset per file, all used by the same reader thread).

Because reads and writes are on different datasets, they can proceed in parallel. If in the future a driver is proven thread-safe, the design can be extended to allow parallel reads of different variables, but the prototype does not rely on that.

## 6. Synchronisation

Two options were considered for the queues:

1. **Mutex + condition variable** (selected): simple, correct, and the tile-granularity lock contention is negligible because each critical section only moves a pointer and calls `notify_one`. The `std::mutex` protects the `std::deque<TileSlot*>` and a `cancelled` flag.
2. **Lock-free SPSC/MPSC queue** (not selected): would reduce latency slightly, but the queues are between exactly one producer and one consumer, and the simpler mutex/CV implementation is easier to reason about and to extend with cancellation on error.

Error handling: each worker wraps its loop in `try/catch`, stores the first exception message in a shared `std::string` under a mutex, calls `cancel()` on all queues, and exits. `cancel()` causes blocked `pop()` calls to return `nullptr`, so the remaining threads drain any queued real tiles and then exit without deadlocking.

## 7. Edge cases

| Case | Handling |
|------|----------|
| **Partial tiles at raster edges** | `xsize`/`ysize` are `min(tile_size_, ncols/nrows - xoff/yoff)`. The slot buffers are resized to the actual `n_pix`; reads and writes use the exact window. No special compute logic is needed because the per-pixel loop only touches `0..n_pix-1`. |
| **Mask handling** | If a mask is configured, the reader thread reads the mask into `slot->mask` before pushing the slot to compute. The compute thread applies the same `isnan(mask) || mask == 0` test as the serial loop, writing NaN for masked pixels. The mask reader lives in the reader thread. |
| **Variable subset** | The engine always computes and writes all 19 bioclimatic variables to the 19-band `bio.tif`. Variable subsetting is still done by the R wrapper after `compute()` returns. The pipeline does not need to know the selected variables. |
| **Tile ordering** | Tiles are produced by the reader in the original row-major order and consumed/written in the same order because the queues are FIFO. If the writer were allowed to reorder, it would require an explicit reorder buffer; the prototype avoids this. |
| **OpenMP inside a worker thread** | The compute thread calls the existing OpenMP `parallel for` with `num_threads(n_threads_)`. This is safe because the loop is pure C++/GDAL and makes no R API calls. |
| **Empty / single-tile rasters** | The pipeline works for any number of tiles. With fewer than three tiles some queues simply go empty; the main thread joins the workers and closes the writer. |

## 8. Expected speedup

The theoretical best-case per-tile wall time is the longest stage:

```
tile_time ≈ max(read_time, compute_time, write_time)
```

For tiles where the three stages are comparable, the overlap can approach a 3× reduction in per-tile latency. More realistically, for current climate rasters:

- `read_time` is usually the largest (4 variables × 12 bands, or 4 multi-band RasterIO calls).
- `compute_time` is smaller but non-negligible.
- `write_time` is typically the smallest (one multi-band RasterIO).

Hence the prototype should primarily reduce the `read + compute` serial sum to closer to `max(read, compute)`, with the write hiding underneath. The overall wall-time speedup will depend on tile size, file format, and whether the input is multi-band or single-band. For small synthetic tests the effect may be modest; for large GeoTIFFs with tiled/compressed I/O it should be more pronounced.

## 9. Prototype implementation notes

The prototype adds:

- A `bool pipeline_` member and `BioclimEngine::set_pipeline(bool)` method.
- A `compute_pipelined()` private method in `src/BioclimEngine.cpp`.
- A `TileSlot` struct and `TileQueue` class in an anonymous namespace.
- A `engine_set_pipeline` Rcpp export with `@keywords internal` (not in `NAMESPACE`), so the public R API is unchanged and the pipeline is opt-in via `xbioclim:::engine_set_pipeline(eng, TRUE)`.
- `-pthread` added to `PKG_LIBS` in `src/Makevars.in` so the C++ `std::thread` code links on Linux/Unix.

`BioclimEngine::compute()` performs the usual input validation, then either returns `compute_pipelined()` (when the internal `pipeline_` flag is true) or continues with the original serial tiled loop. The default behaviour is therefore identical to Phase A, and the new code can be disabled or removed by not setting the flag.

## 10. Future work

- Measure read/compute/write time per tile with `Rprintf` or internal counters to quantify overlap.
- Evaluate whether parallelising the four variable reads (each on its own `GdalReader`) is safe and beneficial for the target GDAL build.
- Adapt the slot layout to the pixel-major buffer redesign planned for Phase B.
- Consider a lock-free queue once the design is proven correct.
