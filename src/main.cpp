#include "xbioclim/bioclim.hpp"
#include "xbioclim/gdal_io.hpp"
#include "xbioclim/version.hpp"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Minimal argument parser — no external dependencies required.
// ---------------------------------------------------------------------------

static void print_usage(const char* progname) {
    std::cerr
        << "xbioclim v" << XBIOCLIM_VERSION_STRING << "\n\n"
        << "Usage: " << progname << " [options]\n\n"
        << "Required options (each expects exactly 12 file paths):\n"
        << "  --tas      <f1> ... <f12>   Monthly mean temperature files\n"
        << "  --tasmax   <f1> ... <f12>   Monthly max temperature files\n"
        << "  --tasmin   <f1> ... <f12>   Monthly min temperature files\n"
        << "  --pr       <f1> ... <f12>   Monthly precipitation files\n\n"
        << "Output options:\n"
        << "  --outdir   <dir>            Output directory (default: .)\n"
        << "  --prefix   <str>            Output filename prefix (default: bio)\n\n"
        << "Processing options:\n"
        << "  --tile-size <n>             Tile height in pixels (default: 512)\n"
        << "  --raw                       Skip CHELSA scale/offset decoding\n\n"
        << "CHELSA decoding (applied unless --raw):\n"
        << "  --tas-scale    <f>          (default: 0.1)\n"
        << "  --tas-offset   <f>          (default: -273.15)\n"
        << "  --tasmax-scale <f>          (default: 0.1)\n"
        << "  --tasmax-offset <f>         (default: -273.15)\n"
        << "  --tasmin-scale <f>          (default: 0.1)\n"
        << "  --tasmin-offset <f>         (default: -273.15)\n"
        << "  --pr-scale     <f>          (default: 0.1)\n"
        << "  --pr-offset    <f>          (default: 0.0)\n\n"
        << "Other:\n"
        << "  --help                      Show this message\n"
        << "  --version                   Show version\n";
}

// Collect up to 12 non-flag tokens after the current position.
static std::vector<std::string> collect_files(int argc, char* argv[], int& i) {
    std::vector<std::string> files;
    while (i + 1 < argc && argv[i + 1][0] != '-' && files.size() < 12) {
        files.push_back(argv[++i]);
    }
    return files;
}

int main(int argc, char* argv[]) {
    std::vector<std::string> tas_files, tasmax_files, tasmin_files, pr_files;
    std::string outdir  = ".";
    std::string prefix  = "bio";
    int tile_h          = 512;
    bool raw            = false;

    // Default CHELSA V2 scale/offset for temperature (K*10 → °C) and precip
    xbioclim::ScaleOffset tas_so    = {0.1f, -273.15f};
    xbioclim::ScaleOffset tasmax_so = {0.1f, -273.15f};
    xbioclim::ScaleOffset tasmin_so = {0.1f, -273.15f};
    xbioclim::ScaleOffset pr_so     = {0.1f, 0.0f};

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--help" || arg == "-h") { print_usage(argv[0]); return 0; }
        if (arg == "--version")             { std::cout << XBIOCLIM_VERSION_STRING << "\n"; return 0; }
        if (arg == "--raw")                 { raw = true; continue; }

        if (arg == "--tas")    { tas_files    = collect_files(argc, argv, i); continue; }
        if (arg == "--tasmax") { tasmax_files  = collect_files(argc, argv, i); continue; }
        if (arg == "--tasmin") { tasmin_files  = collect_files(argc, argv, i); continue; }
        if (arg == "--pr")     { pr_files     = collect_files(argc, argv, i); continue; }

        if (arg == "--outdir"  && i + 1 < argc) { outdir  = argv[++i]; continue; }
        if (arg == "--prefix"  && i + 1 < argc) { prefix  = argv[++i]; continue; }
        if (arg == "--tile-size" && i + 1 < argc) {
            try {
                tile_h = std::stoi(argv[++i]);
                if (tile_h <= 0) {
                    std::cerr << "Error: --tile-size must be a positive integer\n";
                    return 1;
                }
            } catch (const std::exception&) {
                std::cerr << "Error: invalid value for --tile-size\n";
                return 1;
            }
            continue;
        }

        // Helper lambda for parsing float options
        auto parse_float = [&](const char* name, float& dest) -> bool {
            try {
                dest = std::stof(argv[++i]);
                return true;
            } catch (const std::exception&) {
                std::cerr << "Error: invalid value for " << name << "\n";
                return false;
            }
        };

        if (arg == "--tas-scale"     && i + 1 < argc) { if (!parse_float("--tas-scale",     tas_so.scale))     return 1; continue; }
        if (arg == "--tas-offset"    && i + 1 < argc) { if (!parse_float("--tas-offset",    tas_so.offset))    return 1; continue; }
        if (arg == "--tasmax-scale"  && i + 1 < argc) { if (!parse_float("--tasmax-scale",  tasmax_so.scale))  return 1; continue; }
        if (arg == "--tasmax-offset" && i + 1 < argc) { if (!parse_float("--tasmax-offset", tasmax_so.offset)) return 1; continue; }
        if (arg == "--tasmin-scale"  && i + 1 < argc) { if (!parse_float("--tasmin-scale",  tasmin_so.scale))  return 1; continue; }
        if (arg == "--tasmin-offset" && i + 1 < argc) { if (!parse_float("--tasmin-offset", tasmin_so.offset)) return 1; continue; }
        if (arg == "--pr-scale"      && i + 1 < argc) { if (!parse_float("--pr-scale",      pr_so.scale))     return 1; continue; }
        if (arg == "--pr-offset"     && i + 1 < argc) { if (!parse_float("--pr-offset",     pr_so.offset))    return 1; continue; }

        std::cerr << "Unknown option: " << arg << "\n";
        print_usage(argv[0]);
        return 1;
    }

    // Validate inputs
    auto check = [](const std::vector<std::string>& v, const char* name) {
        if (v.size() != 12) {
            std::cerr << "Error: --" << name << " requires exactly 12 files, got "
                      << v.size() << "\n";
            return false;
        }
        return true;
    };
    if (!check(tas_files, "tas") || !check(tasmax_files, "tasmax") ||
        !check(tasmin_files, "tasmin") || !check(pr_files, "pr")) {
        return 1;
    }

    // Override scale/offset when --raw is specified
    if (raw) {
        tas_so = tasmax_so = tasmin_so = pr_so = {1.0f, 0.0f};
    }

    try {
        // Open input rasters
        xbioclim::GdalReader reader(tas_files, tasmax_files, tasmin_files,
                                     pr_files, tas_so, tasmax_so, tasmin_so, pr_so);

        const int x_size = reader.raster_x_size();
        const int y_size = reader.raster_y_size();

        // Create output writer
        xbioclim::GdalWriter writer(outdir, prefix, x_size, y_size,
                                     reader.projection_wkt(),
                                     reader.geotransform());

        // Process tile rows
        const int tile_w = x_size;  // full-width strips
        for (int y_off = 0; y_off < y_size; y_off += tile_h) {
            const int bh = std::min(tile_h, y_size - y_off);

            std::cerr << "Processing rows " << y_off << ".."
                      << (y_off + bh - 1) << " / " << y_size << "\n";

            xbioclim::ClimateBlock block = reader.read_block(0, y_off, tile_w, bh);
            xbioclim::BioBlock bio = xbioclim::compute_bioclim(block);
            writer.write_block(0, y_off, tile_w, bh, bio);
        }

        writer.finalize();
        std::cerr << "Done. Wrote 19 BIO rasters to " << outdir << "/\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
