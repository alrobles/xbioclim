#!/usr/bin/env python3
"""Parse xbioclim benchmark output and emit JSON + Markdown reports.

Usage:
    python3 tools/parse_benchmark.py --input <dir> --output <dir> [--tag <tag>]

The script scans all *.txt and *.log files under --input, extracts timing
information from ``time`` command output and Catch2 benchmark sections,
collects hardware metadata, and writes:

  <output>/benchmark-report.json
  <output>/benchmark-report.md
"""

import argparse
import json
import os
import platform
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


# ---------------------------------------------------------------------------
# Hardware info helpers
# ---------------------------------------------------------------------------

def _run(cmd: list[str], default: str = "unknown") -> str:
    try:
        return subprocess.check_output(cmd, stderr=subprocess.DEVNULL, text=True).strip()
    except Exception:
        return default


def collect_hardware_info() -> dict:
    info: dict = {}

    info["os"] = platform.platform()
    info["arch"] = platform.machine()
    info["python"] = platform.python_version()

    # CPU model
    cpu_model = "unknown"
    try:
        with open("/proc/cpuinfo") as fh:
            for line in fh:
                if line.startswith("model name"):
                    cpu_model = line.split(":", 1)[1].strip()
                    break
    except OSError:
        cpu_model = _run(["sysctl", "-n", "machdep.cpu.brand_string"])
    info["cpu_model"] = cpu_model

    # Physical CPU count
    info["cpu_count_logical"] = os.cpu_count() or 0

    # RAM
    ram_kb = 0
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemTotal:"):
                    ram_kb = int(line.split()[1])
                    break
    except OSError:
        pass
    info["ram_mb"] = ram_kb // 1024

    # Kernel / uname
    info["uname"] = _run(["uname", "-r"])

    # Compiler versions (best-effort)
    gcc_path = _run(["which", "gcc"])
    info["gcc_version"] = _run(["gcc", "--version"]).splitlines()[0] if gcc_path != "unknown" else "n/a"
    clang_path = _run(["which", "clang"])
    info["clang_version"] = _run(["clang", "--version"]).splitlines()[0] if clang_path != "unknown" else "n/a"

    return info


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

# Matches ``time`` built-in or /usr/bin/time output lines like:
#   real  0m3.142s
#   user  0m2.980s
#   sys   0m0.162s
# or GNU time:
#   0:03.14 elapsed
_TIME_RE = re.compile(
    r"^(?P<kind>real|user|sys)\s+(?P<min>\d+)m(?P<sec>[\d.]+)s",
    re.MULTILINE,
)

# Catch2 benchmark line:
#   BIO throughput    10    ...    12345 ns    ...
_CATCH2_BENCH_RE = re.compile(
    r"^(?P<name>\S.*?)\s{2,}(?P<samples>\d+)\s+.*?(?P<mean>[\d.]+)\s+ns",
    re.MULTILINE,
)


def parse_time_output(text: str) -> dict | None:
    """Return {real_s, user_s, sys_s} from ``time`` command output, or None."""
    matches = {}
    for m in _TIME_RE.finditer(text):
        seconds = float(m.group("min")) * 60 + float(m.group("sec"))
        matches[m.group("kind")] = seconds
    if "real" in matches:
        return {k + "_s": v for k, v in matches.items()}
    return None


def parse_catch2_benchmarks(text: str) -> list[dict]:
    """Return list of {name, samples, mean_ns} from Catch2 benchmark output."""
    results = []
    for m in _CATCH2_BENCH_RE.finditer(text):
        results.append({
            "name": m.group("name").strip(),
            "samples": int(m.group("samples")),
            "mean_ns": float(m.group("mean")),
        })
    return results


def parse_file(path: Path) -> dict:
    text = path.read_text(errors="replace")
    result: dict = {"file": path.name}

    timing = parse_time_output(text)
    if timing:
        result["timing"] = timing

    catch2 = parse_catch2_benchmarks(text)
    if catch2:
        result["catch2_benchmarks"] = catch2

    # Keep first 40 lines of raw output for reference
    result["raw_head"] = "\n".join(text.splitlines()[:40])

    return result


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def build_report(input_dir: Path, tag: str) -> dict:
    hardware = collect_hardware_info()

    files_parsed = []
    for suffix in ("*.txt", "*.log"):
        for p in sorted(input_dir.glob(suffix)):
            files_parsed.append(parse_file(p))

    return {
        "tag": tag,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "hardware": hardware,
        "results": files_parsed,
    }


def render_markdown(report: dict) -> str:
    lines = [
        f"# xbioclim Benchmark Report — {report['tag']}",
        "",
        f"**Generated:** {report['generated_at']}",
        "",
        "## Hardware",
        "",
        "| Property | Value |",
        "| --- | --- |",
    ]
    hw = report["hardware"]
    for key, val in hw.items():
        lines.append(f"| {key} | {val} |")

    lines += ["", "## Results", ""]

    results = report.get("results", [])
    if not results:
        lines.append("_No benchmark output files were found._")
    else:
        for res in results:
            lines.append(f"### `{res['file']}`")
            lines.append("")
            timing = res.get("timing")
            if timing:
                lines += [
                    "| Metric | Seconds |",
                    "| --- | --- |",
                ]
                for k, v in timing.items():
                    lines.append(f"| {k} | {v:.3f} |")
                lines.append("")
            catch2 = res.get("catch2_benchmarks", [])
            if catch2:
                lines += [
                    "| Benchmark | Samples | Mean (ns) |",
                    "| --- | --- | --- |",
                ]
                for b in catch2:
                    lines.append(f"| {b['name']} | {b['samples']} | {b['mean_ns']:.1f} |")
                lines.append("")
            raw = res.get("raw_head", "")
            if raw:
                lines += [
                    "<details>",
                    "<summary>Raw output (first 40 lines)</summary>",
                    "",
                    "```",
                    raw,
                    "```",
                    "</details>",
                    "",
                ]

    lines += [
        "---",
        "",
        "_Report generated by `tools/parse_benchmark.py`_",
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Directory with benchmark output files")
    parser.add_argument("--output", required=True, help="Directory to write JSON/MD reports")
    parser.add_argument("--tag", default="dev", help="Release tag (e.g. v0.3.0)")
    args = parser.parse_args()

    input_dir = Path(args.input)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_dir.is_dir():
        print(f"[parse_benchmark] WARNING: input directory '{input_dir}' does not exist; "
              "writing empty report.", file=sys.stderr)
        input_dir = output_dir  # will find no files → empty results

    report = build_report(input_dir, args.tag)

    json_path = output_dir / "benchmark-report.json"
    md_path = output_dir / "benchmark-report.md"

    json_path.write_text(json.dumps(report, indent=2))
    md_path.write_text(render_markdown(report))

    print(f"[parse_benchmark] JSON report  → {json_path}")
    print(f"[parse_benchmark] MD  report   → {md_path}")
    print(f"[parse_benchmark] Files parsed : {len(report['results'])}")


if __name__ == "__main__":
    main()
