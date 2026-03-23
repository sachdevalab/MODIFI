#!/usr/bin/env python3
"""
Aggregate /usr/bin/time -v stats for ipdSummary vs MODIFI across reference sizes.

Paths follow benchmark/compare/batch_ipd.sh and batch_modifi.sh. Missing or empty
.time files (e.g. ipdSummary still running) are skipped for that tool/reference.

Output CSVs use CPU and wall time in hours, peak RSS in GB.
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

# batch_modifi.sh / batch_ipd.sh
REF_DIR = Path("/home/shuaiw/borg/contigs")
IPD_OUT = Path("/home/shuaiw/borg/paper/ipdsummary/soil_1/ipd.out2")
MODIFI_ROOT = Path("/home/shuaiw/borg/paper/ipdsummary/soil_1/modifi.out2")
OUT_DIR = Path("/home/shuaiw/MODIFI/tmp/figures/base_benchmark")
REFS = ["test_100.fa", "test_200.fa", "test_300.fa", "test_500.fa"]


def fasta_total_bases(fasta_path: Path) -> int | None:
    if not fasta_path.is_file():
        return None
    total = 0
    with fasta_path.open() as f:
        for line in f:
            if line.startswith(">"):
                continue
            total += len(line.strip())
    return total


def parse_gnu_time_verbose(path: Path) -> dict[str, float] | None:
    if not path.is_file() or path.stat().st_size == 0:
        return None
    text = path.read_text(errors="replace")

    um = re.search(r"User time \(seconds\):\s*([\d.]+)", text)
    sm = re.search(r"System time \(seconds\):\s*([\d.]+)", text)
    mm = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", text)
    em = re.search(
        r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)", text
    )
    if not (um and sm and em):
        return None

    user_s = float(um.group(1))
    sys_s = float(sm.group(1))
    wall_spec = em.group(1).strip()
    parts = wall_spec.split(":")
    try:
        if len(parts) == 3:
            wall_s = int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
        elif len(parts) == 2:
            wall_s = int(parts[0]) * 60 + float(parts[1])
        else:
            return None
    except ValueError:
        return None

    max_rss_kb = int(mm.group(1)) if mm else -1
    if max_rss_kb < 0:
        return None

    return {
        "cpu_hr": (user_s + sys_s) / 3600.0,
        "wall_hr": wall_s / 3600.0,
        "max_rss_gb": max_rss_kb / (1024.0 * 1024.0),
    }


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    wide_path = OUT_DIR / "ipd_vs_modifi_runtime_wide.csv"
    long_path = OUT_DIR / "ipd_vs_modifi_runtime_long.csv"

    wide_rows: list[dict[str, str]] = []
    long_rows: list[dict[str, str]] = []

    for ref_name in REFS:
        base = Path(ref_name).stem
        ref_path = REF_DIR / ref_name
        ref_bases = fasta_total_bases(ref_path)
        ref_bases_s = str(ref_bases) if ref_bases is not None else ""

        row: dict[str, str] = {
            "reference_fa": ref_name,
            "ref_bases": ref_bases_s,
            "ipd_cpu_hr": "",
            "ipd_wall_hr": "",
            "ipd_max_rss_gb": "",
            "modifi_cpu_hr": "",
            "modifi_wall_hr": "",
            "modifi_max_rss_gb": "",
        }

        ipd_time = IPD_OUT / f"{base}.ipdSummary.time"
        ipd = parse_gnu_time_verbose(ipd_time)
        if ipd:
            row["ipd_cpu_hr"] = f'{ipd["cpu_hr"]:.6f}'
            row["ipd_wall_hr"] = f'{ipd["wall_hr"]:.6f}'
            row["ipd_max_rss_gb"] = f'{ipd["max_rss_gb"]:.6f}'
            long_rows.append(
                {
                    "reference_fa": ref_name,
                    "ref_bases": ref_bases_s,
                    "software": "ipdSummary",
                    "cpu_hr": row["ipd_cpu_hr"],
                    "wall_hr": row["ipd_wall_hr"],
                    "max_rss_gb": row["ipd_max_rss_gb"],
                    "time_file": str(ipd_time),
                }
            )

        mod_time = MODIFI_ROOT / base / "modifi.host.time"
        mod = parse_gnu_time_verbose(mod_time)
        if mod:
            row["modifi_cpu_hr"] = f'{mod["cpu_hr"]:.6f}'
            row["modifi_wall_hr"] = f'{mod["wall_hr"]:.6f}'
            row["modifi_max_rss_gb"] = f'{mod["max_rss_gb"]:.6f}'
            long_rows.append(
                {
                    "reference_fa": ref_name,
                    "ref_bases": ref_bases_s,
                    "software": "MODIFI",
                    "cpu_hr": row["modifi_cpu_hr"],
                    "wall_hr": row["modifi_wall_hr"],
                    "max_rss_gb": row["modifi_max_rss_gb"],
                    "time_file": str(mod_time),
                }
            )

        wide_rows.append(row)

    wide_fields = [
        "reference_fa",
        "ref_bases",
        "ipd_cpu_hr",
        "ipd_wall_hr",
        "ipd_max_rss_gb",
        "modifi_cpu_hr",
        "modifi_wall_hr",
        "modifi_max_rss_gb",
    ]
    with wide_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=wide_fields)
        w.writeheader()
        w.writerows(wide_rows)

    long_fields = [
        "reference_fa",
        "ref_bases",
        "software",
        "cpu_hr",
        "wall_hr",
        "max_rss_gb",
        "time_file",
    ]
    with long_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=long_fields)
        w.writeheader()
        w.writerows(long_rows)

    print(f"Wrote {wide_path}", file=sys.stderr)
    print(f"Wrote {long_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
