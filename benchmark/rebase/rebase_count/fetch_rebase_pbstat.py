#!/usr/bin/env python3

import csv
import urllib.request
from collections import Counter
from pathlib import Path

import matplotlib
from bs4 import BeautifulSoup


matplotlib.use("Agg")
import matplotlib.pyplot as plt


URL = "https://rebase.neb.com/rebase/pbstatlist.html"
OUTPUT_CSV = Path("rebase_pbstatlist.csv")
COUNTS_CSV = Path("unique_recseqs_pacbio_counts.csv")
SUMMARY_TXT = Path("unique_recseqs_pacbio_summary.txt")
PLOT_PNG = Path("unique_recseqs_pacbio_distribution.png")


def normalize_text(text: str) -> str:
    text = text.replace("\xa0", " ")
    return " ".join(text.split()).strip()


def expand_row(tr) -> list[str]:
    cells: list[str] = []
    for cell in tr.find_all(["td", "th"], recursive=False):
        text = normalize_text(cell.get_text(" ", strip=True))
        colspan = int(cell.get("colspan", 1))
        cells.extend([text] * colspan)
    return cells


def combine_headers(group_headers: list[str], sub_headers: list[str]) -> list[str]:
    combined: list[str] = []
    for group, sub in zip(group_headers, sub_headers):
        if group and sub and group != sub:
            combined.append(f"{group} {sub}")
        else:
            combined.append(group or sub or "")

    seen: Counter[str] = Counter()
    unique_headers: list[str] = []
    for header in combined:
        header = header or "separator"
        seen[header] += 1
        if header == "separator":
            unique_headers.append(f"separator_{seen[header]}")
        elif seen[header] == 1:
            unique_headers.append(header)
        else:
            unique_headers.append(f"{header}_{seen[header]}")
    return unique_headers


def fetch_html(url: str) -> bytes:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "Mozilla/5.0"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read()


def numeric_values(rows: list[list[str]], headers: list[str], column_name: str) -> list[int]:
    try:
        column_index = headers.index(column_name)
    except ValueError as exc:
        raise RuntimeError(f"Missing expected column: {column_name}") from exc

    values: list[int] = []
    for row in rows:
        value = row[column_index]
        if value in ("", "-"):
            continue
        values.append(int(value))
    return values


def main() -> None:
    html = fetch_html(URL)
    soup = BeautifulSoup(html, "html.parser")
    tables = soup.find_all("table")
    if len(tables) < 2:
        raise RuntimeError("Could not find the PacBio statistics table.")

    table = tables[1]
    rows = table.find_all("tr", recursive=False)
    if len(rows) < 3:
        raise RuntimeError("PacBio statistics table is missing expected rows.")

    header_group = expand_row(rows[0])
    header_sub = expand_row(rows[1])
    headers = combine_headers(header_group, header_sub)

    total_genomes_from_footer = None
    raw_rows: list[list[str]] = []
    for tr in rows[2:]:
        row = expand_row(tr)
        if len(row) != len(headers):
            continue

        first_cell = row[0]
        if first_cell == "Organisms" or row[13] == "Unique recseqs":
            continue
        if first_cell.startswith("Totals:"):
            total_genomes_from_footer = int(first_cell.split(":")[1].strip())
            continue
        if not first_cell:
            continue

        raw_rows.append(row)

    keep_indices = [idx for idx, name in enumerate(headers) if not name.startswith("separator_")]
    kept_headers = [headers[idx] for idx in keep_indices]
    cleaned_rows = [[row[idx] for idx in keep_indices] for row in raw_rows]

    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(kept_headers)
        writer.writerows(cleaned_rows)

    genome_total = len(cleaned_rows)
    if total_genomes_from_footer is not None and total_genomes_from_footer != genome_total:
        raise RuntimeError(
            f"Genome count mismatch: footer says {total_genomes_from_footer}, parsed {genome_total}."
        )

    pacbio_counts = numeric_values(cleaned_rows, kept_headers, "Unique recseqs PacBio")
    rebase_counts = numeric_values(cleaned_rows, kept_headers, "Unique recseqs REBASE")

    distribution = Counter(pacbio_counts)
    with COUNTS_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["motif_count", "genome_count"])
        for motif_count in sorted(distribution):
            writer.writerow([motif_count, distribution[motif_count]])

    pacbio_minimum = min(pacbio_counts)
    pacbio_maximum = max(pacbio_counts)
    pacbio_average = sum(pacbio_counts) / len(pacbio_counts)
    pacbio_total = sum(pacbio_counts)

    rebase_minimum = min(rebase_counts)
    rebase_maximum = max(rebase_counts)
    rebase_average = sum(rebase_counts) / len(rebase_counts)
    rebase_total = sum(rebase_counts)

    with SUMMARY_TXT.open("w", encoding="utf-8") as handle:
        handle.write(f"Source URL: {URL}\n")
        handle.write(f"Genome total in table: {genome_total}\n")
        handle.write("Rule: genomes with '-' are excluded from all counts below.\n")
        handle.write(f"Unique recseqs PacBio genomes counted: {len(pacbio_counts)}\n")
        handle.write(f"Unique recseqs PacBio total: {pacbio_total}\n")
        handle.write(f"Unique recseqs PacBio min: {pacbio_minimum}\n")
        handle.write(f"Unique recseqs PacBio max: {pacbio_maximum}\n")
        handle.write(f"Unique recseqs PacBio average: {pacbio_average:.6f}\n")
        handle.write(f"Unique recseqs REBASE genomes counted: {len(rebase_counts)}\n")
        handle.write(f"Unique recseqs REBASE total: {rebase_total}\n")
        handle.write(f"Unique recseqs REBASE min: {rebase_minimum}\n")
        handle.write(f"Unique recseqs REBASE max: {rebase_maximum}\n")
        handle.write(f"Unique recseqs REBASE average: {rebase_average:.6f}\n")

    x_values = sorted(distribution)
    y_values = [distribution[value] for value in x_values]

    plt.figure(figsize=(12, 6))
    plt.bar(x_values, y_values, width=0.8, color="#4c78a8")
    plt.title("Distribution of Unique recseqs PacBio Motif Counts")
    plt.xlabel("Motif count per genome")
    plt.ylabel("Number of genomes (excluding '-')")
    plt.xticks(x_values, rotation=90)
    plt.tight_layout()
    plt.savefig(PLOT_PNG, dpi=200)
    plt.close()

    print(f"Wrote {OUTPUT_CSV}")
    print(f"Wrote {COUNTS_CSV}")
    print(f"Wrote {SUMMARY_TXT}")
    print(f"Wrote {PLOT_PNG}")
    print(f"Genome total in table: {genome_total}")
    print(
        "PacBio counted / total / min / max / avg: "
        f"{len(pacbio_counts)} / {pacbio_total} / {pacbio_minimum} / {pacbio_maximum} / {pacbio_average:.6f}"
    )
    print(
        "REBASE counted / total / min / max / avg: "
        f"{len(rebase_counts)} / {rebase_total} / {rebase_minimum} / {rebase_maximum} / {rebase_average:.6f}"
    )


if __name__ == "__main__":
    main()
