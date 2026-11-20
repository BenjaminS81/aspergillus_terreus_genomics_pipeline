#!/usr/bin/env python3
"""Calculate reference-free FASTA assembly metrics from a TSV manifest."""

from __future__ import annotations

import csv
import sys
from pathlib import Path


def fasta_lengths_and_bases(path: Path):
    lengths = []
    current = 0
    gc = 0
    n_count = 0
    total = 0
    with path.open("rt", encoding="ascii", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
                continue
            seq = line.upper()
            current += len(seq)
            total += len(seq)
            gc += seq.count("G") + seq.count("C")
            n_count += seq.count("N")
    if current:
        lengths.append(current)
    return lengths, total, gc, n_count


def nx(lengths, total, fraction):
    threshold = total * fraction
    running = 0
    for index, length in enumerate(sorted(lengths, reverse=True), start=1):
        running += length
        if running >= threshold:
            return length, index
    return 0, 0


def main(manifest_path: str, output_path: str):
    with open(manifest_path, newline="", encoding="utf-8") as source:
        manifest = list(csv.DictReader(source, delimiter="\t"))

    fields = [
        "label", "species_group", "sample", "data_class", "has_raw_reads", "source_fasta", "sequences",
        "total_length_bp", "min_length_bp", "max_length_bp", "mean_length_bp",
        "N50_bp", "L50", "N90_bp", "L90", "GC_percent", "N_percent",
        "sequences_ge_1kb", "sequences_ge_10kb", "sequences_ge_50kb",
    ]
    with open(output_path, "w", newline="", encoding="utf-8") as destination:
        writer = csv.DictWriter(destination, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in manifest:
            path = Path(row["source_fasta"])
            lengths, total, gc, n_count = fasta_lengths_and_bases(path)
            n50, l50 = nx(lengths, total, 0.5)
            n90, l90 = nx(lengths, total, 0.9)
            writer.writerow({
                "label": row["label"],
                "species_group": row["species_group"],
                "sample": row["sample"],
                "data_class": row.get("data_class", "SPAdes read-backed"),
                "has_raw_reads": row.get("has_raw_reads", "yes"),
                "source_fasta": row["source_fasta"],
                "sequences": len(lengths),
                "total_length_bp": total,
                "min_length_bp": min(lengths, default=0),
                "max_length_bp": max(lengths, default=0),
                "mean_length_bp": f"{(total / len(lengths)) if lengths else 0:.2f}",
                "N50_bp": n50,
                "L50": l50,
                "N90_bp": n90,
                "L90": l90,
                "GC_percent": f"{(100 * gc / total) if total else 0:.4f}",
                "N_percent": f"{(100 * n_count / total) if total else 0:.6f}",
                "sequences_ge_1kb": sum(x >= 1_000 for x in lengths),
                "sequences_ge_10kb": sum(x >= 10_000 for x in lengths),
                "sequences_ge_50kb": sum(x >= 50_000 for x in lengths),
            })


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: assembly_metrics.py MANIFEST.tsv OUTPUT.tsv")
    main(sys.argv[1], sys.argv[2])
