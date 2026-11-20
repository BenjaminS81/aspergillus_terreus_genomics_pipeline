#!/usr/bin/env python3
"""Calculate basic nucleotide transcript-set metrics."""

from __future__ import annotations

import csv
import sys
from pathlib import Path


def fasta_stats(path: Path):
    lengths, current, total, gc, n_count = [], 0, 0, 0, 0
    with path.open("rt", encoding="ascii", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
            else:
                seq = line.upper()
                current += len(seq)
                total += len(seq)
                gc += seq.count("G") + seq.count("C")
                n_count += seq.count("N")
    if current:
        lengths.append(current)
    return lengths, total, gc, n_count


def nx(lengths, total, fraction):
    running = 0
    for index, length in enumerate(sorted(lengths, reverse=True), 1):
        running += length
        if running >= total * fraction:
            return length, index
    return 0, 0


def main(manifest_text, output_text):
    with open(manifest_text, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    fields = ["label", "species_group", "sample", "data_class", "source_transcripts",
              "transcripts", "total_transcript_bp", "min_transcript_bp", "max_transcript_bp",
              "mean_transcript_bp", "transcript_N50_bp", "transcript_L50", "GC_percent", "N_percent"]
    with open(output_text, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            lengths, total, gc, n_count = fasta_stats(Path(row["source_transcripts"]))
            n50, l50 = nx(lengths, total, 0.5)
            writer.writerow({
                "label": row["label"], "species_group": row["species_group"], "sample": row["sample"],
                "data_class": row["data_class"], "source_transcripts": row["source_transcripts"],
                "transcripts": len(lengths), "total_transcript_bp": total,
                "min_transcript_bp": min(lengths, default=0), "max_transcript_bp": max(lengths, default=0),
                "mean_transcript_bp": f"{total / len(lengths) if lengths else 0:.2f}",
                "transcript_N50_bp": n50, "transcript_L50": l50,
                "GC_percent": f"{100 * gc / total if total else 0:.4f}",
                "N_percent": f"{100 * n_count / total if total else 0:.6f}",
            })


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: transcript_metrics.py MANIFEST.tsv OUTPUT.tsv")
    main(sys.argv[1], sys.argv[2])
