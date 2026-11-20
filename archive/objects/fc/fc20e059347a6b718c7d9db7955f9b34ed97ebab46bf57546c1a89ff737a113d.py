#!/usr/bin/env python3
"""Combine per-sample BUSCO and read-mapping results into one QC table."""

from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path


def read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_flagstat(path):
    values = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            passed, failed, label = fields[0], fields[1], fields[2]
            if not passed.isdigit() or not failed.isdigit():
                continue
            values[label] = int(passed) + int(failed)
    total = values.get("total (QC-passed reads + QC-failed reads)", 0)
    mapped = values.get("mapped", 0)
    proper = values.get("properly paired", 0)
    return {
        "reads_total": total,
        "reads_mapped": mapped,
        "mapped_percent": 100 * mapped / total if total else 0,
        "reads_properly_paired": proper,
        "properly_paired_percent": 100 * proper / total if total else 0,
    }


def parse_busco_json(path):
    with path.open(encoding="utf-8") as handle:
        obj = json.load(handle)
    results = obj.get("results", obj)
    keys = {
        "Complete percentage": "busco_complete_percent",
        "Single copy percentage": "busco_single_percent",
        "Multi copy percentage": "busco_duplicated_percent",
        "Fragmented percentage": "busco_fragmented_percent",
        "Missing percentage": "busco_missing_percent",
        "n_markers": "busco_markers",
        "internal_stop_codon_percent": "busco_internal_stop_percent",
    }
    return {out: results.get(source, "") for source, out in keys.items()}


def main(root_text):
    root = Path(root_text)
    metrics = {row["label"]: row for row in read_tsv(root / "metrics" / "assembly_metrics.tsv")}
    rows = []
    for label, row in sorted(metrics.items()):
        depth_rows = read_tsv(root / "read_mapping" / "per_sample" / label / "depth_summary.tsv")
        depth = depth_rows[0] if depth_rows else {}
        mapping = parse_flagstat(root / "read_mapping" / "per_sample" / label / "flagstat.tsv")
        combined = dict(row)
        combined.update(mapping)
        combined.update({
            "coverage_breadth_percent": depth.get("coverage_breadth_percent", ""),
            "mean_depth": depth.get("mean_depth", ""),
        })

        busco_sample_dir = root / "busco_results" / "terreus_genome_busco" / f"{label}.fasta"
        busco_candidates = list(busco_sample_dir.glob("short_summary.specific.*.json"))
        if not busco_candidates:
            busco_candidates = list(busco_sample_dir.rglob("short_summary*.json"))
        if busco_candidates:
            combined.update(parse_busco_json(busco_candidates[0]))
        rows.append(combined)

    all_fields = []
    for row in rows:
        for key in row:
            if key not in all_fields:
                all_fields.append(key)
    output = root / "terreus_assembly_qc_summary.tsv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=all_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: summarize_qc.py QC_ROOT")
    main(sys.argv[1])
