#!/usr/bin/env python3
"""Merge complete genome, read-support, BUSCO, and transcript QC results."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path


def read_tsv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_busco(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        obj = json.load(handle)
    result = obj["results"]
    return {
        "busco_complete_percent": result.get("Complete percentage", ""),
        "busco_single_percent": result.get("Single copy percentage", ""),
        "busco_duplicated_percent": result.get("Multi copy percentage", ""),
        "busco_fragmented_percent": result.get("Fragmented percentage", ""),
        "busco_missing_percent": result.get("Missing percentage", ""),
        "busco_markers": result.get("n_markers", ""),
        "busco_internal_stop_percent": result.get("internal_stop_codon_percent", ""),
    }


def first_specific_json(folder: Path) -> Path | None:
    matches = sorted(folder.glob("short_summary.specific.*.json"))
    if matches:
        return matches[0]
    matches = sorted(folder.rglob("short_summary*.json")) if folder.exists() else []
    return matches[0] if matches else None


def main(root_text: str) -> None:
    root = Path(root_text)
    complete = root / "complete_data"
    assembly_rows = read_tsv(complete / "complete_assembly_metrics.tsv")
    old_qc = {r["label"]: r for r in read_tsv(root / "terreus_assembly_qc_summary.tsv")}
    transcript_metrics = {r["label"]: r for r in read_tsv(complete / "complete_transcript_metrics.tsv")}

    mapping_fields = ["reads_total", "reads_mapped", "mapped_percent", "reads_properly_paired",
                      "properly_paired_percent", "coverage_breadth_percent", "mean_depth"]
    genome_output = []
    availability = []
    for row in assembly_rows:
        label = row["label"]
        combined = dict(row)
        mapping = old_qc.get(label, {})
        combined.update({field: mapping.get(field, "") for field in mapping_fields})

        old_folder = root / "busco_results" / "terreus_genome_busco" / f"{label}.fasta"
        new_folder = complete / "busco_genomes" / label
        busco_json = first_specific_json(old_folder) or first_specific_json(new_folder)
        if busco_json:
            combined.update(parse_busco(busco_json))
        genome_output.append(combined)

        transcript_folder = complete / "busco_transcripts" / label
        availability.append({
            "label": label, "species_group": row["species_group"], "sample": row["sample"],
            "data_class": row["data_class"], "genome_fasta": "yes",
            "genome_busco": "yes" if busco_json else "no",
            "raw_reads": row["has_raw_reads"],
            "read_mapping_qc": "yes" if mapping else "no",
            "final_transcripts": "yes" if label in transcript_metrics else "no",
            "transcript_busco": "yes" if first_specific_json(transcript_folder) else "no",
        })

    transcript_output = []
    for label, row in sorted(transcript_metrics.items()):
        combined = dict(row)
        busco_json = first_specific_json(complete / "busco_transcripts" / label)
        if busco_json:
            combined.update(parse_busco(busco_json))
        transcript_output.append(combined)

    write_tsv(complete / "complete_genome_qc_summary.tsv", genome_output)
    write_tsv(complete / "complete_transcript_qc_summary.tsv", transcript_output)
    write_tsv(complete / "data_availability.tsv", availability)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: summarize_complete_qc.py QC_ROOT")
    main(sys.argv[1])
