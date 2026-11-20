#!/usr/bin/env python3
"""Build deduplicated genome/transcript inventories for complete Terreus QC."""

from __future__ import annotations

import csv
import sys
from pathlib import Path


def normalize_sample(name: str) -> str:
    return "ATCC20541" if name.upper() == "ATCC_20541" else name


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main(project_text: str, output_text: str) -> None:
    project = Path(project_text)
    output = Path(output_text)
    output.mkdir(parents=True, exist_ok=True)

    fastp_root = project / "results" / "0_0_fastp"
    read_pairs: dict[tuple[str, str], tuple[Path, Path]] = {}
    for read1 in fastp_root.glob("*/*_1.fastq.gz"):
        species = read1.parent.name
        sample = read1.name.removesuffix("_1.fastq.gz")
        read2 = read1.with_name(f"{sample}_2.fastq.gz")
        if read2.is_file() and read1.stat().st_size and read2.stat().st_size:
            read_pairs[(species, sample)] = (read1, read2)

    genomes: dict[tuple[str, str], dict] = {}
    spades_root = project / "results" / "0_1_spades"
    for fasta in sorted(spades_root.glob("*/*/scaffolds.fasta")):
        species, sample = fasta.parent.parent.name, fasta.parent.name
        key = (species, sample)
        genomes[key] = {
            "label": f"{species}__{sample}", "species_group": species, "sample": sample,
            "data_class": "SPAdes read-backed", "has_raw_reads": "yes" if key in read_pairs else "no",
            "source_fasta": str(fasta),
        }

    assembly_root = project / "data" / "assemblies"
    for annotation_class in ("annotated", "unannotated"):
        for sample_dir in sorted((assembly_root / annotation_class).iterdir()):
            if not sample_dir.is_dir():
                continue
            candidates = [p for p in sample_dir.glob("*_genomic.fna") if p.name != "cds_from_genomic.fna"]
            if len(candidates) != 1:
                raise RuntimeError(f"Expected one genomic FNA in {sample_dir}, found {len(candidates)}")
            sample = normalize_sample(sample_dir.name)
            key = ("terreus", sample)
            genomes[key] = {
                "label": f"terreus__{sample}", "species_group": "terreus", "sample": sample,
                "data_class": f"Reference {annotation_class}", "has_raw_reads": "no",
                "source_fasta": str(candidates[0]),
            }

    fun_root = project / "results" / "0_2_repeat_mask_funannotate"
    terreus_fun = fun_root / "terreus"
    for fasta in sorted(terreus_fun.glob("*/contigs_renamed.fasta")):
        sample = normalize_sample(fasta.parent.name)
        key = ("terreus", sample)
        if key in genomes:
            continue
        genomes[key] = {
            "label": f"terreus__{sample}", "species_group": "terreus", "sample": sample,
            "data_class": "FASTA-only local", "has_raw_reads": "no", "source_fasta": str(fasta),
        }

    genome_rows = sorted(genomes.values(), key=lambda r: (r["data_class"], r["species_group"], r["sample"]))
    if len(genome_rows) != 66:
        raise RuntimeError(f"Expected 66 distinct genome assemblies, found {len(genome_rows)}")

    transcript_candidates: dict[tuple[str, str], list[Path]] = {}
    for path in fun_root.glob("*/*/*/annotate_results/*transcripts.fa"):
        section = path.parents[3].name
        raw_sample = path.parents[2].name
        sample = normalize_sample(raw_sample)
        species = "hortae" if section == "hortae" else "terreus"
        transcript_candidates.setdefault((species, sample), []).append(path)

    transcripts = []
    for key, candidates in sorted(transcript_candidates.items()):
        if key not in genomes:
            continue
        species, sample = key
        preferred_names = {f"{sample}_transcripts.fa".lower(), f"{sample.replace('.', '.')}_transcripts.fa".lower()}
        preferred = [p for p in candidates if p.name.lower() in preferred_names]
        if not preferred:
            preferred = [p for p in candidates if p.name.lower().endswith(".transcripts.fa")]
        chosen = sorted(preferred or candidates, key=lambda p: (len(p.name), p.name))[0]
        genome = genomes[key]
        transcripts.append({
            "label": genome["label"], "species_group": species, "sample": sample,
            "data_class": genome["data_class"], "source_transcripts": str(chosen),
            "matched_genome_fasta": genome["source_fasta"],
        })

    if len(transcripts) != 59:
        raise RuntimeError(f"Expected 59 distinct final transcript sets, found {len(transcripts)}")

    write_tsv(output / "complete_genome_manifest.tsv", genome_rows,
              ["label", "species_group", "sample", "data_class", "has_raw_reads", "source_fasta"])
    write_tsv(output / "complete_transcript_manifest.tsv", transcripts,
              ["label", "species_group", "sample", "data_class", "source_transcripts", "matched_genome_fasta"])

    genome_keys = set(genomes)
    reads_without_assemblies = sorted(set(read_pairs) - genome_keys)
    with (output / "inventory_audit.txt").open("w", encoding="utf-8") as handle:
        handle.write(f"Distinct genome assemblies included: {len(genome_rows)}\n")
        for cls in sorted({r['data_class'] for r in genome_rows}):
            handle.write(f"  {cls}: {sum(r['data_class'] == cls for r in genome_rows)}\n")
        handle.write(f"Distinct final transcript sets included: {len(transcripts)}\n")
        handle.write(f"Trimmed read pairs found: {len(read_pairs)}\n")
        handle.write(f"Assemblies with matching trimmed reads: {sum(r['has_raw_reads'] == 'yes' for r in genome_rows)}\n")
        handle.write("Read pairs omitted because no assembly was found:\n")
        for species, sample in reads_without_assemblies:
            handle.write(f"  {species} / {sample}\n")
        handle.write("Known duplicate aliases normalized:\n  ATCC_20541 -> ATCC20541\n")
        handle.write("Known duplicate genome excluded:\n  data/reference/GCA_000149615.1_ASM14961v1_genomic.fna (duplicate NIH2624)\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: build_complete_inventory.py PROJECT_ROOT OUTPUT_DIR")
    main(sys.argv[1], sys.argv[2])
