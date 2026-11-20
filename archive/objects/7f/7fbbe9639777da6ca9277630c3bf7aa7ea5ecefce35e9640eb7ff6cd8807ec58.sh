#!/usr/bin/env bash
set -euo pipefail
ulimit -v 1572864
BASE='/mnt/c/Users/q287js/Documents/ChatGPT/Aspergillus genomics'
cd "$BASE"
/usr/bin/python3 audit_assembly_variants.py
/usr/bin/python3 extract_terreus_genomic_markers.py variant_search_manifest.json variant_genomic_markers
/usr/bin/python3 prepare_historical_merge.py
cd taxonomy_reference_build/2026-11-10/analysis
for gene in BenA CaM; do
  /home/benjamin/miniforge3/envs/klebsiella_analysis/bin/mafft --addfragments "${gene}_historical_2011.fasta" --keeplength --mapout --thread 1 "${gene}_genomic_combined_aligned.fasta" > "${gene}_with_historical_aligned.fasta" 2> "logs/${gene}_historical_mafft.log"
done
cd "$BASE"
bash infer_terreus_trees_capped.sh
