#!/usr/bin/env bash
set -euo pipefail
ulimit -v 1572864
BASE='/mnt/c/Users/q287js/Documents/ChatGPT/Aspergillus genomics'
cd "$BASE"
/usr/bin/python3 rescue_divergent_terreus_markers.py
/usr/bin/python3 extract_terreus_genomic_markers.py genome_search_rescued_manifest.json
/usr/bin/python3 prepare_genomic_combined_alignments.py
cd taxonomy_reference_build/2026-11-10/analysis
for gene in BenA CaM RPB2; do
  /home/benjamin/miniforge3/envs/klebsiella_analysis/bin/mafft --localpair --maxiterate 1000 --adjustdirectionaccurately --thread 1 "${gene}_genomic_combined_unaligned.fasta" > "${gene}_genomic_combined_aligned.fasta" 2> "logs/${gene}_genomic_mafft.log"
done
