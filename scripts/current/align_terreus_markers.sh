#!/usr/bin/env bash
set -euo pipefail
cd '/mnt/c/Users/q287js/Documents/ChatGPT/Aspergillus genomics/taxonomy_reference_build/2026-11-10/analysis'
MAFFT=/home/benjamin/miniforge3/envs/klebsiella_analysis/bin/mafft
for gene in BenA CaM RPB2; do
  "$MAFFT" --localpair --maxiterate 1000 --adjustdirectionaccurately --thread 4 "${gene}_combined_unaligned.fasta" > "${gene}_combined_aligned.fasta" 2> "logs/${gene}_mafft.log"
done
