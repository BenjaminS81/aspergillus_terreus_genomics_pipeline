#!/usr/bin/env bash
set -euo pipefail
ulimit -v 1572864
BASE='/mnt/c/Users/q287js/Documents/ChatGPT/Aspergillus genomics'
cd "$BASE"
/usr/bin/python3 finish_terreus_alignments.py
cd taxonomy_reference_build/2026-11-10/analysis
mkdir -p trees
IQ=/home/benjamin/miniforge3/envs/klebsiella_recombination/bin/iqtree2
"$IQ" -s tree_inputs/combined_all_strains.fasta -p tree_inputs/partitions.nex -m MFP+MERGE -B 1000 -alrt 1000 -bnni -seed 20261110 -nt 1 -o REF_436_Aspergillus_flavipes_NRRL302 --prefix trees/combined_all_strains > logs/combined_iqtree.stdout 2>&1
for gene in BenA CaM RPB2; do
  "$IQ" -s "tree_inputs/${gene}_all_strains.fasta" -m MFP -B 1000 -alrt 1000 -bnni -seed 20261110 -nt 1 -mem 1G -o REF_436_Aspergillus_flavipes_NRRL302 --prefix "trees/${gene}" > "logs/${gene}_iqtree.stdout" 2>&1
done
