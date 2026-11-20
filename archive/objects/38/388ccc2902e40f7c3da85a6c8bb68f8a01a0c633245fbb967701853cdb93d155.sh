#!/usr/bin/env bash
set -euo pipefail
ulimit -v 1572864
BASE='/mnt/c/Users/q287js/Documents/ChatGPT/Aspergillus genomics'
cd "$BASE"
/usr/bin/python3 prepare_type_sensitivity.py
cd taxonomy_reference_build/2026-11-10/analysis
/home/benjamin/miniforge3/envs/klebsiella_recombination/bin/iqtree2 -s tree_inputs/combined_type_focused_sensitivity.fasta -p trees/combined_all_strains.best_scheme.nex -B 1000 -alrt 1000 -bnni -seed 20261110 -nt 1 -o REF_436_Aspergillus_flavipes_NRRL302 --prefix trees/type_focused_sensitivity > logs/type_focused_sensitivity.stdout 2>&1
cd "$BASE"
bash finish_variant_rescue_capped.sh
