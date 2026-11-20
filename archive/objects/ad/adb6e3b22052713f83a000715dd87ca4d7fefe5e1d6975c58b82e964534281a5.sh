#!/bin/bash
set -euo pipefail
ulimit -v 1572864
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd '/mnt/c/Users/q287js/Documents/ChatGPT/Aspergillus genomics'
exec /usr/bin/python3 infer_core66_combinations.py
