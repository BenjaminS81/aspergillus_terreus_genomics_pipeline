#!/bin/bash
set -euo pipefail
ulimit -v 1572864
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd '/mnt/c/Users/q287js/Documents/ChatGPT/Aspergillus genomics'
/usr/bin/python3 rescue_lov_proteins.py
/usr/bin/python3 summarize_lov_genomic_screen.py
