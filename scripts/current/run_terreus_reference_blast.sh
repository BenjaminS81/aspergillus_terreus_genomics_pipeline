#!/usr/bin/env bash
set -euo pipefail
cd '/mnt/c/Users/q287js/Documents/ChatGPT/Aspergillus genomics/taxonomy_reference_build/2026-11-10/analysis'
mkdir -p logs
for gene in BenA CaM RPB2; do
  /usr/bin/blastn -query "${gene}_type_queries.fasta" -subject "${gene}_local_ungapped.fasta" -task blastn -evalue 1e-10 -dust no -num_threads 4 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -out "${gene}_reference_to_local.blast.tsv"
done
