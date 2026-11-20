#!/usr/bin/env bash
set -euo pipefail
ulimit -v 1572864
cd '/mnt/c/Users/q287js/Documents/ChatGPT/Aspergillus genomics'
/usr/bin/python3 finish_variant_rescue.py
/usr/bin/python3 extract_terreus_genomic_markers.py variant_search_rescued_manifest.json variant_genomic_markers
/usr/bin/python3 compare_assembly_marker_versions.py
