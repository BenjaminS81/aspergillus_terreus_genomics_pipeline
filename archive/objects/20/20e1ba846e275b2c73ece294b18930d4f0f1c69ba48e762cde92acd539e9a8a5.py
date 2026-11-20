from pathlib import Path
import json
P=Path(__file__).resolve().parent/'taxonomy_reference_build/2026-11-10/analysis'
refs=json.loads((P/'validated_reference_metadata.json').read_text())
for l in ['BenA','CaM','RPB2']:
    (P/(l+'_genomic_combined_unaligned.fasta')).write_text(''.join('>'+d['reference_id']+'\n'+d[l+'_analysis_sequence']+'\n' for d in refs if d[l+'_analysis_sequence'])+(P/'genomic_markers'/(l+'_genomic.fasta')).read_text())
