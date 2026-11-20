# Scientific status and remaining gates

This document is a dated interpretation guide for the 2026-11-20 source import. Execution snapshots are not live monitoring. The source publication does not pause, resume or change the existing BUSCO or Panaroo jobs.

| Gate | Available evidence | Remaining requirement |
|---|---|---|
| G0: cohort and identity | 77 assembly labels, 72 canonical sequences, explicit 66-tip core cohort | Establish biological culture equivalence; reconcile aliases and freeze analysis-specific eligibility |
| G1: reads and assembly | Read/assembly inventories, available fastp and assembly metrics | Verify exact read-to-assembly linkage and historical subsampling parameters; omit read-only checks where no reads exist |
| G2: completeness and annotation | BUSCO workflow, capped queue, annotation and copy-state audits | Complete and validate remaining QC; harmonize fungal gene models and resolve model exceptions |
| G3: locus identity | Genomic BLAST, reference windows, orthogroup checks and copy evidence | Resolve missing, truncated, divergent or duplicated loci without choosing an arbitrary best hit |
| G4: phylogenetic interpretation | Bootstrap-supported marker analyses and exact-cohort comparisons | Core branch support/concordance, taxon-sampling sensitivity and unresolved type-qualified conflicts |

The granular original gate records are retained in `manifests/initial_gate_snapshot.json`; their recorded preparation status is historical and does not override later evidence. The saved execution status is separately dated. Metadata for phenotype assay units/endpoints, codes -1/-2 and R22/R23/R24 remains unknown where unsupported, as requested; no values have been invented.

## Taxonomic transfer

Use the largest sampled three-gene clade containing the verified type anchors of one species and no other species' type anchor. Preserve reported labels, original strain IDs and evidence separately. The analysis threshold for a supported enclosing marker clade is SH-aLRT >=80 and UFBoot >=95; topological membership below that threshold remains provisional. A singleton type does not automatically classify its neighbors.

The terreus/montoensis assignments are intermixed in the current core tree. Neither corresponding type genome is in that mixed 40-tip region. Their type references occur in the three-gene analysis. Citrinoterreus, hortae and pseudoterreus retain group boundaries between the current core and three-gene trees but differ internally. Type genomes actually present in the core remain explicitly identified; absent type references must not be inserted as artificial core tips.

## cgMLST gates

The current 66-genome matrix provides 4,661 strict single-copy candidates, or 6,129 at >=95% occupancy with no duplicated models. These counts do not define a validated scheme.

1. Freeze a species-balanced panel, provenance and independent holdout set; add missing verified type genomes.
2. Validate fungal exon-joined CDS, translation and representative isoforms against genomic coordinates.
3. Check orthology, duplication histories, synteny and candidate stability across annotation versions.
4. Pilot nucleotide allele calling with explicit missing/truncated/ambiguous/multicopy states.
5. Define a stable section-wide locus catalog and species-specific extensions with permanent identifiers.
6. Validate reproducibility, callability, discrimination and interpretation on independent genomes.
7. Release immutable allele definitions and schema versions; establish clinical thresholds only with separate evidence.

Citrinoterreus, hortae and pseudoterreus are plausible first pilots, but their small current panels do not establish general species-wide schemes. Terreus/montoensis needs conflict resolution first. One- and two-genome panels cannot establish population-wide core status. Full proposed gates are in `manifests/cgmlst_gates.json`.

## Evo-devo extension

Keep lovastatin architecture as an accessory trait. Validate genomic absence versus failed detection, reconcile biological units, fit trait-appropriate evolutionary models with phylogenetic uncertainty, and evaluate phenotype associations after population structure and ecology. Avoid treating a missing phenotype code as a quantitative value. The detailed proposal is recorded in `manifests/evodevo_gates.json`.
