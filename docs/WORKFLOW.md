# Workflow, dependencies and outputs

Solid arrows represent implemented paths. Dashed arrows lead to proposals or additional validation. A successful computational run does not automatically close its scientific gate.

```mermaid
flowchart TD
  A[Public NCBI and ENA metadata / study assemblies] --> B[Accession, version and checksum inventory]
  R[Raw reads when available] --> C[fastp / assembly QC / SPAdes history]
  B --> D[Canonical assemblies and alias review]
  C --> D
  D --> E[BUSCO genome QC / assembly metrics]
  D --> F[Fungal annotation and gene model checks]
  F --> G[Protein orthogroups / single-copy candidate matrix]
  G --> H[MAFFT / core concatenation / IQ-TREE]
  P[Taxonomic studies and verified type references] --> M[Curated BenA / CaM / RPB2 library]
  D --> X[BLAST and genomic marker recovery]
  M --> X
  F --> O[Orthogroup, paralog and coordinate cross-checks]
  X --> O
  O --> T[Single genes / pairs / combined three-gene trees]
  T --> N[Type-qualified monophyly and provisional assignments]
  H --> V[Exact shared-genome tree comparisons]
  T --> V
  N --> V
  V --> Q[RF / NMDS / clade recovery / quartets / tanglegrams]
  D --> L[Lovastatin BLAST / protein rescue / genomic architecture]
  L --> Z[Linear clusters and circular genome views]
  E --> K[QC and interpretation gates]
  O --> K
  Q --> K
  Z --> K
  K -.-> S[Broader type-genome sampling / core gCF and sCF / topology tests]
  G -.-> CG[Validated section-wide cgMLST backbone and species extensions]
  L -.-> EV[Structure-aware evo-devo and ancestral-state analyses]
```

## Source and version graph

```mermaid
flowchart LR
  W[Observed workspace source] --> HASH[SHA-256 identity]
  D[Original project commands and notebooks] --> HASH
  HASH --> OBJ[Immutable source objects]
  HASH --> IDX[Version and category manifest]
  W --> CUR[Latest observed script view]
  IDX --> DRY[Dry-run restoration plan]
  OBJ --> DRY
  DRY --> EXT[Separate external analysis workspace]
  INPUT[Inventory-listed data kept outside Git] --> EXT
  EXT --> RUN[Reviewed and resource-capped execution]
  RUN --> OUT[New versioned outputs outside Git]
```

## What each stage returns

| Stage | Required inputs | Implemented output | Interpretation gate |
|---|---|---|---|
| Provenance | Public records, study labels, file versions | Accession and source manifests | Sequence identity is not culture identity |
| QC | Assembly FASTA; reads when present | Assembly statistics and BUSCO evidence | Read-dependent checks are not applicable without reads |
| Annotation | Genome, repeats, fungal gene models | CDS/protein/transcript models | Multiple isoforms are not automatically separate genes |
| Marker extraction | Curated references and genome sequence | Candidate loci, coordinates and copy evidence | BLAST similarity alone does not establish orthology |
| Phylogeny | Validated aligned loci and partitions | Individual, combined and core trees | Missing loci, supports and reference sampling remain explicit |
| Taxonomy | Verified type reference identities | Type-exclusive sampled clades and name-transfer audit | Do not propagate names across mixed core clades |
| Congruence | Identical leaf identifiers in both trees | RF, NMDS, boundary recovery, quartets and tanglegrams | Core is a comparison baseline, not automatic ground truth |
| BGC | Reference genes and genomic coordinates | Cluster architecture, flanks and circular views | Failed detection differs from biological absence |
| cgMLST | Validated single-copy genomic CDS | Candidate inventories only at present | Alleles, callability and holdout validation are still required |

The individual source files and exact historical variants are linked from the script catalog and version manifest. These graphs describe the implemented workflow and its open extensions; they do not imply a single fully automated orchestrator exists.
