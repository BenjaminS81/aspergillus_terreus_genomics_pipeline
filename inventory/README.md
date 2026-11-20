# Input inventory, without data payloads

The inventory is an auditable record of the inputs available, fetched or considered at the 2026-11-20 audit. It does not imply that every fetched item was included in every analysis or that every local file passed QC.

| Catalog | Meaning |
|---|---|
| [Genome cohort](genome_cohort.tsv) | 77 labels, provenance, accessions, canonical sequence hashes and tree inclusion |
| [Assembly assets](assembly_assets.json) / [variants](assembly_variants.json) | Alternative local assembly assets and source paths |
| [Public assembly cross-check](public_assembly_crosscheck.json) | Availability and relevance of queried public assembly records |
| [Verified public additions](public_additions.json) | New downloads, checksums and staged status; not silently inserted into existing trees |
| [Raw-read files](local_raw_read_files.tsv) | Local file presence and sizes; small/failed files are not thereby certified usable |
| [ENA runs queried](ena_runs_queried.tsv) | Metadata and remote retrieval locations; local presence is a separate field |
| [Reference strains](reference_strains.tsv) | Library records, type flags, locus accessions and explicit tree membership |
| [Fetched marker accessions](fetched_marker_accessions.tsv) | GenBank accession versions and organisms; sequence payloads omitted |
| [Historical reference strains](historical_reference_strains.json) | Additional historical study records, separately scoped from the current reference tree |
| [Lovastatin proteins](lovastatin_protein_references.json) | Nine fetched UniProt reference proteins; amino-acid sequences excluded |
| [Local public assembly files](local_public_assembly_files.tsv) | Observed public assembly file versions, sizes and locations |
| [Publications](publications.json) / [downloads](publication_downloads.json) | Bibliographic sources and retrieval provenance |
| [Taxa](taxa.tsv) / [taxon evidence](taxon_evidence.json) | Source and analysis names, queried-only versus fetched/used roles |
| [Eligibility](analysis_eligibility.json) | Analysis-specific eligibility snapshot |

`PROJECT/` denotes the external `terreus_genotyping` project root. `WORKSPACE/` denotes the external analysis workspace. These are inventory location tokens, not paths to files supplied by this repository. Raw sample metadata measurements are excluded; study labels remain as identifiers needed to audit inclusion.

“Sequence fetched” is supported by a recorded sequence retrieval, not merely a database search. A source name and a provisional type-derived assignment can differ. Taxon names are not silently collapsed into accepted synonyms. Catalog-only names appear because they were assessed as nomenclatural context, not because a genome of that taxon was downloaded.

The primary genotype panel, downloaded references, staged new genomes and metadata-only search results have separate scopes. Counts are not additive across those catalogs. Identical sequence hashes do not prove that labels refer to the same biological culture. WT/sector aliases supported by read/base counts remain provisional until sequence identity is established.
