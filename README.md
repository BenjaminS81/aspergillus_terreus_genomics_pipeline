# Aspergillus section Terrei genomics

[![Repository checks](https://github.com/BenjaminS81/aspergillus_terreus_genomics_pipeline/actions/workflows/validate.yml/badge.svg)](https://github.com/BenjaminS81/aspergillus_terreus_genomics_pipeline/actions/workflows/validate.yml)
[![Source-only policy](https://github.com/BenjaminS81/aspergillus_terreus_genomics_pipeline/actions/workflows/source-policy.yml/badge.svg)](https://github.com/BenjaminS81/aspergillus_terreus_genomics_pipeline/actions/workflows/source-policy.yml)

Versioned methods, scripts and provenance for assembly QC, fungal marker orthology, type-anchored taxonomy, genome phylogeny, lovastatin-cluster architecture and clade-level congruence in *Aspergillus* section Terrei.

This repository preserves the available implemented source versions and inventories their external inputs. **Sequence data, alignments, trees, figures, raw analytical results and notebook outputs are not distributed here.** It is a research workflow archive with tested repository utilities, not a claim that every historical script is portable or every scientific gate is closed. Badge checks validate repository integrity and source-only policy, not biological conclusions.

## Start here

- [Workflow and analysis graph](docs/WORKFLOW.md): inputs, tools, outputs and decision gates.
- [Current script catalog](docs/SCRIPT_CATALOG.md): scripts grouped by purpose.
- [Version and restoration guide](docs/VERSIONS.md): immutable snapshots, provenance and external-workspace restoration.
- [Data inventory](inventory/README.md): genomes, reads, marker accessions, publications and download status.
- [Genera and species](docs/TAXA.md): complete names within the audited scope, with fetched/used/catalog-only evidence.
- [Scientific status and gates](docs/GATES.md): completed work, open issues and cgMLST/evo-devo proposals.
- [Methods and exact parameters](docs/METHODS.md): recorded settings and interpretation limits.
- [Reproducibility](docs/REPRODUCIBILITY.md), [security](SECURITY.md), [contribution rules](CONTRIBUTING.md), [rights and attribution](NOTICE.md).

## Repository structure

```text
archive/objects/       Distinct historical source versions, addressed by SHA-256
scripts/current/       Latest observed top-level workspace scripts
tools/                 Tested, data-free inspection and restoration utilities
tests/                 Integrity, safe restoration and phylogeny helper tests
config/                Validated configuration example; no sample measurements
inventory/             Metadata-only input, accession and taxon catalogs
manifests/             Source versions, hashes, analysis versions and dated gates
docs/                  Workflow graphs, methods, limitations and operating guide
.github/               Pinned CI workflows, Dependabot and review templates
```

## Validate without downloading data

Python 3.11 or later is sufficient for repository checks; no third-party Python packages are required.

```bash
python tools/validate_repository.py
python -m unittest discover -s tests -v
python tools/catalog.py summary
python tools/catalog.py taxa --scope fetched
python tools/restore_sources.py --destination /path/to/separate/workspace --dry-run
```

Restoration only writes source files. It does not launch analyses, download data, edit existing files or change running jobs. Inspect legacy absolute paths and environments before execution. Scientific scripts require the external files listed in the inventories.

## Scientific scope

The stored core tree has 66 genome labels. The current complete BenA-CaM-RPB2 tree has 230 tips: 75 genome-derived and 155 reference tips. Older 283-tip and type-focused 108-tip analyses are retained as distinct inference contexts. Matched comparisons prune each tree to the exact same 66 genome labels; pruning a reference-rich tree is not equivalent to inferring a new tree on those genomes.

The current orthogroup screen contains 4,661 groups with one gene model in every core genome. These are cgMLST **candidates**, not a validated scheme. Type-derived assignments remain conditional on the sampled references; the terreus/montoensis core disagreement is open. See the gates before interpreting results as species delineation or epidemiological typing.

## Versions and provenance

The import represents the available files at the 2026-11-20 audit. Snapshot dates and filesystem modification times are evidence labels, not fabricated Git history or proof of an execution date. Identical script copies share one archival object; every observed source location is retained in the version manifest. Earlier edits that were never saved cannot be recovered.

Use the repository URL and a commit hash when citing a specific source version. Publication sources and their accession identifiers are recorded in the inventory. A software license has not been selected by the rights holders; this import does not relicense third-party tools or source material.
