# Reproducibility and environments

## Supported repository operations

Repository utilities use Python >=3.11 and the standard library. CI validates archive hashes, cohort/reference relationships, source-only policy, local documentation links, Python syntax and safe restoration. Unit tests exercise actual Newick pruning/support behavior and refusal of unsafe restore paths. These checks do not run the scientific pipeline.

## Scientific execution

Restore scripts into a separate workspace with the expected hierarchy. Review all absolute paths before running them; historical files retain original paths as provenance. Restore required external inputs from the inventory and verify recorded checksums where available. Missing inputs must cause a visible failure, not an invented result.

The original environment YAML snapshots are indexed in the source manifest. They can contain platform-specific packages and prefixes and are historical evidence rather than a validated portable lockfile. The original `config/config.yaml` has malformed mapping syntax and is archived as evidence only. A valid repository example is provided in `config/project.example.json`; it does not claim that every legacy script consumes it.

Known local figure-runtime versions at publication were NumPy 2.4.3, Matplotlib 3.10.8 and PyMuPDF 1.28.2. A separately vendored Biopython 1.88 was used for independent tree checks. Those vendor directories and binaries are excluded. Other tool versions should be read from their recorded logs/environment snapshots; unknown versions remain unknown.

Do not import historical scripts to discover functions: many execute work at module top level. Inspect them statically or restore and invoke them deliberately. Several figure scripts expect prior output directories or use `__file__` to derive paths, so running directly from `scripts/current/` may target the wrong location.

## Resource rules

Run jobs sequentially where they compete with the existing analyses. Use per-process/cgroup caps, not just a tool's requested thread count. The publication and tree-review utilities do not inspect, kill or resize other jobs. A nominal `-mem` option or `OMP_NUM_THREADS` value is not a universal hard memory ceiling. Review the archived capped launchers and verify limits before executing heavy tools.

## Version changes

Record exact input accession versions, hashes, tool versions, commands, seeds, partition files and exclusion reasons. Write outputs to a new dated version folder outside Git. Preserve independent biological sample identity separately from byte- or sequence-identical assemblies. Use matched cohorts for distance comparisons and disclose changes in reference sampling.
