# Available source versions

`manifests/script_versions.json` maps every observed source location to a SHA-256-addressed object under `archive/objects/`. It records origin, original relative path, content hashes, snapshot label, source mtime and any export transformation. Identical bytes are stored once. `scripts/current/` is the latest observed top-level workspace view, not necessarily the most recently modified snapshot of every historical entry point.

`manifests/analysis_versions.json` lists dated analysis folders. `manifests/tree_versions.json` identifies the core, complete 230-tip, historical 283-tip and type-focused 108-tip tree sources by hash; tree payloads remain outside Git. The separately inferred exact-66 three-gene combination is catalogued through its script and analysis-version provenance.

Original scripts may have absolute Windows or WSL paths, sibling imports and output paths based on `__file__`. Moving them between categories would silently change those paths. The source archive therefore uses a category index and a flat current view, with restoration to the original workspace-relative hierarchy when needed.

```bash
python tools/restore_sources.py --destination /outside/repository/analysis --dry-run
python tools/restore_sources.py --destination /outside/repository/analysis
python tools/restore_sources.py --destination /outside/repository/history --include-history --dry-run
```

The default restores current workspace sources only. `--include-history` also restores workspace-relative saved versions. Original-project command snippets are retained in the archive and catalog, not automatically executed or installed. Existing different files, path traversal, symlinks and destinations inside this repository are refused. Identical existing files are skipped.

Notebook exports contain code cells only; they are historical text, not automatically valid stand-alone Python programs. Language-tagged executable fences are extracted from the mixed lab notebook. Full notebook outputs and mixed analytical narratives remain local. Original and transformed hashes distinguish these exports from byte-identical source snapshots.

No commits have been backdated and no unavailable intermediate edits are invented. Future code changes should create ordinary Git commits and update the manifest intentionally; do not silently replace an archived object.
