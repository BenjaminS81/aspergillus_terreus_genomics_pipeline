# Contributing

Use a focused branch and pull request. Describe the concrete scientific or operational problem, resulting behavior, validation and material limitations. Preserve meaningful commits; do not create empty commits, self-reviews, artificial issues or extra PRs to inflate activity or achievements.

Before a merge:

1. Run `python tools/validate_repository.py` and `python -m unittest discover -s tests -v`.
2. Confirm no data payloads, credentials, notebook outputs or generated figures are staged.
3. Update manifests and documentation when an input version, method or gate changes.
4. Keep original identifiers and taxonomic evidence separate from inferred assignments.
5. Use the same leaves for tree comparisons; explain missing data and copy ambiguity.
6. Preserve immutable archive objects. New versions get new hashes and new manifest records.

Use explicit seeds and bounded resources where supported. Avoid invoking heavy scientific workflows in repository CI. A passing repository check is not a substitute for scientific review. License changes require agreement from the relevant rights holders.
