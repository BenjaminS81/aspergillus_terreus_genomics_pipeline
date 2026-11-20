# Security and source-only publication policy

Report vulnerabilities privately through the repository's GitHub security reporting feature if enabled, or contact the repository owner privately. Do not disclose credentials, private samples or exploitable details in public issues.

CI has read-only repository permissions, pinned third-party action commits, no persisted checkout credentials and no `pull_request_target` execution. It checks the tracked tree for credential patterns, private keys, large files, binary payloads, sequence-data extensions, embedded long nucleotide strings and unsafe workflow configuration. This is a targeted publication guard, not a guarantee that all security weaknesses are detectable.

Historical scripts are archived research source. They may call external programs, make downloads, contain absolute paths or run at import time. They are not executed by CI. Review commands before running them. External inputs must be verified independently, and a download must not overwrite an existing validated input silently.

No tokens or credentials are needed for local checks. GitHub authentication belongs in the platform credential mechanism, never a script, manifest, notebook or workflow file. Downloaded papers and third-party binaries are excluded. The repository does not grant permission to redistribute third-party data or software.
