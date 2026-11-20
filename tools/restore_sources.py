"""Restore archived source without executing it or overwriting existing work."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def safe_child(root: Path, relative: str) -> Path:
    rel = Path(relative)
    if rel.is_absolute() or '..' in rel.parts or ':' in relative or '\\' in relative:
        raise ValueError(f'Unsafe relative path: {relative}')
    target = (root / rel).resolve()
    if not target.is_relative_to(root.resolve()):
        raise ValueError(f'Path escapes destination: {relative}')
    return target

def plan_restore(repo: Path, destination: Path, include_history: bool = False) -> list[dict]:
    repo, destination = repo.resolve(), destination.resolve()
    if destination == repo or destination.is_relative_to(repo):
        raise ValueError('Choose a separate destination outside this repository.')
    records = json.loads((repo / 'manifests/script_versions.json').read_text(encoding='utf-8'))
    plan, seen = [], set()
    for row in records:
        if row['origin'] != 'workspace' or not (include_history or row['latest_workspace_entry']):
            continue
        src = safe_child(repo, row['stored_path'])
        dest = safe_child(destination, row['original_relative_path'])
        if dest in seen:
            raise ValueError(f'Duplicate restore target: {dest}')
        seen.add(dest)
        content = src.read_bytes()
        if hashlib.sha256(content).hexdigest() != row['stored_sha256']:
            raise ValueError(f'Archive hash mismatch: {src}')
        action = 'write'
        if dest.exists():
            if not dest.is_file() or dest.read_bytes() != content:
                raise ValueError(f'Refusing to overwrite existing different content: {dest}')
            action = 'skip-identical'
        plan.append({'source': str(src), 'target': str(dest), 'action': action})
    return plan

def restore(repo: Path, destination: Path, include_history: bool = False, dry_run: bool = True) -> list[dict]:
    plan = plan_restore(repo, destination, include_history)
    if not dry_run:
        for row in plan:
            if row['action'] == 'write':
                target = Path(row['target'])
                target.parent.mkdir(parents=True, exist_ok=True)
                # Exclusive creation also prevents a race from replacing newly created files.
                with target.open('xb') as handle:
                    handle.write(Path(row['source']).read_bytes())
    return plan

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--destination', type=Path, required=True)
    p.add_argument('--include-history', action='store_true')
    p.add_argument('--dry-run', action='store_true')
    a = p.parse_args()
    print(json.dumps(restore(ROOT, a.destination, a.include_history, a.dry_run), indent=2))

if __name__ == '__main__':
    main()
