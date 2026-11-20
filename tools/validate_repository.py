"""Data-free integrity and publication guard; never execute archived workflows."""
from __future__ import annotations
import argparse
import ast
import csv
import hashlib
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DENIED = {'.fa', '.fasta', '.fna', '.faa', '.fastq', '.fq', '.gz', '.zip', '.bam', '.sam', '.cram', '.vcf', '.bcf', '.gff', '.gff3', '.treefile', '.contree', '.ufboot', '.pdf', '.png', '.jpg', '.xlsx', '.ipynb', '.pyc'}
SECRETS = [
    re.compile(r'gh[pousr]_[A-Za-z0-9]{30,}'),
    re.compile(r'github_pat_[A-Za-z0-9_]{40,}'),
    re.compile(r'AKIA[A-Z0-9]{16}'),
    re.compile(r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----'),
    re.compile(r'(?i)(?:api[_-]?key|access[_-]?token|password)\s*[:=]\s*[\"\x27]([A-Za-z0-9_+/=-]{24,})[\"\x27]')
]

def repository_files(root: Path) -> list[Path]:
    top = subprocess.run(['git', '-C', str(root), 'rev-parse', '--show-toplevel'], capture_output=True, text=True)
    if top.returncode == 0 and Path(top.stdout.strip()).resolve() == root.resolve():
        res = subprocess.run(['git', '-C', str(root), 'ls-files', '-z'], capture_output=True, check=True)
        paths = [root / x.decode('utf-8') for x in res.stdout.split(b'\0') if x]
        if paths:
            return paths
    return [p for p in root.rglob('*') if p.is_file() and not {'.git', '__pycache__', '.pytest_cache'} & set(p.parts)]

def publication_issues(path: Path, root: Path) -> list[str]:
    rel = path.relative_to(root).as_posix()
    issues = []
    if path.is_symlink():
        issues.append(f'{rel}: symlink not allowed')
    if path.suffix.lower() in DENIED or path.name == '.env':
        issues.append(f'{rel}: forbidden payload extension')
    raw = path.read_bytes()
    if len(raw) > 5 * 1024 * 1024:
        issues.append(f'{rel}: exceeds 5 MiB source/metadata limit')
    try:
        text = raw.decode('utf-8-sig')
    except UnicodeDecodeError:
        return issues + [f'{rel}: non-UTF8/binary payload']
    if '\x00' in text:
        issues.append(f'{rel}: binary NUL bytes')
    if re.search(r'(?<![A-Za-z])[ACGTNacgtn]{200,}(?![A-Za-z])', text):
        issues.append(f'{rel}: embedded long nucleotide payload')
    for i, pattern in enumerate(SECRETS):
        if pattern.search(text):
            issues.append(f'{rel}: possible credential pattern {i}; value suppressed')
    if rel.startswith('.github/workflows/'):
        if 'pull_request_target' in text:
            issues.append(f'{rel}: privileged PR trigger is prohibited')
        for action in re.findall(r'uses:\s*([^\s#]+)', text):
            if not re.fullmatch(r'[\w./-]+@[0-9a-f]{40}', action):
                issues.append(f'{rel}: action is not SHA-pinned')
        if 'contents: read' not in text or 'persist-credentials: false' not in text:
            issues.append(f'{rel}: missing least-privilege checkout policy')
    return issues

def validate(root: Path, security_only: bool = False) -> dict:
    root = root.resolve()
    files = repository_files(root)
    problems = [issue for p in files for issue in publication_issues(p, root)]
    if not security_only:
        citation = (root / 'CITATION.cff').read_text(encoding='utf-8')
        for field in ('authors', 'cff-version', 'message', 'title'):
            if not re.search(r'^' + re.escape(field) + r':', citation, re.M):
                problems.append(f'CITATION.cff: missing required field {field}')
        records = json.loads((root / 'manifests/script_versions.json').read_text(encoding='utf-8'))
        summary = json.loads((root / 'manifests/publication_summary.json').read_text(encoding='utf-8'))
        if summary['script_source_records'] != len(records):
            problems.append('Script source count differs from publication summary')
        objects = set()
        for r in records:
            p = (root / r['stored_path']).resolve()
            if not p.is_relative_to(root) or not p.is_file():
                problems.append(f'Missing/unsafe archive object: {r["stored_path"]}')
                continue
            if hashlib.sha256(p.read_bytes()).hexdigest() != r['stored_sha256']:
                problems.append(f'Archive hash mismatch: {r["stored_path"]}')
            objects.add(r['stored_path'])
            if r['latest_workspace_entry']:
                current = root / 'scripts/current' / Path(r['original_relative_path']).name
                if not current.is_file() or current.read_bytes() != p.read_bytes():
                    problems.append(f'Current source differs from manifest: {current.name}')
        actual_objects = {p.relative_to(root).as_posix() for p in files if p.is_relative_to(root / 'archive/objects')}
        if actual_objects != objects:
            problems.append('Archived object list and source-version manifest differ')
        if summary['distinct_archived_source_objects'] != len(objects):
            problems.append('Distinct archive count differs from publication summary')
        for p in files:
            if p.suffix == '.json':
                try:
                    json.loads(p.read_text(encoding='utf-8-sig'))
                except (ValueError, UnicodeError) as e:
                    problems.append(f'{p.relative_to(root)}: invalid JSON: {e}')
            if p.suffix == '.py':
                try:
                    ast.parse(p.read_text(encoding='utf-8-sig'))
                except (SyntaxError, UnicodeError) as e:
                    problems.append(f'{p.relative_to(root)}: invalid Python: {e}')
            if p.suffix == '.md':
                for target in re.findall(r'(?<!!)\[[^\]]+\]\(([^)]+)\)', p.read_text(encoding='utf-8')):
                    if '://' in target or target.startswith('#'):
                        continue
                    target = target.split('#')[0].strip('<>')
                    if target and not (p.parent / target).exists():
                        problems.append(f'{p.relative_to(root)}: broken local link {target}')
        with (root / 'inventory/genome_cohort.tsv').open(encoding='utf-8', newline='') as f:
            cohort = list(csv.DictReader(f, delimiter='\t'))
        if len({r['sample'] for r in cohort}) != len(cohort):
            problems.append('Duplicate genome cohort identifiers')
        if sum(r['core66'] == 'True' for r in cohort) != 66:
            problems.append('Core cohort does not contain 66 labels')
        with (root / 'inventory/reference_strains.tsv').open(encoding='utf-8', newline='') as f:
            refs = list(csv.DictReader(f, delimiter='\t'))
        if sum(r['in_complete230'] == 'True' for r in refs) != 155:
            problems.append('Complete-tree reference cohort does not contain 155 tips')
        taxa = json.loads((root / 'inventory/taxon_evidence.json').read_text(encoding='utf-8'))
        names = (root / 'inventory/species_and_infraspecific_names.txt').read_text(encoding='utf-8').splitlines()
        if names != [r['taxon_name'] for r in taxa] or len(set(names)) != len(names):
            problems.append('Taxon list and evidence table differ')
        with (root / 'inventory/fetched_marker_accessions.tsv').open(encoding='utf-8', newline='') as f:
            accessions = list(csv.DictReader(f, delimiter='\t'))
        if len(accessions) != summary['fetched_marker_accessions'] or len({r['accession'] for r in accessions}) != len(accessions):
            problems.append('Fetched-accession count or uniqueness check failed')
        if any(r['sequence_fetched'] and not any(k in r['evidence'] for k in ('marker_sequence_fetched', 'genome_download_verified', 'historical_sequence_fetched', 'protein_reference_fetched')) for r in taxa):
            problems.append('Fetched taxon lacks a sequence retrieval evidence class')
    if problems:
        raise ValueError('\n'.join(problems))
    return {'status': 'passed', 'files_checked': len(files), 'security_only': security_only}

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--security-only', action='store_true')
    a = p.parse_args()
    print(json.dumps(validate(ROOT, a.security_only), indent=2))

if __name__ == '__main__':
    main()
