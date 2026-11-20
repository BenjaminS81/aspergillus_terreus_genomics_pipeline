import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

restore = module('restore', ROOT / 'tools/restore_sources.py')
guard = module('guard', ROOT / 'tools/validate_repository.py')
tree = module('tree', ROOT / 'scripts/current/terreus_tree_helpers.py')

class RestoreTests(unittest.TestCase):
    def fixture(self, root):
        repo = root / 'repo'
        (repo / 'archive').mkdir(parents=True)
        (repo / 'manifests').mkdir()
        content = b'print("source only")\n'
        (repo / 'archive/source.py').write_bytes(content)
        row = {'origin': 'workspace', 'latest_workspace_entry': True, 'stored_path': 'archive/source.py', 'original_relative_path': 'a.py', 'stored_sha256': hashlib.sha256(content).hexdigest()}
        (repo / 'manifests/script_versions.json').write_text(json.dumps([row]))
        return repo, row

    def test_dry_run_and_idempotent_restore(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); repo, row = self.fixture(root); dest = root / 'external'
            self.assertEqual(restore.restore(repo, dest)[0]['action'], 'write')
            self.assertFalse(dest.exists())
            restore.restore(repo, dest, dry_run=False)
            self.assertEqual(restore.restore(repo, dest)[0]['action'], 'skip-identical')

    def test_conflicting_file_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); repo, row = self.fixture(root); dest = root / 'external'; dest.mkdir()
            (dest / 'a.py').write_text('user work')
            with self.assertRaises(ValueError):restore.restore(repo, dest, dry_run=False)
            self.assertEqual((dest / 'a.py').read_text(), 'user work')

    def test_traversal_hash_mismatch_and_in_repo_destination_refused(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); repo, row = self.fixture(root)
            with self.assertRaises(ValueError):restore.restore(repo, repo / 'nested')
            for unsafe in ('../escape.py', '/absolute.py', 'C:/escape.py', 'a\\b.py'):
                with self.assertRaises(ValueError):restore.safe_child(root, unsafe)
            (repo / 'archive/source.py').write_text('tampered')
            with self.assertRaises(ValueError):restore.restore(repo, root / 'external')

class PublicationGuardTests(unittest.TestCase):
    def test_sequence_and_binary_payload_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); p = root / 'payload.txt'; p.write_text('ACGT' * 70)
            self.assertTrue(guard.publication_issues(p, root))
            p.write_bytes(b'\x00\xff')
            self.assertTrue(guard.publication_issues(p, root))

    def test_token_pattern_rejected_without_echoing_token(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); p = root / 'config.txt'; value = 'gh' + 'p_' + 'a' * 36; p.write_text(value)
            issues = guard.publication_issues(p, root)
            self.assertTrue(issues)
            self.assertNotIn(value, '\n'.join(issues))

class TreeHelperTests(unittest.TestCase):
    def test_pruning_preserves_path_lengths_and_unique_tips(self):
        r = tree.parse_newick('((a:1,b:2)90/99:3,c:4,d:5);')
        p = tree.prune(r, {'a', 'c', 'd'})
        self.assertEqual({n.name for n in p.leaves()}, {'a', 'c', 'd'})
        self.assertEqual(next(n.length for n in p.leaves() if n.name == 'a'), 4)
        self.assertEqual({n.name for n in r.leaves()}, {'a', 'b', 'c', 'd'})

    def test_support_threshold_and_collapsed_support(self):
        self.assertTrue(tree.strong(tree.Node('80/95')))
        self.assertFalse(tree.strong(tree.Node('79.9/100')))
        self.assertFalse(tree.strong(tree.Node('100')))
        r = tree.parse_newick('(((a:1,b:1)90/99:2,c:1)95/100:3,d:1);')
        p = tree.prune(r, {'a', 'b', 'd'})
        self.assertEqual(p.children[0].name, '')
        self.assertEqual(p.children[0].length, 5)

if __name__ == '__main__':
    unittest.main()
