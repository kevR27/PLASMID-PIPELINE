"""Synthetic repeat screening, including the bundled design sequences."""
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from check_repeat_regions import audit, candidate_regions
from plasmid_utils import read_single_sequence


class RepeatScreenTests(unittest.TestCase):
    def test_all_designs_have_screenable_repeat_intervals(self):
        for design_path in sorted((ROOT / 'designs').glob('*.fa')):
            design = read_single_sequence(design_path)
            with self.subTest(design=design_path.name):
                self.assertTrue(candidate_regions(design))
                rows = audit(design, design, sample='control')
                self.assertTrue(all(r['status'] == 'SEQUENCE_MATCH' for r in rows))

    def test_102_bp_deletion_in_right_repeat_is_detected(self):
        path = ROOT / 'designs' / 'TAL MERRF 024130 right clone 12-YES.fa'
        design = read_single_sequence(path)
        start, end = candidate_regions(design)[0]
        mutated = design[:start + 200] + design[start + 302:]
        rows = audit(design, mutated, sample='deletion')
        self.assertEqual(rows[0]['status'], 'SEQUENCE_DIFFERENT')
        self.assertEqual(rows[0]['assembly_bp'], end - start - 102)
        self.assertGreaterEqual(rows[0]['gap_bases'], 102)


if __name__ == '__main__':
    unittest.main()
