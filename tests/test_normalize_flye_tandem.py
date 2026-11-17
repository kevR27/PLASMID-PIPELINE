"""Regression tests for the Flye tandem normalization script.

Run: python -m unittest discover -s tests -v
"""
import csv
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts' / 'normalize_flye_tandem.py'


class NormalizationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        random.seed(77)
        cls.design = ''.join(random.choices('ACGT', k=6164))
        cls.unit = cls.design[613:] + cls.design[:613]

    def run_case(self, contig, expect_success=True):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'design.fasta').write_text('>design\n' + self.design + '\n')
            (root/'raw.fasta').write_text('>contig\n' + contig + '\n')
            args = [sys.executable, str(SCRIPT), '--assembly', str(root/'raw.fasta'),
                    '--design', str(root/'design.fasta'), '--sample', 'sample',
                    '--output-fasta', str(root/'normalized.fasta'),
                    '--report', str(root/'normalization.tsv'),
                    '--boundaries', str(root/'boundaries.tsv'),
                    '--pairwise', str(root/'pairwise.tsv'),
                    '--copies-fasta', str(root/'copies.fasta')]
            p = subprocess.run(args, capture_output=True, text=True)
            if not expect_success:
                self.assertNotEqual(p.returncode, 0)
                self.assertFalse((root/'normalized.fasta').exists())
                return p.stderr
            self.assertEqual(p.returncode, 0, p.stderr)
            with (root/'normalized.fasta').open() as handle:
                seq = str(next(SeqIO.parse(handle, 'fasta')).seq)
            with (root/'normalization.tsv').open() as h:
                report = next(csv.DictReader(h, delimiter='\t'))
            return seq, report

    def test_single_contig_passes_unmodified(self):
        seq, rep = self.run_case(self.unit)
        self.assertEqual(seq, self.unit)
        self.assertEqual(rep['estimated_copies'], '1')

    def test_two_copies_with_indel(self):
        second = self.unit[:3100] + self.unit[3101:]
        seq, rep = self.run_case(self.unit + second)
        self.assertEqual(seq, self.design)
        self.assertEqual(rep['estimated_copies'], '2')
        self.assertGreaterEqual(float(rep['min_pairwise_identity']), .999)

    def test_reverse_strand(self):
        second = self.unit[:3100] + self.unit[3101:]
        seq, rep = self.run_case(str(Seq(self.unit + second).reverse_complement()))
        self.assertEqual(seq, self.design)
        self.assertEqual(rep['estimated_copies'], '2')
        self.assertEqual(rep['raw_orientation_for_detection'], '-')

    def test_three_units(self):
        seq, rep = self.run_case(self.unit * 3)
        self.assertEqual(seq, self.design)
        self.assertEqual(rep['estimated_copies'], '3')

    def test_shared_variant_is_preserved_not_overwritten_by_design(self):
        # Both copies carry an identical SNV relative to the design.
        mutated = list(self.unit)
        mutated[1500] = 'A' if mutated[1500] != 'A' else 'C'
        mutated = ''.join(mutated)
        observed, rep = self.run_case(mutated * 2)
        expected = mutated[-613:] + mutated[:-613]
        self.assertEqual(observed, expected)
        self.assertNotEqual(observed, self.design)
        self.assertEqual(rep['estimated_copies'], '2')

    def test_divergent_units_rejected(self):
        divergent = list(self.unit)
        for j in range(120):
            pos = 200 + 40*j
            divergent[pos] = 'C' if divergent[pos] != 'C' else 'A'
        stderr = self.run_case(self.unit + ''.join(divergent), expect_success=False)
        self.assertIn('Refusing to collapse', stderr)


if __name__ == '__main__':
    unittest.main()
