import random
import tempfile
import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from reconstruct_flye_graph import reconstruct, equivalence


class GraphReconstructionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        root = Path(self.tmp.name)
        random.seed(42)
        self.seq = "".join(random.choices("ACGT", k=6163))
        self.nano = root / "nano.fa"
        self.nano.write_text(">nano\n" + self.seq + "\n")
        self.contigs = {
            "contig_1": self.seq[5273:5946],
            "contig_2": self.seq[5946:] + self.seq[:3167],
            "contig_3": self.seq[3167:5325]
        }
        self.graph = root / "graph.gfa"
        self.graph.write_text(
            "P\tcontig_1\tedge_1+\t*\n"
            "P\tcontig_2\tedge_2+\t*\n"
            "P\tcontig_3\tedge_3+\t*\n"
            "L\tedge_1\t+\tedge_2\t+\t0M\n"
            "L\tedge_2\t+\tedge_3\t+\t0M\n"
            "L\tedge_1\t-\tedge_3\t-\t0M\n"
        )

    def test_three_contig_candidate(self):
        (seq, order, _, overlaps), _ = reconstruct(
            self.contigs, self.graph, self.nano, 6164)
        self.assertEqual(len(seq), 6163)
        self.assertTrue(equivalence(seq, self.seq))
        self.assertEqual(sum(overlaps), 52)
        self.assertEqual(set(order), set(self.contigs))

    def test_missing_link_rejected(self):
        self.graph.write_text(self.graph.read_text().replace(
            "L\tedge_1\t-\tedge_3\t-\t0M\n", ""))
        with self.assertRaises(ValueError):
            reconstruct(self.contigs, self.graph, self.nano, 6164)

    def test_disagreement_rejected(self):
        changed = dict(self.contigs)
        old = changed["contig_2"][0]
        changed["contig_2"] = ("T" if old != "T" else "A") + changed["contig_2"][1:]
        with self.assertRaises(ValueError):
            reconstruct(changed, self.graph, self.nano, 6164)


if __name__ == "__main__":
    unittest.main()
