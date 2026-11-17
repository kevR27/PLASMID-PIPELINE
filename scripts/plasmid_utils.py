"""Small helpers shared by several scripts (kept here so each one exists only once)."""
import csv
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


def read_single_sequence(fasta_path):
    """Return the sequence (upper case) of a FASTA file that must hold exactly one record."""
    records = list(SeqIO.parse(fasta_path, 'fasta'))
    if len(records) != 1:
        raise ValueError(f'Expected exactly one FASTA record in {fasta_path}, found {len(records)}')
    return str(records[0].seq).upper()


def reverse_complement(sequence):
    return str(Seq(sequence).reverse_complement())


def is_same_circular_sequence(sequence_a, sequence_b):
    """True if the two sequences are the same circular molecule (any rotation, either strand)."""
    if len(sequence_a) != len(sequence_b):
        return False
    doubled = sequence_a + sequence_a
    return sequence_b in doubled or reverse_complement(sequence_b) in doubled


def read_table(path):
    """Read a tab-separated file with a header into a list of dicts."""
    with open(path, newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_table(path, rows, fields):
    """Write a list of dicts as a tab-separated file (creates the folder if needed)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def write_fasta(path, name, sequence, width=80):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w') as handle:
        handle.write(f'>{name}\n')
        for start in range(0, len(sequence), width):
            handle.write(sequence[start:start + width] + '\n')
