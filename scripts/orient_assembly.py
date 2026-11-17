#!/usr/bin/env python3
"""Step 3: flip and rotate an assembled plasmid so it starts where the design starts.

An assembler can return the circular plasmid on either strand and starting anywhere.
To compare it base-by-base with the design we need the same strand and the same start.
We find a 31-bp piece from the start of the design that occurs exactly once in the
assembly (tried on both strands) and rotate the assembly to that point.

If the assembly has several contigs the LONGEST one is used and the sample is marked
for review. If no unique 31-bp piece is found the sequence is left as it is and marked
COULD_NOT_ORIENT. The original assembly file is never changed.
"""
import argparse

from Bio import SeqIO

from plasmid_utils import (is_same_circular_sequence, read_single_sequence, reverse_complement,
                           write_fasta, write_table)

ANCHOR_LENGTH = 31
ANCHOR_SEARCH_WINDOW = 500   # look for the anchor among the first 500 design positions


def orient(design, sequence):
    """Return (oriented_sequence, strand, rotation_start_1based, anchor_in_design_1based, status)."""
    for strand, candidate in (('+', sequence), ('-', reverse_complement(sequence))):
        for design_pos in range(min(len(design) - ANCHOR_LENGTH + 1, ANCHOR_SEARCH_WINDOW)):
            anchor = design[design_pos:design_pos + ANCHOR_LENGTH]
            if set(anchor) - set('ACGT'):
                continue
            if design.count(anchor) != 1 or candidate.count(anchor) != 1:
                continue
            found_at = candidate.find(anchor)
            rotation_start = (found_at - design_pos) % len(candidate)
            rotated = candidate[rotation_start:] + candidate[:rotation_start]
            return rotated, strand, rotation_start + 1, design_pos + 1, 'OK'
    return sequence, '.', 1, 0, 'COULD_NOT_ORIENT'


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    for name in ('assembly', 'design', 'name', 'output-fasta', 'output-report'):
        parser.add_argument('--' + name, required=True)
    args = parser.parse_args()

    contigs = list(SeqIO.parse(args.assembly, 'fasta'))
    if not contigs:
        raise ValueError(f'No sequences in {args.assembly}')
    design = read_single_sequence(args.design)
    longest = max(contigs, key=lambda c: len(c.seq))
    sequence = str(longest.seq).upper()

    oriented, strand, rotation_start, anchor, status = orient(design, sequence)
    write_fasta(args.output_fasta, f'{args.name}_oriented_to_design', oriented)

    report = dict(
        sample=args.name, number_of_contigs=len(contigs),
        contig_lengths=';'.join(str(len(c.seq)) for c in contigs),
        longest_contig_length=len(sequence), design_length=len(design),
        strand=strand, status=status,
        identical_circular_sequence='YES' if is_same_circular_sequence(design, sequence) else 'NO',
        needs_review='YES' if len(contigs) != 1 or status != 'OK' else 'NO')
    write_table(args.output_report, [report], list(report))


if __name__ == '__main__':
    main()
