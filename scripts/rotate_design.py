#!/usr/bin/env python3
"""Make a copy of the design rotated by half its length.

Reads are aligned to a linear FASTA, so reads crossing the FASTA's start/end (the
"origin" of the circle) get split. In the rotated copy that point sits in the middle,
so we can check it separately. The origin's new position is saved for the origin check.
"""
import argparse
from pathlib import Path

from plasmid_utils import read_single_sequence, write_fasta


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    for name in ('design', 'output-fasta', 'origin-file'):
        parser.add_argument('--' + name, required=True)
    args = parser.parse_args()

    design = read_single_sequence(args.design)
    shift = len(design) // 2
    write_fasta(args.output_fasta, 'design_rotated', design[shift:] + design[:shift])
    # position (0-based) in the rotated copy where the original start/end now joins
    Path(args.origin_file).write_text(f'{len(design) - shift}\n')


if __name__ == '__main__':
    main()
