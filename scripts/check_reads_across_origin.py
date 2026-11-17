#!/usr/bin/env python3
"""Count reads that run across the design's artificial start/end point.

Uses the reads aligned to the rotated design (where that point is in the middle).
Reads crossing it show the two ends of the design really are joined in your plasmid
(circular). A long concatemer or chimeric read can also cross, so this is supporting
evidence, not proof.
"""
import argparse
from pathlib import Path

import pysam

from plasmid_utils import write_table


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    for name in ('bam', 'origin-file', 'output'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--min-flank', type=int, default=50, help='bp required on each side of the point')
    parser.add_argument('--min-mapq', type=int, default=20)
    args = parser.parse_args()

    origin = int(Path(args.origin_file).read_text().strip())
    primary_alignments = crossing = 0
    with pysam.AlignmentFile(args.bam, 'rb') as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate:
                continue
            primary_alignments += 1
            covers_point = (read.reference_start <= origin - args.min_flank
                            and read.reference_end >= origin + args.min_flank)
            if covers_point and read.mapping_quality >= args.min_mapq:
                crossing += 1

    row = dict(origin_position_in_rotated_design=origin, min_flank_bp=args.min_flank,
               min_mapping_quality=args.min_mapq, primary_alignments=primary_alignments,
               alignments_crossing_origin=crossing)
    write_table(args.output, [row], list(row))


if __name__ == '__main__':
    main()
