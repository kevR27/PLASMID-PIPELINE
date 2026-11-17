#!/usr/bin/env python3
"""For every assembly-vs-design difference, count how many raw reads agree with the assembly.

Looks at the reads aligned to the design:
  SNV : reads showing the assembly base at that position (plus counts of A/C/G/T there)
  INS : reads with the same inserted bases right after the position
  DEL : reads with a deletion of the same length right after the previous base
MNVs and differences at the artificial circular start/end are only marked for manual review.
These are quick checks, not formal variant calls: look at indels in IGV.
"""
import argparse

import pysam

from plasmid_utils import read_table, write_table

BASES = 'ACGT'
COLUMNS = ['type', 'position_1based', 'design_sequence', 'assembly_sequence',
           'reads_covering_site', 'reads_supporting_assembly', 'support_percent',
           'reads_with_A', 'reads_with_C', 'reads_with_G', 'reads_with_T', 'note']


def read_is_usable(read, min_mapq):
    return not (read.is_unmapped or read.is_secondary or read.is_supplementary
                or read.is_duplicate or read.mapping_quality < min_mapq)


def check_difference(bam, difference, min_mapq, min_base_quality):
    """Return (covering, supporting, percent, base_counts, note) for one difference."""
    kind = difference['type']
    pos = int(difference['position_1based'])
    design_part = difference['design_sequence'].upper()
    assembly_part = difference['assembly_sequence'].upper()
    no_counts = ['NA'] * 4

    if kind not in ('SNV', 'INS', 'DEL'):
        return 'NA', 'NA', 'NA', no_counts, 'MNV_check_manually'
    if (kind == 'INS' and pos == 0) or (kind == 'DEL' and pos == 1):
        return 'NA', 'NA', 'NA', no_counts, 'at_circular_origin_check_manually'

    # 0-based position of the base we look at: the SNV base itself, or for an indel
    # the base just before it (CIGAR indels are attached to the preceding base)
    anchor = pos - 1 if kind in ('SNV', 'INS') else pos - 2
    covering = supporting = 0
    base_counts = dict.fromkeys(BASES, 0)

    for column in bam.pileup(bam.references[0], anchor, anchor + 1, truncate=True, max_depth=1000000,
                             min_base_quality=0, min_mapping_quality=0, stepper='nofilter'):
        if column.reference_pos != anchor:
            continue
        for pileup_read in column.pileups:
            read = pileup_read.alignment
            q = pileup_read.query_position
            if not read_is_usable(read, min_mapq) or q is None or pileup_read.is_refskip or pileup_read.is_del:
                continue
            if read.query_qualities is not None and read.query_qualities[q] < min_base_quality:
                continue
            covering += 1
            base = read.query_sequence[q].upper()
            if base in base_counts:
                base_counts[base] += 1
            if kind == 'SNV':
                supporting += base == assembly_part
            elif kind == 'INS':
                after = read.query_sequence[q + 1:q + 1 + len(assembly_part)].upper()
                if pileup_read.indel == len(assembly_part) and after == assembly_part:
                    supporting += 1
            else:   # DEL
                supporting += pileup_read.indel == -len(design_part)

    percent = round(100 * supporting / covering, 4) if covering else 'NA'
    counts = [base_counts[b] for b in BASES] if kind == 'SNV' else no_counts
    note = 'deletion_length_only' if kind == 'DEL' else 'quick_pileup_check'
    return covering, supporting, percent, counts, note


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    for name in ('bam', 'differences', 'output'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--min-mapq', type=int, default=20)
    parser.add_argument('--min-base-quality', type=int, default=10)
    args = parser.parse_args()

    rows = []
    with pysam.AlignmentFile(args.bam, 'rb') as bam:
        for difference in read_table(args.differences):
            covering, supporting, percent, counts, note = check_difference(
                bam, difference, args.min_mapq, args.min_base_quality)
            row = {key: difference[key] for key in COLUMNS[:4]}
            row.update(reads_covering_site=covering, reads_supporting_assembly=supporting,
                       support_percent=percent, note=note)
            row.update({f'reads_with_{b}': c for b, c in zip(BASES, counts)})
            rows.append(row)
    write_table(args.output, rows, COLUMNS)


if __name__ == '__main__':
    main()
