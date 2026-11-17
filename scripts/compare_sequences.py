#!/usr/bin/env python3
"""Compare two plasmid sequences base by base (global alignment).

Used for: assembly vs design, Flye vs NANOPlasty, and barcode vs barcode.
Both sequences should already start at the same place (see orient_assembly.py).

Writes
  --summary      one row: lengths, matching/substituted/inserted/deleted bases,
                 identity %, number of differences, identical YES/NO
  --differences  one row per difference (SNV = 1 base changed, MNV = several bases,
                 INS = bases only in the second sequence, DEL = bases only in the first)

Column names use the labels you pass, e.g. --reference-label design --other-label assembly
gives columns "design_sequence" and "assembly_sequence".
"""
import argparse

from Bio import Align

from plasmid_utils import read_single_sequence, write_table


def compare(reference, other):
    """Globally align `other` to `reference`. Returns (summary_dict, list_of_difference_dicts)."""
    aligner = Align.PairwiseAligner(mode='global', match_score=2, mismatch_score=-4,
                                    open_gap_score=-7, extend_gap_score=-2)
    alignment = aligner.align(reference, other)[0]
    ref_index, other_index = alignment.indices   # -1 marks a gap

    matching = substituted = inserted = deleted = 0
    differences = []
    i = 0
    while i < len(ref_index):
        r, o = int(ref_index[i]), int(other_index[i])
        j = i + 1
        if r >= 0 and o >= 0:
            if reference[r] == other[o]:
                matching += 1
                i += 1
                continue
            # a run of consecutive mismatches
            while (j < len(ref_index) and ref_index[j] >= 0 and other_index[j] >= 0
                   and reference[int(ref_index[j])] != other[int(other_index[j])]):
                j += 1
            ref_part = reference[r:int(ref_index[j - 1]) + 1]
            other_part = ''.join(other[int(x)] for x in other_index[i:j])
            substituted += len(ref_part)
            differences.append(dict(type='SNV' if len(ref_part) == 1 else 'MNV', position_1based=r + 1,
                                    reference_sequence=ref_part, other_sequence=other_part,
                                    length_bp=len(ref_part)))
        elif r < 0:                                   # bases only in `other`: insertion
            while j < len(ref_index) and ref_index[j] < 0:
                j += 1
            inserted_part = ''.join(other[int(x)] for x in other_index[i:j])
            inserted += len(inserted_part)
            # position = the reference base BEFORE the insertion (0 = before the first base)
            after_base = int(ref_index[i - 1]) + 1 if i and ref_index[i - 1] >= 0 else 0
            differences.append(dict(type='INS', position_1based=after_base, reference_sequence='',
                                    other_sequence=inserted_part, length_bp=len(inserted_part)))
        else:                                         # bases only in `reference`: deletion
            while j < len(ref_index) and other_index[j] < 0:
                j += 1
            deleted_part = ''.join(reference[int(x)] for x in ref_index[i:j])
            deleted += len(deleted_part)
            differences.append(dict(type='DEL', position_1based=r + 1, reference_sequence=deleted_part,
                                    other_sequence='', length_bp=len(deleted_part)))
        i = j

    total = matching + substituted + inserted + deleted
    summary = dict(
        reference_bp=len(reference), other_bp=len(other), matching_bp=matching,
        substituted_bp=substituted, inserted_bp=inserted, deleted_bp=deleted,
        identity_percent=round(100 * matching / total, 5) if total else 0,
        reference_covered_percent=round(100 * (matching + substituted) / len(reference), 5),
        differences=len(differences),
        identical='YES' if matching == len(reference) == len(other) and not differences else 'NO')
    return summary, differences


def rename(row, reference_label, other_label):
    """Swap the generic 'reference'/'other' in column names for the real labels."""
    return {key.replace('reference', reference_label).replace('other', other_label): value
            for key, value in row.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    for name in ('reference', 'other', 'summary', 'differences'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--reference-label', default='reference')
    parser.add_argument('--other-label', default='other')
    args = parser.parse_args()

    summary, differences = compare(read_single_sequence(args.reference), read_single_sequence(args.other))
    summary = rename(summary, args.reference_label, args.other_label)
    write_table(args.summary, [summary], list(summary))
    columns = rename(dict.fromkeys(['type', 'position_1based', 'reference_sequence', 'other_sequence', 'length_bp']),
                     args.reference_label, args.other_label)
    write_table(args.differences, [rename(d, args.reference_label, args.other_label) for d in differences],
                list(columns))


if __name__ == '__main__':
    main()
