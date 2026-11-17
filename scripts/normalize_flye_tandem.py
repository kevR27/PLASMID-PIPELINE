#!/usr/bin/env python3
"""Conservatively detect and normalize multi-unit Flye plasmid contigs.

This script does NOT prove a monomeric biological plasmid: a concatemer-like
Flye contig can be an assembly representation or reflect real tandem molecules.

Requires Biopython (Bio.SeqIO and Bio.Align). The designed plasmid is used only
for expected length and to locate unique short sequence anchors. Sequence bases
in the output are taken from a single observed Flye unit, never copied from the
reference. Disagreement among units is reported and not resolved by reference.
"""

import argparse
import csv
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq


def read_one(path):
    records = list(SeqIO.parse(path, 'fasta'))
    if len(records) != 1:
        raise ValueError(f"Expected exactly one FASTA sequence in {path}; found {len(records)}")
    seq = str(records[0].seq).upper()
    if not seq or set(seq) - set('ACGTN'):
        raise ValueError(f"Invalid or empty DNA FASTA in {path}")
    return records[0].id, seq


def write_fasta(path, records):
    with open(path, 'w') as out:
        for name, seq in records:
            out.write(f">{name}\n")
            for pos in range(0, len(seq), 80):
                out.write(seq[pos:pos + 80] + '\n')


def write_tsv(path, rows, fields):
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def positions(haystack, needle):
    result = []
    at = haystack.find(needle)
    while at >= 0:
        result.append(at)
        at = haystack.find(needle, at + 1)
    return result


def candidate_anchor(reference, raw, copies, k, length_tolerance):
    """Find one reference-unique kmer appearing once per plasmid unit.

    The raw contig can begin at any circular origin and in either orientation.
    An anchor must also recur at nearly the expected unit periodicity, including
    the interval across the raw FASTA end and beginning.
    """
    ref_circular = reference + reference[:k - 1]
    best = None
    for orientation, oriented_raw in [('+', raw), ('-', str(Seq(raw).reverse_complement()))]:
        # Scan all offsets, not merely a sparse selection: variants may affect
        # the first few plausible anchors.
        for refpos in range(len(reference)):
            anchor = ref_circular[refpos:refpos + k]
            if 'N' in anchor or len(positions(ref_circular, anchor)) != 1:
                continue
            occurrences = positions(oriented_raw, anchor)
            if len(occurrences) != copies:
                continue
            periods = [occurrences[i+1] - occurrences[i] for i in range(copies - 1)]
            periods += [len(oriented_raw) - occurrences[-1] + occurrences[0]]
            deviations = [abs(x-len(reference)) for x in periods]
            if any(x > length_tolerance for x in deviations):
                continue
            # Prefer near-perfect periodicity, then lowest reference offset.
            score = (max(deviations), sum(deviations), refpos, 0 if orientation == '+' else 1)
            if best is None or score < best[0]:
                best = (score, orientation, oriented_raw, refpos, anchor, occurrences, periods)
            if score[0] == 0 and score[1] == 0 and refpos == 0:
                return best
    return best


def pairwise_counts(a, b, aligner):
    alignment = aligner.align(a, b)[0]
    c = alignment.counts()
    aligned_total = c.identities + c.mismatches + c.gaps
    return {
        'identities': c.identities, 'mismatches': c.mismatches,
        'gap_bases': c.gaps,
        'identity': c.identities / aligned_total if aligned_total else 0.0,
        'differences': c.mismatches + c.gaps,
        'aligned_columns': aligned_total,
    }


def run(args):
    raw_id, raw = read_one(args.assembly)
    _, reference = read_one(args.design)
    length = len(reference)
    n = max(1, round(len(raw) / length))
    if n > args.max_copies or abs(len(raw) / length - n) > args.copy_number_tolerance:
        raise ValueError(f'Raw contig length {len(raw)} is not near 1..{args.max_copies} '
                         f'integral units of design length {length}. Do not automatically collapse it.')

    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -3
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1

    output = Path(args.output_fasta)
    output.parent.mkdir(parents=True, exist_ok=True)
    for path in [args.report, args.boundaries, args.pairwise, args.copies_fasta]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    if n == 1:
        # A single length-consistent unit can continue; this does not prove it is circular.
        oriented_raw = raw
        intervals = [{'copy': 1, 'start_0based': 0, 'end_0based_exclusive': len(raw),
                      'wraps_raw_origin': 'no', 'length_bp': len(raw)}]
        units = [raw]
        pair_rows = []
        strand, refpos, anchor = '+', 'NA', 'NA'
        selected = 0
        min_identity = 1.0
        status = 'single_unit_not_circularity_proof'
    else:
        anchor_result = candidate_anchor(reference, raw, n, args.anchor_length,
                                         round(length * args.unit_length_tolerance))
        if anchor_result is None:
            raise ValueError('No unique reference anchor recurs once per estimated unit '
                             'at near the expected spacing. Refusing to collapse.')
        _, strand, oriented_raw, refpos, anchor, occ, periods = anchor_result
        units, intervals = [], []
        for i in range(n):
            start = occ[i]
            end = occ[i+1] if i < n-1 else len(oriented_raw) + occ[0]
            unit = (oriented_raw + oriented_raw)[start:end] if end > len(oriented_raw) else oriented_raw[start:end]
            units.append(unit)
            intervals.append({
                'copy': i+1, 'start_0based': start,
                'end_0based_exclusive': end % len(oriented_raw),
                'wraps_raw_origin': 'yes' if end > len(oriented_raw) else 'no',
                'length_bp': len(unit),
            })
        pair_rows = []
        distance_sums = [0] * n
        for i in range(n):
            for j in range(i+1, n):
                row = pairwise_counts(units[i], units[j], aligner)
                distance_sums[i] += row['differences']
                distance_sums[j] += row['differences']
                pair_rows.append({'copy_a': i+1, 'copy_b': j+1, **row})
        min_identity = min(row['identity'] for row in pair_rows)
        if min_identity < args.min_copy_identity:
            raise ValueError(f'Minimum copy identity {min_identity:.5%} is below '
                             f'{args.min_copy_identity:.5%}; the copies may differ. '
                             'Refusing to collapse.')
        # The observed unit with fewest pairwise differences is retained.
        # Ties prefer the unit with length nearest expectation, then the first unit.
        # Sequence identity to the design is never used for copy selection.
        selected = min(range(n), key=lambda i: (distance_sums[i], abs(len(units[i])-length), i))
        status = 'putative_tandem_normalized_requires_read_review'

    name = args.sample if args.sample else raw_id
    write_fasta(args.copies_fasta, [(f'{name}_flye_copy{i+1}', seq) for i,seq in enumerate(units)])
    write_fasta(output, [(f'{name}_flye_normalized', units[selected])])
    write_tsv(args.boundaries, intervals, ['copy', 'start_0based', 'end_0based_exclusive',
                                            'wraps_raw_origin', 'length_bp'])
    write_tsv(args.pairwise, pair_rows, ['copy_a', 'copy_b', 'identities', 'mismatches',
                                         'gap_bases', 'identity', 'differences', 'aligned_columns'])
    write_tsv(args.report, [{
        'sample': name, 'status': status, 'raw_contig_id': raw_id,
        'raw_length_bp': len(raw), 'design_length_bp': length,
        'estimated_copies': n, 'anchor_length_bp': args.anchor_length if n > 1 else 'NA',
        'anchor_reference_start_0based': refpos, 'raw_orientation_for_detection': strand,
        'copy_lengths_bp': ','.join(map(str, [len(s) for s in units])),
        'selected_copy': selected + 1, 'normalized_length_bp': len(units[selected]),
        'min_pairwise_identity': f'{min_identity:.8f}',
        'notes': 'Copy selection is NOT reference-guided; manual review of reads is required '
                 'to distinguish an assembly artifact from a biological concatemer.'
    }], ['sample', 'status', 'raw_contig_id', 'raw_length_bp', 'design_length_bp',
         'estimated_copies', 'anchor_length_bp', 'anchor_reference_start_0based',
         'raw_orientation_for_detection', 'copy_lengths_bp', 'selected_copy',
         'normalized_length_bp', 'min_pairwise_identity', 'notes'])
    print(f'{name}: {status}; raw={len(raw)}bp, units={n}, '
          f'lengths={[len(s) for s in units]}, selected={selected+1}, '
          f'min_identity={min_identity:.5%}', file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--assembly', required=True, help='Original Flye assembly.fasta')
    p.add_argument('--design', required=True, help='Designed plasmid FASTA')
    p.add_argument('--sample', default='')
    p.add_argument('--output-fasta', required=True)
    p.add_argument('--report', required=True)
    p.add_argument('--boundaries', required=True)
    p.add_argument('--pairwise', required=True)
    p.add_argument('--copies-fasta', required=True)
    p.add_argument('--anchor-length', type=int, default=31)
    p.add_argument('--min-copy-identity', type=float, default=0.995)
    p.add_argument('--copy-number-tolerance', type=float, default=0.05,
                   help='Allowed difference between length/design and nearest integer')
    p.add_argument('--unit-length-tolerance', type=float, default=0.05,
                   help='Allowed relative deviation of each inferred unit length')
    p.add_argument('--max-copies', type=int, default=6)
    args = p.parse_args()
    if not 0 < args.min_copy_identity <= 1:
        p.error('--min-copy-identity must lie in (0,1]')
    if args.anchor_length < 15:
        p.error('--anchor-length must be >=15')
    run(args)

if __name__ == '__main__':
    main()
