#!/usr/bin/env python3
"""Screen design-derived repeat-rich intervals against an oriented plasmid assembly.

Candidate intervals are based on recurrent 30-mers. They may overlap TALE/TALEN
arrays, but do not constitute curated TALEN annotations. Unique exact flanking
anchors are required; spanning-read alignments assess structure, not every base.
"""
import argparse
from collections import Counter
from pathlib import Path

from Bio.Align import PairwiseAligner
from plasmid_utils import read_single_sequence, write_table

FIELDS = ['sample', 'candidate', 'design_start_1based', 'design_end_1based',
          'design_bp', 'assembly_bp', 'alignment_identity_percent',
          'substitutions', 'gap_bases', 'design_marker_30mer_copies',
          'assembly_marker_30mer_copies', 'reference_spanning_reads',
          'status', 'note']


def candidate_regions(sequence, k=30, min_occ=3, join_gap=100, min_span=200):
    counts = Counter(sequence[i:i + k] for i in range(len(sequence) - k + 1))
    starts = [i for i in range(len(sequence) - k + 1)
              if 'N' not in sequence[i:i + k]
              and counts[sequence[i:i + k]] >= min_occ]
    if not starts:
        return []
    blocks = []
    first = prev = starts[0]
    for i in starts[1:]:
        if i - prev > join_gap:
            if prev + k - first >= min_span:
                blocks.append((first, prev + k))
            first = i
        prev = i
    if prev + k - first >= min_span:
        blocks.append((first, prev + k))
    return blocks


def anchored_interval(design, assembly, start, end, anchor_bp=35):
    """Extract an assembled interval between unique design flanks, if possible."""
    if start < anchor_bp or end + anchor_bp > len(design):
        return None
    left = design[start - anchor_bp:start]
    right = design[end:end + anchor_bp]
    if (set(left + right) - set('ACGT')
            or design.count(left) != 1 or design.count(right) != 1
            or assembly.count(left) != 1 or assembly.count(right) != 1):
        return None
    lpos = assembly.find(left) + anchor_bp
    rpos = assembly.find(right)
    observed = assembly[lpos:rpos] if lpos <= rpos else assembly[lpos:] + assembly[:rpos]
    if len(observed) > len(design) * 0.9:
        return None
    return observed


def spanning_reads(bam, start, end, flank):
    if bam is None:
        return 'NOT_ASSESSED'
    import pysam
    count = 0
    with pysam.AlignmentFile(bam, 'rb') as handle:
        if len(handle.references) != 1:
            raise ValueError('Expected one design reference in BAM')
        ref = handle.references[0]
        for read in handle.fetch(ref, start, end):
            if (read.is_unmapped or read.is_secondary or read.is_supplementary
                    or read.is_duplicate or read.mapping_quality < 20):
                continue
            if read.reference_start <= start - flank and read.reference_end >= end + flank:
                count += 1
    return count


def audit(design, assembly, sample='', bam=None, flank=75):
    aligner = PairwiseAligner(mode='global', match_score=2,
                              mismatch_score=-3, open_gap_score=-5,
                              extend_gap_score=-1)
    rows = []
    for idx, (start, end) in enumerate(candidate_regions(design), 1):
        target = design[start:end]
        counts = Counter(target[i:i + 30] for i in range(len(target) - 29))
        marker = max(counts, key=lambda x: (counts[x], x))
        row = {
            'sample': sample, 'candidate': f'candidate_repeat_{idx:02d}',
            'design_start_1based': start + 1, 'design_end_1based': end,
            'design_bp': end - start, 'assembly_bp': 'NOT_ASSESSED',
            'alignment_identity_percent': 'NA', 'substitutions': 'NA',
            'gap_bases': 'NA', 'design_marker_30mer_copies': counts[marker],
            'assembly_marker_30mer_copies': 'NA',
            'reference_spanning_reads': spanning_reads(bam, start, end, flank),
            'status': 'UNRESOLVED', 'note': ''
        }
        observed = anchored_interval(design, assembly, start, end)
        if observed is None:
            row['status'] = 'FLANK_ANCHORS_UNRESOLVED'
            row['note'] = 'Unique exact flanks not found; investigate rearrangement, mismatches or orientation'
        else:
            comparison = aligner.align(target, observed)[0].counts()
            n = comparison.identities + comparison.mismatches + comparison.gaps
            row.update(
                assembly_bp=len(observed),
                alignment_identity_percent=round(100 * comparison.identities / n, 5) if n else 0,
                substitutions=comparison.mismatches,
                gap_bases=comparison.gaps,
                assembly_marker_30mer_copies=sum(
                    observed[i:i + 30] == marker for i in range(len(observed) - 29)),
                status='SEQUENCE_MATCH' if target == observed else 'SEQUENCE_DIFFERENT',
                note=('Inferred repeat-rich interval; not a curated TALEN feature. '
                      'Spanning alignments do not prove internal repeat-module order.')
            )
        rows.append(row)
    if not rows:
        rows.append(dict(sample=sample, candidate='NONE',
                         status='NO_REPEAT_RICH_CANDIDATES',
                         note='No interval passed the repeat-discovery thresholds'))
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--design', required=True)
    p.add_argument('--assembly', required=True)
    p.add_argument('--sample', default='')
    p.add_argument('--bam', default='')
    p.add_argument('--output', required=True)
    a = p.parse_args()
    design = read_single_sequence(a.design)
    assembly = read_single_sequence(a.assembly)
    rows = audit(design, assembly, a.sample, a.bam or None)
    Path(a.output).parent.mkdir(parents=True, exist_ok=True)
    write_table(a.output, rows, FIELDS)
    print(f'{a.sample}: {len(rows)} repeat candidate row(s) written to {a.output}')


if __name__ == '__main__':
    main()
