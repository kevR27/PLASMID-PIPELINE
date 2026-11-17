#!/usr/bin/env python3
"""Summarize the whole-plasmid BLAST table (assembly = query, design = subject).

BLAST reports local HSPs, never a global identity, and a plasmid is often split into
several HSPs. Besides the single best HSP we therefore report how much of the query
and of the subject is covered by the UNION of all HSPs. Use the global alignment
(compare_consensus.py) as the primary identity metric.
"""
import argparse

from plasmid_utils import write_table

COLUMNS = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
           'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']
OUT_FIELDS = ['best_hsp_identity_percent', 'best_hsp_query_coverage_percent', 'best_hsp_length',
              'hit_count', 'query_union_coverage_percent', 'subject_union_coverage_percent']


def union_length(intervals):
    """Total bases covered by 1-based inclusive intervals (overlaps counted once)."""
    total, current_end = 0, 0
    for start, end in sorted((min(a, b), max(a, b)) for a, b in intervals):
        if end <= current_end:
            continue
        total += end - max(start, current_end + 1) + 1
        current_end = end
    return total


def summarize(hits):
    if not hits:
        return dict(best_hsp_identity_percent='NA', best_hsp_query_coverage_percent='0',
                    best_hsp_length=0, hit_count=0,
                    query_union_coverage_percent=0, subject_union_coverage_percent=0)
    span = lambda r: abs(int(r['qend']) - int(r['qstart'])) + 1
    best = max(hits, key=lambda r: (span(r), float(r['pident'])))
    qlen, slen = int(best['qlen']), int(best['slen'])
    q_union = union_length([(int(r['qstart']), int(r['qend'])) for r in hits])
    s_union = union_length([(int(r['sstart']), int(r['send'])) for r in hits])
    return dict(best_hsp_identity_percent=best['pident'],
                best_hsp_query_coverage_percent=round(100 * span(best) / qlen, 4),
                best_hsp_length=best['length'], hit_count=len(hits),
                query_union_coverage_percent=round(100 * q_union / qlen, 4),
                subject_union_coverage_percent=round(100 * s_union / slen, 4))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--blast', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    hits = []
    with open(args.blast) as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) == len(COLUMNS):
                hits.append(dict(zip(COLUMNS, fields)))
    write_table(args.output, [summarize(hits)], OUT_FIELDS)


if __name__ == '__main__':
    main()
