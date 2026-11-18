#!/usr/bin/env python3
"""Collect the results of all samples into one summary table, plus a barcode-vs-barcode
comparison for replicates of the same construct.

Reads the per-sample files in results/<sample>/ (see the README for the folder layout).
The column `needs_review` is YES when any of these is true:
  - the assembly has several contigs or could not be oriented to the design
  - the assembly is not identical to the design
  - Flye (if run) disagrees with NANOPlasty
  - an insert you asked about is not found as a full-length exact match
"""
import argparse
import csv
from itertools import combinations
from pathlib import Path

from compare_sequences import compare
from plasmid_utils import is_same_circular_sequence, read_single_sequence, read_table, write_table


def first_row(path):
    return read_table(path)[0]


def read_coverage_table(path):
    """`samtools coverage` output: header starts with '#', first data row is our plasmid."""
    with open(path) as handle:
        lines = [line.rstrip('\n').split('\t') for line in handle if line.strip()]
    if len(lines) < 2:
        return {}
    return dict(zip([name.lstrip('#') for name in lines[0]], lines[1]))


def count_vcf_records(path):
    with open(path) as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith('#'))


def summarize_sample(sample_row, results, use_sniffles, use_flye):
    sample = sample_row['sample']
    folder = results / sample
    reads = first_row(folder / 'reads' / 'read_stats.tsv')
    nanoplasty_source = first_row(folder / 'nanoplasty' / 'assembly_source.tsv')
    orientation = first_row(folder / 'assembly' / 'orientation.tsv')
    versus_design = first_row(folder / 'assembly_vs_design' / 'summary.tsv')
    origin = first_row(folder / 'assembly_vs_design' / 'reads_across_origin.tsv')
    coverage = read_coverage_table(folder / 'alignments' / 'read_coverage.tsv')
    blast = first_row(folder / 'blast' / 'whole_plasmid_summary.tsv')

    insert_rows = [r for r in read_table(folder / 'blast' / 'inserts_summary.tsv')
                   if r.get('insert_id') != 'NOT_ANNOTATED']
    inserts_found = sum(1 for r in insert_rows if r['status'] == 'FULL_LENGTH_EXACT_SEQUENCE')
    repeat_rows = [r for r in read_table(folder / 'repeat_integrity' / 'summary.tsv')
                   if r.get('candidate') != 'NONE']
    repeat_matches = sum(r['status'] == 'SEQUENCE_MATCH' for r in repeat_rows)
    repeat_review = not repeat_rows or repeat_matches != len(repeat_rows)
    spanning = [int(r['reference_spanning_reads']) for r in repeat_rows
                if r.get('reference_spanning_reads', '').isdigit()]
    support_assessed = bool(repeat_rows) and len(spanning) == len(repeat_rows)
    min_spanning = min(spanning) if support_assessed else 'NOT_ASSESSED'
    limited_repeat_support = not support_assessed or min_spanning < 3

    row = dict(
        sample=sample, barcode=sample_row['barcode'], construct=sample_row['construct'],
        assembler='NANOPlasty 0.4', nanoplasty_source=nanoplasty_source['source_fasta'],
        expected_length_bp=sample_row['expected_length_bp'],
        assembly_length_bp=versus_design['assembly_bp'],
        contigs=orientation['number_of_contigs'], orientation_status=orientation['status'],
        identical_circular_sequence=orientation['identical_circular_sequence'],
        reads_crossing_origin=origin['alignments_crossing_origin'],
        estimated_depth=reads['estimated_depth'],
        mapped_coverage_percent=coverage.get('coverage', 'NA'),
        mapped_mean_depth=coverage.get('meandepth', 'NA'),
        assembly_identity_percent=versus_design['identity_percent'],
        substituted_bp=versus_design['substituted_bp'], inserted_bp=versus_design['inserted_bp'],
        deleted_bp=versus_design['deleted_bp'], number_of_differences=versus_design['differences'],
        vcf_records=count_vcf_records(folder / 'assembly_vs_design' / 'differences.vcf'),
        blast_best_hit_identity_percent=blast['best_hsp_identity_percent'],
        blast_best_hit_query_coverage_percent=blast['best_hsp_query_coverage_percent'],
        blast_all_hits_query_coverage_percent=blast['query_union_coverage_percent'],
        inserts_checked=len(insert_rows) if insert_rows else 'NOT_PROVIDED',
        inserts_found_full_length=inserts_found if insert_rows else 'NOT_CHECKED',
        insert_validation_status=(
            ('DESIGN_MATCH_INSERT_NOT_INDEPENDENTLY_LOCALIZED'
             if versus_design['identical'] == 'YES' else
             'INSERT_NOT_LOCALIZED_DESIGN_DIFFERS')
            if not insert_rows else
            'ALL_EXACT' if inserts_found == len(insert_rows) else 'REVIEW_REQUIRED'),
        repeat_candidate_regions=len(repeat_rows),
        repeat_candidate_regions_matching=repeat_matches,
        repeat_min_spanning_alignments=min_spanning,
        repeat_integrity_status=('NO_CANDIDATES' if not repeat_rows else
                                 'REVIEW_REQUIRED' if repeat_review else
                                 'DESIGN_MATCH_LIMITED_READ_SUPPORT' if limited_repeat_support
                                 else 'DESIGN_MATCH_WITH_SPANNING_ALIGNMENTS'))

    if use_sniffles:
        sniffles = first_row(folder / 'sniffles' / 'summary.tsv')
        row.update(sniffles_calls_on_design=sniffles['calls_on_design'],
                   sniffles_calls_on_rotated_design=sniffles['calls_on_rotated_design'])
    else:
        row.update(sniffles_calls_on_design='NOT_RUN', sniffles_calls_on_rotated_design='NOT_RUN')

    flye_agrees = True
    if use_flye:
        flye_orientation = first_row(folder / 'flye' / 'orientation.tsv')
        flye_vs_design = first_row(folder / 'flye' / 'vs_design.tsv')
        flye_vs_nanoplasty = first_row(folder / 'flye' / 'vs_nanoplasty.tsv')
        flye_agrees = flye_vs_nanoplasty['identical'] == 'YES'
        row.update(flye_contigs=flye_orientation['number_of_contigs'],
                   flye_identity_to_design_percent=flye_vs_design['identity_percent'],
                   flye_identity_to_nanoplasty_percent=flye_vs_nanoplasty['identity_percent'],
                   flye_agrees_with_nanoplasty='YES' if flye_agrees else 'NO')

    problem = (orientation['needs_review'] == 'YES' or versus_design['identical'] != 'YES'
               or not flye_agrees or (insert_rows and inserts_found < len(insert_rows))
               or repeat_review or limited_repeat_support)
    row['whole_construct_match'] = ('YES' if versus_design['identical'] == 'YES' else 'NO')
    row['retention_screen_result'] = ('CONSISTENT_WITH_DESIGN_AND_REPEATS'
                                      if not problem else 'REVIEW_REQUIRED')
    row['needs_review'] = 'YES' if problem else 'NO'
    return row


def compare_replicates(samples, results):
    """Compare every pair of barcodes that carry the same construct."""
    rows = []
    for a, b in combinations(samples, 2):
        if a['construct'] != b['construct']:
            continue
        seq_a = read_single_sequence(results / a['sample'] / 'assembly' / 'assembly_oriented.fasta')
        seq_b = read_single_sequence(results / b['sample'] / 'assembly' / 'assembly_oriented.fasta')
        summary, _ = compare(seq_a, seq_b)
        rows.append(dict(construct=a['construct'], barcode_a=a['barcode'], barcode_b=b['barcode'],
                         identical_circular_sequence='YES' if is_same_circular_sequence(seq_a, seq_b) else 'NO',
                         identity_percent=summary['identity_percent'],
                         number_of_differences=summary['differences']))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    for name in ('sample-sheet', 'results-folder', 'summary-out', 'replicates-out'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--sniffles', choices=['yes', 'no'], default='yes')
    parser.add_argument('--flye', choices=['yes', 'no'], default='no')
    args = parser.parse_args()

    results = Path(args.results_folder)
    with open(args.sample_sheet, newline='') as handle:
        samples = list(csv.DictReader(handle))

    rows = [summarize_sample(s, results, args.sniffles == 'yes', args.flye == 'yes') for s in samples]
    write_table(args.summary_out, rows, list(rows[0]))
    replicates = compare_replicates(samples, results)
    write_table(args.replicates_out, replicates,
                ['construct', 'barcode_a', 'barcode_b', 'identical_circular_sequence',
                 'identity_percent', 'number_of_differences'])


if __name__ == '__main__':
    main()
