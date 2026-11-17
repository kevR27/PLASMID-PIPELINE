#!/usr/bin/env python3
"""Check that your inserted sequences are present in the assembled plasmid (BLAST, run locally).

Two ways to say what the inserts are (they can be combined):
  1. --coordinates insert_coordinates.csv : positions (1-based, inclusive) in the DESIGNED
                                  plasmid. The insert is cut out of the design FASTA.
  2. --insert-sequences inserts.fasta : the insert sequences themselves. FASTA header is
                                  ">construct|insert_id" (only for that construct) or
                                  ">insert_id" (checked in every construct).

Each insert is BLASTed against the assembled plasmid written twice head-to-tail, so a hit
may cross the artificial FASTA origin. For every insert we report whether a full-length
100% match exists, how many distinct full-length copies there are, and where they sit
in the (design-oriented) assembly.

A full-length exact hit proves the SEQUENCE is present. It does not by itself prove
the junctions to the vector are correct; that is what the whole-assembly comparison
against the design is for.
"""
import argparse
import csv
import subprocess
from pathlib import Path

from Bio import SeqIO

from plasmid_utils import read_single_sequence, write_table

BLAST_FIELDS = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']
OUT_FIELDS = ['construct', 'insert_id', 'source', 'start_1based', 'end_1based', 'query_bp',
              'best_identity_percent', 'best_query_coverage_percent',
              'exact_copies', 'assembly_positions_1based', 'status']
MIN_INSERT_BP = 11  # BLASTN word size


def check_id(name):
    if not name or not all(c.isalnum() or c in '_-' for c in name):
        raise ValueError(f'Insert ID may contain only letters, digits, underscore or hyphen: {name!r}')


def inserts_from_coordinates(design_seq, csv_path, construct):
    """Cut the insert fragments out of the designed plasmid."""
    inserts = []
    with open(csv_path, newline='') as handle:
        for row in csv.DictReader(handle):
            if row['construct'].strip() != construct or not row.get('insert_id', '').strip():
                continue
            name = row['insert_id'].strip()
            check_id(name)
            start, end = int(row['start_1based']), int(row['end_1based'])
            if not (1 <= start <= len(design_seq) and 1 <= end <= len(design_seq)):
                raise ValueError(f'Invalid insert interval: {row}')
            # start > end means the insert crosses the artificial origin of the FASTA
            seq = design_seq[start - 1:end] if end >= start else design_seq[start - 1:] + design_seq[:end]
            inserts.append(dict(insert_id=name, source='design_coordinates',
                                start=start, end=end, seq=seq))
    return inserts


def inserts_from_fasta(fasta_path, construct):
    """Read insert sequences supplied directly by the user."""
    inserts = []
    for record in SeqIO.parse(fasta_path, 'fasta'):
        owner, _, name = record.id.partition('|')
        if not name:                      # no construct prefix: applies to every construct
            owner, name = '', record.id
        if owner and owner != construct:
            continue
        check_id(name)
        inserts.append(dict(insert_id=name, source='fasta', start='', end='',
                            seq=str(record.seq).upper()))
    return inserts


def run_blast(blastn, query_fasta, subject_fasta):
    result = subprocess.run(
        [blastn, '-task', 'blastn', '-dust', 'no', '-query', str(query_fasta),
         '-subject', str(subject_fasta), '-outfmt', '6 ' + ' '.join(BLAST_FIELDS)],
        text=True, capture_output=True, check=True)
    hits = {}
    for line in result.stdout.splitlines():
        parts = line.split('\t')
        if len(parts) == len(BLAST_FIELDS):
            row = dict(zip(BLAST_FIELDS, parts))
            hits.setdefault(row['qseqid'], []).append(row)
    return hits


def evaluate(insert, hits, assembly_len):
    """Best hit, and the distinct full-length exact copies, for one insert."""
    n = len(insert['seq'])
    best = (False, 0.0, 0.0)           # (exact, query coverage %, identity %)
    copies = set()
    for hit in hits:
        qs, qe = int(hit['qstart']), int(hit['qend'])
        coverage = 100 * (abs(qe - qs) + 1) / n
        identity = float(hit['pident'])
        exact = (qs == 1 and qe == n and int(hit['mismatch']) == 0 and int(hit['gapopen']) == 0
                 and int(hit['length']) == n and abs(identity - 100.0) < 1e-6)
        if exact:
            s0, s1 = int(hit['sstart']), int(hit['send'])
            # modulo the assembly length: hits in the two halves of the doubled subject are one copy
            copies.add(((min(s0, s1) - 1) % assembly_len, '+' if s1 >= s0 else '-'))
        best = max(best, (exact, coverage, identity))
    positions = ';'.join(f'{p + 1}-{(p + n - 1) % assembly_len + 1}({strand})'
                         for p, strand in sorted(copies))
    status = 'FULL_LENGTH_EXACT_SEQUENCE' if best[0] else ('PARTIAL_OR_DIFFERENT' if hits else 'NO_HIT')
    return best, len(copies), positions, status


def run(args):
    design = read_single_sequence(args.design)
    assembled = read_single_sequence(args.assembly)
    construct = args.construct
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    query_path = Path(args.queries_out)
    query_path.parent.mkdir(parents=True, exist_ok=True)

    inserts = inserts_from_coordinates(design, args.coordinates, construct)
    if args.insert_sequences:
        inserts += inserts_from_fasta(args.insert_sequences, construct)
    names = [i['insert_id'] for i in inserts]
    if len(names) != len(set(names)):
        raise ValueError(f'Duplicate insert IDs for {construct}: {names}')
    for item in inserts:
        if len(item['seq']) < MIN_INSERT_BP:
            raise ValueError(f"Insert {item['insert_id']} is shorter than {MIN_INSERT_BP} bp")

    if not inserts:   # still write a table so the report says clearly that nothing was assessed
        query_path.write_text('')
        write_table(out, [dict(construct=construct, insert_id='NOT_ANNOTATED',
                               status='NO_INSERTS_PROVIDED')], OUT_FIELDS)
        return

    query_path.write_text(''.join(f">{i['insert_id']}\n{i['seq']}\n" for i in inserts))
    doubled = out.with_suffix('.circular_subject.fasta')
    doubled.write_text(f'>assembly_circular_doubled\n{assembled * 2}\n')
    try:
        hits = run_blast(args.blastn, query_path, doubled)
    finally:
        doubled.unlink(missing_ok=True)

    rows = []
    for item in inserts:
        (exact, coverage, identity), copies, positions, status = evaluate(
            item, hits.get(item['insert_id'], []), len(assembled))
        rows.append(dict(construct=construct, insert_id=item['insert_id'], source=item['source'],
                         start_1based=item['start'], end_1based=item['end'],
                         query_bp=len(item['seq']), best_identity_percent=round(identity, 4),
                         best_query_coverage_percent=round(coverage, 4), exact_copies=copies,
                         assembly_positions_1based=positions, status=status))
    write_table(out, rows, OUT_FIELDS)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    for name in ('design', 'assembly', 'coordinates', 'construct', 'output', 'queries-out'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--insert-sequences', default='',
                        help='optional FASTA with the insert sequences themselves')
    parser.add_argument('--blastn', default='blastn')
    run(parser.parse_args())


if __name__ == '__main__':
    main()
