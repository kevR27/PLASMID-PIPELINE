"""Quick offline tests of the pipeline's scripts (no real reads, NANOPlasty or BLAST needed).

Run from the project folder:   python3 tests/test_scripts.py
"""
import csv
import gzip
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / 'scripts'
sys.path.insert(0, str(SCRIPTS))

from compare_sequences import compare                              # noqa: E402
from differences_to_vcf import write_vcf                           # noqa: E402
from orient_assembly import orient                                 # noqa: E402
from plasmid_utils import (is_same_circular_sequence, read_single_sequence,   # noqa: E402
                           read_table, reverse_complement, write_table)


def run(script, *args):
    subprocess.run([sys.executable, str(SCRIPTS / script), *map(str, args)], check=True)


def write_fasta_file(path, name, sequence):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(f'>{name}\n{sequence}\n')


def fake_blast(path, lines):
    """A tiny stand-in for blastn that prints fixed hits (tests our parsing, not BLAST)."""
    printf_args = ''.join(f'{line}\\n' for line in lines)
    Path(path).write_text(f"#!/bin/sh\nprintf '{printf_args}'\n")
    Path(path).chmod(0o755)


def test_project_files():
    with (ROOT / 'samples.csv').open(newline='') as handle:
        samples = list(csv.DictReader(handle))
    assert len(samples) == 9
    for row in samples:
        design = read_single_sequence(ROOT / row['design_fasta'])
        seed = read_single_sequence(ROOT / row['seed_fasta'])
        assert len(design) == int(row['expected_length_bp'])
        assert len(seed) == 800 and design.count(seed) == 1


def test_orient_and_compare():
    design = read_single_sequence(ROOT / 'designs' / 'TAL MERRF 024128 left clone 9-YES.fa')
    rotated = design[2000:] + design[:2000]
    oriented, strand, _, _, status = orient(design, rotated)
    assert oriented == design and strand == '+' and status == 'OK'
    flipped = reverse_complement(rotated)
    oriented, strand, _, _, status = orient(design, flipped)
    assert oriented == design and strand == '-' and status == 'OK'
    assert is_same_circular_sequence(design, flipped)

    summary, differences = compare('ACGTTAGC', 'ACGTCAGC')
    assert [d['type'] for d in differences] == ['SNV'] and differences[0]['position_1based'] == 5
    assert summary['identity_percent'] == 87.5
    assert compare('ACGTTAGC', 'ACGTTTAGC')[0]['inserted_bp'] == 1
    assert compare('ACGTTAGC', 'ACGTAGC')[0]['deleted_bp'] == 1
    assert compare('ACGTTAGC', 'ACGTTAGC')[0]['identical'] == 'YES'


def test_vcf(tmp):
    design = tmp / 'design.fasta'
    write_fasta_file(design, 'ref', 'ACGTTAGC')
    columns = ['type', 'position_1based', 'design_sequence', 'assembly_sequence', 'length_bp']
    cases = [('ACGTCAGC', '\t5\t.\tT\tC\t'), ('ACGTTTAGC', '\t3\t.\tG\tGT\t'),
             ('ACGTAGC', '\t3\t.\tGT\tG\t'), ('TACGTTAGC', 'CIRCULAR_EDGE')]
    for i, (assembly, expected) in enumerate(cases):
        _, differences = compare('ACGTTAGC', assembly)
        rows = [dict(type=d['type'], position_1based=d['position_1based'],
                     design_sequence=d['reference_sequence'], assembly_sequence=d['other_sequence'],
                     length_bp=d['length_bp']) for d in differences]
        write_table(tmp / f'{i}.tsv', rows, columns)
        write_vcf(design, tmp / f'{i}.tsv', tmp / f'{i}.vcf')
        assert expected in (tmp / f'{i}.vcf').read_text(), (assembly, (tmp / f'{i}.vcf').read_text())
    write_table(tmp / 'none.tsv', [], columns)
    assert write_vcf(design, tmp / 'none.tsv', tmp / 'none.vcf') == 0


def test_filter_reads_and_config(tmp):
    design = read_single_sequence(ROOT / 'designs' / 'TAL MERRF 024128 left clone 9-YES.fa')
    fastq = tmp / 'reads.fastq.gz'
    with gzip.open(fastq, 'wt') as handle:
        for i, sequence in enumerate((design[:300], design[:1300], design[1200:2600])):
            handle.write(f'@read{i}\n{sequence}\n+\n{"I" * len(sequence)}\n')
    run('filter_reads.py', '--fastq-pattern', fastq, '--output-fastq', tmp / 'filtered.fastq',
        '--stats', tmp / 'stats.tsv', '--length-histogram', tmp / 'hist.tsv', '--sample', 'test',
        '--barcode', 'barcode01', '--design-length', 6655, '--min-length', 1000, '--min-quality', 10)
    stats = read_table(tmp / 'stats.tsv')[0]
    assert stats['all_reads'] == '3' and stats['reads_kept'] == '2' and stats['bases_kept'] == '2700'

    readsdir = tmp / 'filtered_reads'
    readsdir.mkdir(exist_ok=True)
    (readsdir / 'reads.fastq').write_bytes((tmp / 'filtered.fastq').read_bytes())
    run('write_nanoplasty_config.py', '--sample', 'test',
        '--seed', ROOT / 'seeds' / 'tal_merrf_common_seed_800bp_clean.fasta',
        '--design-length', 6655, '--reads-folder', readsdir,
        '--output-folder', tmp / 'nanoplasty' / 'native',
        '--config-file', tmp / 'config.txt', '--threads', 4,
        '--min-read-length', 3000, '--range-fraction', 0.15,
        '--use-full-dataset', 'no', '--extend-seed-directly', 'no',
        '--save-assembled-reads', 'no')
    text = (tmp / 'config.txt').read_text()
    assert 'Type                  = mito' in text
    assert 'R10                   = yes' in text
    assert 'Extend seed directly  = no' in text
    assert 'Genome Range          = ' in text
    assert not (tmp / 'nanoplasty' / 'native').exists()



def test_orient_script_and_plots(tmp):
    design = read_single_sequence(ROOT / 'designs' / 'TAL MERRF 024128 left clone 9-YES.fa')
    flipped = reverse_complement(design[2000:] + design[:2000])
    write_fasta_file(tmp / 'assembly.fasta', 'contig', flipped)
    write_fasta_file(tmp / 'design.fasta', 'ref', design)
    run('orient_assembly.py', '--assembly', tmp / 'assembly.fasta', '--design', tmp / 'design.fasta', '--name', 'test',
        '--output-fasta', tmp / 'oriented.fasta', '--output-report', tmp / 'orientation.tsv')
    assert read_single_sequence(tmp / 'oriented.fasta') == design
    report = read_table(tmp / 'orientation.tsv')[0]
    assert report['status'] == 'OK' and report['needs_review'] == 'NO' and report['identical_circular_sequence'] == 'YES'
    run('compare_sequences.py', '--reference', tmp / 'design.fasta', '--other', tmp / 'oriented.fasta',
        '--reference-label', 'design', '--other-label', 'assembly',
        '--summary', tmp / 'summary.tsv', '--differences', tmp / 'differences.tsv')
    summary = read_table(tmp / 'summary.tsv')[0]
    assert summary['identical'] == 'YES' and 'design_bp' in summary and 'assembly_bp' in summary

    (tmp / 'depth.tsv').write_text(''.join(f'ref\t{i}\t{20 + i % 3}\n' for i in range(1, len(design) + 1)))
    # histogram/stats come from the filter_reads test above (same temp folder)
    run('plot_sample.py', '--sample', 'test', '--read-stats', tmp / 'stats.tsv', '--length-histogram', tmp / 'hist.tsv',
        '--depth', tmp / 'depth.tsv', '--differences', tmp / 'differences.tsv', '--output-folder', tmp / 'plots')
    for plot in ('read_lengths', 'depth', 'differences'):
        for extension in ('png', 'pdf'):
            assert (tmp / 'plots' / f'{plot}.{extension}').stat().st_size > 1000


def test_inserts_blast_and_sniffles(tmp):
    sequence = 'ACGTTAGC' * 10                      # 80 bp "plasmid"
    design, assembly = tmp / 'design.fasta', tmp / 'assembly.fasta'
    write_fasta_file(design, 'ref', sequence)
    write_fasta_file(assembly, 'asm', sequence)
    coordinates = tmp / 'coordinates.csv'

    def check(*extra, blast=None):
        run('check_inserts.py', '--design', design, '--assembly', assembly, '--coordinates', coordinates,
            '--construct', 'tal_merrf_left', '--output', tmp / 'inserts.tsv', '--queries-out', tmp / 'queries.fa',
            *extra, *(['--blastn', blast] if blast else []))
        return {r['insert_id']: r for r in read_table(tmp / 'inserts.tsv')}

    coordinates.write_text('construct,insert_id,start_1based,end_1based\n')
    assert check()['NOT_ANNOTATED']['status'] == 'NO_INSERTS_PROVIDED'

    coordinates.write_text('construct,insert_id,start_1based,end_1based\ntal_merrf_left,insertA,1,40\n')
    blast_hits = fake_blast(tmp / 'blast1', ['insertA\\tassembly_circular_doubled\\t100.000\\t40\\t0\\t0\\t1\\t40\\t1\\t40\\t1e-15\\t50\\t40\\t160'])
    assert check(blast=tmp / 'blast1')['insertA']['status'] == 'FULL_LENGTH_EXACT_SEQUENCE'

    inserts_fasta = tmp / 'inserts.fasta'
    inserts_fasta.write_text(f'>tal_merrf_left|fa_insert\n{sequence[:30]}\n>shared\n{sequence[8:38]}\n'
                             f'>other_construct|skipme\n{sequence[:30]}\n')
    fake_blast(tmp / 'blast2', [
        'fa_insert\\tassembly_circular_doubled\\t100.000\\t30\\t0\\t0\\t1\\t30\\t45\\t16\\t1e-15\\t50\\t30\\t160',
        'fa_insert\\tassembly_circular_doubled\\t100.000\\t30\\t0\\t0\\t1\\t30\\t125\\t96\\t1e-15\\t50\\t30\\t160'])
    results = check('--insert-sequences', inserts_fasta, blast=tmp / 'blast2')
    assert set(results) == {'insertA', 'fa_insert', 'shared'}               # 'skipme' is another construct
    assert results['fa_insert']['status'] == 'FULL_LENGTH_EXACT_SEQUENCE' and results['fa_insert']['source'] == 'fasta'
    assert results['fa_insert']['exact_copies'] == '1'                       # hits 80 bp apart = one circular copy
    assert results['shared']['status'] == 'NO_HIT'

    # BLAST summary: best hit and the union of overlapping hits
    (tmp / 'blast.tsv').write_text('asm\tref\t100\t40\t0\t0\t1\t40\t1\t40\t1e-15\t50\t80\t80\n'
                                   'asm\tref\t100\t40\t0\t0\t30\t70\t30\t70\t1e-15\t50\t80\t80\n')
    run('summarize_blast.py', '--blast', tmp / 'blast.tsv', '--output', tmp / 'blast_summary.tsv')
    summary = read_table(tmp / 'blast_summary.tsv')[0]
    assert float(summary['best_hsp_query_coverage_percent']) == 51.25 and float(summary['query_union_coverage_percent']) == 87.5

    # Sniffles summary: calls near the artificial ends are counted separately
    header = '##fileformat=VCFv4.2\n##contig=<ID=p,length=1000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    (tmp / 'design.vcf').write_text(header + 'p\t10\t.\tN\t<DEL>\t.\tPASS\tEND=300;SVTYPE=DEL\n'
                                    'p\t500\t.\tN\t<INS>\t.\tPASS\tEND=500;SVTYPE=INS\n'
                                    'p\t600\t.\tN\t<DEL>\t.\tGT\tEND=980;SVTYPE=DEL\n')
    (tmp / 'rotated.vcf').write_text(header + 'p\t500\t.\tN\t<INS>\t.\tPASS\tEND=500;SVTYPE=INS\n')
    run('summarize_sniffles.py', '--design-vcf', tmp / 'design.vcf', '--rotated-vcf', tmp / 'rotated.vcf',
        '--output', tmp / 'sniffles_summary.tsv')
    row = read_table(tmp / 'sniffles_summary.tsv')[0]
    assert (row['calls_on_design'], row['pass_on_design'], row['near_ends_of_design']) == ('3', '2', '2'), row
    assert (row['calls_on_rotated_design'], row['near_ends_of_rotated_design']) == ('1', '0'), row


def test_summary_table(tmp):
    """Fabricate the per-sample result files and check the summary + review flag logic."""
    sheet = tmp / 'samples.csv'
    sheet.write_text('barcode,sample,construct,expected_length_bp,design_fasta,seed_fasta\n'
                     'barcode01,s1,c1,8,x,y\nbarcode02,s2,c1,8,x,y\n')
    results = tmp / 'results'
    for sample in ('s1', 's2'):
        folder = results / sample
        write_table(folder / 'reads' / 'read_stats.tsv', [dict(estimated_depth=50)], ['estimated_depth'])
        write_table(folder / 'nanoplasty' / 'assembly_source.tsv',
                    [dict(source_fasta='native/Circularized_assembly.fasta')], ['source_fasta'])
        write_fasta_file(folder / 'assembly' / 'assembly_oriented.fasta', 'a', 'ACGTTAGC')
        write_table(folder / 'assembly' / 'orientation.tsv',
                    [dict(number_of_contigs=1, status='OK', identical_circular_sequence='YES', needs_review='NO')],
                    ['number_of_contigs', 'status', 'identical_circular_sequence', 'needs_review'])
        write_table(folder / 'assembly_vs_design' / 'summary.tsv',
                    [dict(assembly_bp=8, identity_percent=100, substituted_bp=0, inserted_bp=0, deleted_bp=0,
                          differences=0, identical='YES')],
                    ['assembly_bp', 'identity_percent', 'substituted_bp', 'inserted_bp', 'deleted_bp', 'differences', 'identical'])
        (folder / 'assembly_vs_design' / 'differences.vcf').write_text('#CHROM\tPOS\n')
        write_table(folder / 'assembly_vs_design' / 'reads_across_origin.tsv', [dict(alignments_crossing_origin=5)],
                    ['alignments_crossing_origin'])
        (folder / 'alignments').mkdir(parents=True)
        (folder / 'alignments' / 'read_coverage.tsv').write_text('#rname\tcoverage\tmeandepth\nref\t100\t42\n')
        write_table(folder / 'blast' / 'whole_plasmid_summary.tsv',
                    [dict(best_hsp_identity_percent=100, best_hsp_query_coverage_percent=100, query_union_coverage_percent=100)],
                    ['best_hsp_identity_percent', 'best_hsp_query_coverage_percent', 'query_union_coverage_percent'])
        status = 'FULL_LENGTH_EXACT_SEQUENCE' if sample == 's1' else 'NO_HIT'
        write_table(folder / 'blast' / 'inserts_summary.tsv', [dict(insert_id='ins', status=status)], ['insert_id', 'status'])
        write_table(folder / 'sniffles' / 'summary.tsv', [dict(calls_on_design=0, calls_on_rotated_design=0)],
                    ['calls_on_design', 'calls_on_rotated_design'])
    run('make_summary_table.py', '--sample-sheet', sheet, '--results-folder', results,
        '--summary-out', tmp / 'all.tsv', '--replicates-out', tmp / 'rep.tsv')
    rows = {r['sample']: r for r in read_table(tmp / 'all.tsv')}
    assert rows['s1']['needs_review'] == 'NO' and rows['s2']['needs_review'] == 'YES'   # s2 insert missing
    assert rows['s1']['vcf_records'] == '0' and rows['s1']['mapped_mean_depth'] == '42'
    assert read_table(tmp / 'rep.tsv')[0]['identical_circular_sequence'] == 'YES'


def main():
    test_project_files()
    test_orient_and_compare()
    with tempfile.TemporaryDirectory() as folder:
        tmp = Path(folder)
        test_vcf(tmp)
        test_filter_reads_and_config(tmp)
        test_orient_script_and_plots(tmp)
        test_inserts_blast_and_sniffles(tmp)
        test_summary_table(tmp)
    print('PASS: project files, orientation, comparison, VCF, read filtering, NANOPlasty config, '
          'plots, insert check, BLAST and Sniffles summaries, summary table')


if __name__ == '__main__':
    main()
