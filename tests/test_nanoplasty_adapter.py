"""Offline adapter tests: no NANOPlasty installation required."""
import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from run_nanoplasty import capture, find_candidates

class NanoplastyAdapterTests(unittest.TestCase):
    def test_selects_only_full_length_circular_assembly(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); native=root/'nanoplasty'/'native'; native.mkdir(parents=True)
            (native/'Assembled_reads_NP_test.fasta').write_text('>reads\n'+('A'*6000)+'\n')
            (native/'contigs_tmp_test.fasta').write_text('>seed\n'+('C'*800)+'\n')
            (native/'Circularized_assembly_1_test.fasta').write_text('>p\n'+('ACGT'*1500)+'\n')
            assert len(find_candidates(native, 6164)) == 1
            capture(native,6164,'test', root/'output.fa',root/'provenance.tsv')
            with (root/'output.fa').open() as handle:
                assert len(next(SeqIO.parse(handle, 'fasta')).seq) == 6000
            with (root/'provenance.tsv').open() as f:
                assert 'Circularized_assembly_1_test.fasta' in next(csv.DictReader(f,delimiter='\t'))['source_fasta']

    def test_rejects_missing_or_ambiguous_contig(self):
        with tempfile.TemporaryDirectory() as td:
            native=Path(td)/'nanoplasty'/'native'; native.mkdir(parents=True)
            for f in ('Assembled_reads_NP_test.fasta','contigs_tmp_test.fasta'):
                (native/f).write_text('>q\n' + 'A'*6164+'\n')
            with self.assertRaises(RuntimeError):
                capture(native,6164,'test',Path(td)/'out.fa',Path(td)/'s.tsv')
            for f in ('Circularized_assembly_a.fasta','Circularized_assembly_b.fasta'):
                (native/f).write_text('>p\n' + 'A'*6164+'\n')
            with self.assertRaises(RuntimeError):
                capture(native,6164,'test',Path(td)/'out.fa',Path(td)/'s.tsv')

    def test_perl_wrapper_and_isolated_native_output(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            script=root/'fake_nanoplasty.pl'
            # fake Perl assembler exercises the exact launch and output adapter
            script.write_text('''use strict; use warnings;
my $config_file = $ARGV[1];
open my $fh, '<', $config_file or die $!;
my $out; while (<$fh>) { $out = $1 if /^Output path\\s*=\\s*(\\S+)/; }
die "no out" unless $out;
mkdir $out unless -d $out;
open my $fa, '>', "$out/Circularized_assembly_test.fasta" or die $!;
print $fa ">assembled\\n", "ACGT" x 1500, "\\n";
close $fa;
''')
            native=root/'nanoplasty'/'native'; native.mkdir(parents=True)
            (native/'stale.fasta').write_text('>s\nAAA\n')
            config=root/'test.cfg';config.write_text(f'Output path           = {native}/\n')
            cmd=[sys.executable, str(ROOT/'scripts/run_nanoplasty.py'), '--script',str(script),
                 '--config-file',str(config),'--native-folder',str(native),
                 '--sample','test','--design-length','6164',
                 '--output-fasta',str(root/'result.fa'), '--source-report',str(root/'source.tsv')]
            p=subprocess.run(cmd, capture_output=True,text=True)
            self.assertEqual(p.returncode,0,p.stderr)
            self.assertFalse((native/'stale.fasta').exists())
            with (root/'result.fa').open() as handle:
                self.assertEqual(len(next(SeqIO.parse(handle,'fasta')).seq),6000)

if __name__=='__main__':unittest.main()
