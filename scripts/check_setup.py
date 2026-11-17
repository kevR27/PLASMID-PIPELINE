#!/usr/bin/env python3
"""Check your setup before running: sample sheet, designs, seeds, FASTQ files, insert
files and installed programs. Prints every problem it finds, then stops with an error.

Run from the project folder:  python3 scripts/check_setup.py
"""
import csv
import glob
import shutil
import sys
from pathlib import Path

import yaml
from Bio import SeqIO

from plasmid_utils import read_single_sequence, reverse_complement

ROOT = Path(__file__).resolve().parents[1]
config = yaml.safe_load((ROOT / 'config' / 'config.yaml').read_text())
with (ROOT / 'samples.csv').open(newline='') as handle:
    samples = list(csv.DictReader(handle))
problems = []

if not samples:
    problems.append('samples.csv has no rows')

# --- each sample: design, seed, FASTQ ---------------------------------------------------
reads_folder = Path(str(config['reads_folder'])).expanduser()
for s in samples:
    try:
        design = read_single_sequence(ROOT / s['design_fasta'])
        seed = read_single_sequence(ROOT / s['seed_fasta'])
    except (ValueError, FileNotFoundError) as error:
        problems.append(f"{s['sample']}: {error}")
        continue
    if len(design) != int(s['expected_length_bp']):
        problems.append(f"{s['sample']}: design is {len(design)} bp but samples.csv says {s['expected_length_bp']}")
    if len(seed) < 500 or design.count(seed) + design.count(reverse_complement(seed)) != 1:
        problems.append(f"{s['sample']}: seed must be >= 500 bp and occur exactly once in the design")
    pattern = str(reads_folder / config['reads_filename'].format(barcode=s['barcode'], sample=s['sample']))
    fastq_files = [f for f in glob.glob(pattern) if Path(f).is_file()]
    if not fastq_files:
        problems.append(f'no FASTQ files found for {s["sample"]}. Pattern searched: {pattern}')
    elif all(Path(f).stat().st_size == 0 for f in fastq_files):
        problems.append(f'all FASTQ files are empty for {s["sample"]}: {pattern}')
    else:
        megabytes = sum(Path(f).stat().st_size for f in fastq_files) / 1e6
        print(f'found {s["sample"]}: {len(fastq_files)} FASTQ file(s), {megabytes:.1f} MB')

# --- insert files -----------------------------------------------------------------------
constructs = {s['construct']: int(s['expected_length_bp']) for s in samples}
with (ROOT / 'config' / 'insert_coordinates.csv').open(newline='') as handle:
    for row in csv.DictReader(handle):
        length = constructs.get(row['construct'])
        try:
            start, end = int(row['start_1based']), int(row['end_1based'])
            if length is None or not (1 <= start <= length and 1 <= end <= length):
                raise ValueError
        except ValueError:
            problems.append(f'insert_coordinates.csv: invalid row {row}')

insert_fasta = str(config.get('insert_sequences_fasta') or '')
if insert_fasta:
    path = ROOT / insert_fasta
    if not path.is_file():
        problems.append(f'insert_sequences_fasta not found: {path}')
    else:
        for record in SeqIO.parse(path, 'fasta'):
            owner, separator, _ = record.id.partition('|')
            if separator and owner not in constructs:
                problems.append(f'insert FASTA header names an unknown construct: {record.id}')

# --- programs ---------------------------------------------------------------------------
programs = ['snakemake', 'perl', 'kalign', 'minimap2', 'samtools', 'bcftools', 'blastn']
if config.get('run_sniffles', True):
    programs.append('sniffles')
if config.get('run_flye', False):
    programs.append('flye')
# NANOPlasty is an upstream Perl script, not a Conda executable.
script_path = Path(str(config.get('nanoplasty_script', 'tools/NANOPlasty/NANOPlasty0.4.pl')))
if not script_path.is_absolute():
    script_path = ROOT / script_path
if not script_path.is_file():
    problems.append(f'NANOPlasty script not found: {script_path} (clone official repository first)')
for module in ('MCE::Child', 'MCE::Channel', 'Parallel::ForkManager'):
    import subprocess
    check = subprocess.run(['perl', '-M' + module, '-e', '1'], capture_output=True, text=True) if shutil.which('perl') else None
    if check is None or check.returncode != 0:
        problems.append(f'Perl module unavailable: {module}')

for program in programs:
    if not shutil.which(program):
        problems.append(f'program not installed (or environment not activated): {program}')

if problems:
    for problem in problems:
        print('PROBLEM:', problem, file=sys.stderr)
    sys.exit(1)
print(f'All good: {len(samples)} samples, designs, seeds, FASTQ files, insert files and programs.')
