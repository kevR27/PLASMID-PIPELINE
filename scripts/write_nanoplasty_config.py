#!/usr/bin/env python3
"""Write a NANOPlasty 0.4 config for one demultiplexed ONT library.

The official NANOPlasty configuration has separate long-read sections; only
Nanopore input is populated here. The `mito` type is its documented small-
genome mode, NOT a classification of plasmid sequence as mitochondrial DNA.
"""
import argparse
from pathlib import Path


def yes_no(v):
    if v.lower() not in ('yes', 'no'):
        raise argparse.ArgumentTypeError('Must be yes or no')
    return v.lower()


def run(a):
    if a.design_length <= 0 or not (0.02 <= a.range_fraction <= 0.40):
        raise ValueError('Invalid expected length or genome range fraction')
    if a.min_read_length < 100:
        raise ValueError('Minimum read length must be >= 100 bp')
    seed = Path(a.seed).expanduser().resolve(strict=True)
    reads = Path(a.reads_folder).expanduser().resolve(strict=True)
    if not (reads / 'reads.fastq').is_file():
        raise FileNotFoundError(f'Expected a single prefiltered FASTQ: {reads / "reads.fastq"}')
    native = Path(a.output_folder).expanduser().absolute()
    if native.name != 'native':
        raise ValueError('NANOPlasty output must be a dedicated native/ folder')
    low = round(a.design_length * (1 - a.range_fraction))
    high = round(a.design_length * (1 + a.range_fraction))
    text = f'''Project:
-----------------------
Project name          = {a.sample}
Type                  = mito
Genome Range          = {low}-{high}
Save assembled reads  = {a.save_assembled_reads}
Seed Input            = {seed}
Extend seed directly  = {a.extend_seed_directly}
Threads               = {a.threads}
Chloroplast sequence  =
Output path           = {native}/

Nanopore reads:
-----------------------
Nanopore reads        = {reads}/
Local DB and NP reads =
R10                   = yes
Min read length NP    = {a.min_read_length}
Use full dataset      = {a.use_full_dataset}

PacBio reads:
-----------------------
PacBio reads          =
Local DB and PB reads =
Min read length PB    =
Use full dataset      =
'''
    output = Path(a.config_file)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(text)
    print(f'Wrote NANOPlasty configuration: {output} (genome range {low}-{high} bp)')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample', required=True)
    p.add_argument('--seed', required=True)
    p.add_argument('--design-length', required=True, type=int)
    p.add_argument('--reads-folder', required=True)
    p.add_argument('--output-folder', required=True)
    p.add_argument('--config-file', required=True)
    p.add_argument('--threads', required=True, type=int)
    p.add_argument('--min-read-length', type=int, default=3000)
    p.add_argument('--range-fraction', type=float, default=0.15)
    p.add_argument('--use-full-dataset', type=yes_no, default='no')
    p.add_argument('--extend-seed-directly', type=yes_no, default='no')
    p.add_argument('--save-assembled-reads', type=yes_no, default='no')
    run(p.parse_args())


if __name__ == '__main__':
    main()
