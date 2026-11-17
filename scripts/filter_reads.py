#!/usr/bin/env python3
"""Step 1: filter the reads of one barcode and write simple read statistics.

Input: --fastq-pattern, a file name or a pattern such as  /data/run/fastq_pass/barcode01/*.fastq.gz
All files matching the pattern are read one after another and treated as ONE set of reads
(so a barcode folder with hundreds of FASTQ files needs no manual concatenation).
Files can be .fastq or .fastq.gz, even mixed.

Reads are kept if they are at least --min-length bp long and (optionally) have a mean
quality of at least --min-quality. The filtered reads are written as one plain
(uncompressed) FASTQ because that is the assembly tools read reliably.

Outputs
  --output-fastq     the filtered reads
  --stats            one-line table: read counts, bases, N50, mean quality, estimated depth
  --length-histogram read lengths in 100 bp bins (all reads vs kept reads), used for plots
"""
import argparse
import glob
import gzip
import os
from collections import Counter
from pathlib import Path

from plasmid_utils import write_table


def open_reads(path):
    """Open a FASTQ file, compressed or not (decided by the file's first two bytes)."""
    with open(path, 'rb') as handle:
        is_gzip = handle.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt') if is_gzip else open(path, 'rt')


def find_read_files(pattern):
    """All files matching the pattern (~ allowed), sorted so the result is reproducible."""
    files = sorted(f for f in glob.glob(os.path.expanduser(pattern)) if os.path.isfile(f))
    if not files:
        raise FileNotFoundError(f'No FASTQ files found for pattern: {pattern}')
    return files


def read_records(files):
    """Yield (header, sequence, plus_line, quality) from every file in turn."""
    for path in files:
        with open_reads(path) as handle:
            while True:
                header = handle.readline()
                if not header:
                    break
                sequence = handle.readline().strip()
                plus = handle.readline()
                quality = handle.readline().strip()
                if not header.startswith('@') or not plus.startswith('+') or len(sequence) != len(quality):
                    raise ValueError(f'Malformed FASTQ record in {path} (near header {header.strip()[:60]})')
                yield header, sequence, plus, quality


def n50(length_counts, total_bases):
    """N50 from a {read_length: count} table."""
    running = 0
    for length, count in sorted(length_counts.items(), reverse=True):
        running += length * count
        if running >= total_bases / 2:
            return length
    return 0


def filter_reads(args):
    reads_total = bases_total = reads_kept = bases_kept = 0
    quality_sum = 0.0
    all_bins, kept_bins, kept_lengths = Counter(), Counter(), Counter()
    Path(args.output_fastq).parent.mkdir(parents=True, exist_ok=True)

    files = find_read_files(args.fastq_pattern)
    with open(args.output_fastq, 'w') as reads_out:
        for header, sequence, plus, quality in read_records(files):
            length = len(sequence)
            reads_total += 1
            bases_total += length
            length_bin = length // 100 * 100
            all_bins[length_bin] += 1
            mean_quality = sum(ord(c) - 33 for c in quality) / length if length else 0

            if length >= args.min_length and mean_quality >= args.min_quality:
                reads_out.write(f'{header}{sequence}\n{plus}{quality}\n')
                reads_kept += 1
                bases_kept += length
                quality_sum += mean_quality
                kept_lengths[length] += 1
                kept_bins[length_bin] += 1

    if reads_kept == 0:
        raise ValueError(f'No reads of {args.sample} passed the filters; lower min_read_length / min_mean_read_quality')

    stats = dict(
        sample=args.sample, barcode=args.barcode, design_length_bp=args.design_length,
        input_files=len(files), all_reads=reads_total, all_bases=bases_total,
        reads_kept=reads_kept, bases_kept=bases_kept,
        kept_percent=round(100 * reads_kept / reads_total, 2),
        mean_length_kept=round(bases_kept / reads_kept, 1),
        n50_kept=n50(kept_lengths, bases_kept),
        mean_read_quality_kept=round(quality_sum / reads_kept, 2),
        estimated_depth=round(bases_kept / args.design_length, 2),   # bases / plasmid size
        min_read_length=args.min_length, min_mean_read_quality=args.min_quality)
    write_table(args.stats, [stats], list(stats))

    histogram = [dict(length_bin_start_bp=b, all_reads=all_bins[b], reads_kept=kept_bins[b])
                 for b in sorted(set(all_bins) | set(kept_bins))]
    write_table(args.length_histogram, histogram, ['length_bin_start_bp', 'all_reads', 'reads_kept'])


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    for name in ('fastq-pattern', 'output-fastq', 'stats', 'length-histogram', 'sample', 'barcode'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--design-length', type=int, required=True)
    parser.add_argument('--min-length', type=int, default=1000)
    parser.add_argument('--min-quality', type=float, default=0)
    filter_reads(parser.parse_args())


if __name__ == '__main__':
    main()
