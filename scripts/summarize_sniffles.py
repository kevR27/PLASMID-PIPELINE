#!/usr/bin/env python3
"""Count Sniffles2 structural-variant candidates and flag those at an artificial end.

The designed plasmid is a linear FASTA of a circular molecule, so reads crossing the
FASTA origin are split and can create spurious calls there. Sniffles2 is therefore run on
two references (original and half-rotated). A candidate within --edge-window bp of either
end of the reference it was called on is counted as 'near an artificial end'; real SVs
should appear away from the ends in at least one of the two runs.
Candidates are not classified as true SVs.
"""
import argparse
import re

from plasmid_utils import write_table


def scan(vcf_path, length, window):
    """Return (total, pass_count, near_end_count) for a Sniffles2 VCF."""
    total = passed = near = 0
    with open(vcf_path) as handle:
        for line in handle:
            if line.startswith('##contig') and not length:
                match = re.search(r'length=(\d+)', line)
                length = int(match.group(1)) if match else 0
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            total += 1
            if len(fields) > 6 and fields[6] in ('PASS', '.'):
                passed += 1
            if not length:
                continue
            pos = int(fields[1])
            info = dict(kv.split('=', 1) for kv in (fields[7].split(';') if len(fields) > 7 else []) if '=' in kv)
            end = int(info['END']) if info.get('END', '').isdigit() else pos
            if min(pos, end) <= window or max(pos, end) >= length - window:
                near += 1
    return total, passed, near


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--design-vcf', required=True, help='Sniffles2 VCF from reads aligned to the design')
    parser.add_argument('--rotated-vcf', default='', help='Sniffles2 VCF from reads aligned to the rotated design')
    parser.add_argument('--length', type=int, default=0, help='plasmid length (default: from VCF header)')
    parser.add_argument('--edge-window', type=int, default=100)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    total, passed, near = scan(args.design_vcf, args.length, args.edge_window)
    row = dict(calls_on_design=total, pass_on_design=passed, near_ends_of_design=near)
    if args.rotated_vcf:
        r_total, r_passed, r_near = scan(args.rotated_vcf, args.length, args.edge_window)
        row.update(calls_on_rotated_design=r_total, pass_on_rotated_design=r_passed,
                   near_ends_of_rotated_design=r_near)
    row['note'] = 'RAW_READ_SV_CANDIDATES_NOT_CONFIRMED_VARIANTS'
    write_table(args.output, [row], list(row))


if __name__ == '__main__':
    main()
