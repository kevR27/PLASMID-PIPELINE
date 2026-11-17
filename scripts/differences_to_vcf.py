#!/usr/bin/env python3
"""Write the assembly-vs-design differences as a VCF file.

This VCF describes how the ASSEMBLY differs from the DESIGN. It is not a read-based
variant call: there is no QUAL and no genotype. The Snakefile then left-aligns it with
`bcftools norm` so equivalent indel positions (homopolymers, repeats) are written once.
Differences touching the artificial start/end of the circular FASTA get the flag
CIRCULAR_EDGE: check them by eye.
"""
import argparse
from pathlib import Path

from plasmid_utils import read_single_sequence, read_table
from Bio import SeqIO


def vcf_record(design, event):
    """Convert one difference into (position, REF, ALT, is_edge) using VCF conventions."""
    kind = event['type']
    pos = int(event['position_1based'])
    design_part = event['design_sequence'].upper()
    assembly_part = event['assembly_sequence'].upper()
    edge = False

    if kind in ('SNV', 'MNV'):
        if not (design_part and assembly_part and 1 <= pos <= len(design)
                and design[pos - 1:pos - 1 + len(design_part)] == design_part):
            raise ValueError(f'Invalid substitution: {event}')
        return pos, design_part, assembly_part, False

    if kind == 'INS':
        if not assembly_part or not 0 <= pos <= len(design):
            raise ValueError(f'Invalid insertion: {event}')
        if pos == 0:   # inserted before the first base: anchor on the first base instead
            return 1, design[0], assembly_part + design[0], True
        return pos, design[pos - 1], design[pos - 1] + assembly_part, pos == len(design)

    if kind == 'DEL':
        if not design_part or not 1 <= pos <= len(design) or design[pos - 1:pos - 1 + len(design_part)] != design_part:
            raise ValueError(f'Invalid deletion: {event}')
        edge = pos + len(design_part) - 1 == len(design)
        if pos > 1:    # normal case: anchor on the base before the deletion
            return pos - 1, design[pos - 2] + design_part, design[pos - 2], edge
        if len(design_part) < len(design):   # deletion at the very start: anchor on the base after it
            after = design[len(design_part)]
            return 1, design_part + after, after, True
        raise ValueError('The whole design is deleted; cannot be written as a VCF indel')

    raise ValueError(f'Unknown difference type {kind}')


def write_vcf(design_fasta, differences_tsv, vcf_path):
    record = next(SeqIO.parse(design_fasta, 'fasta'))
    design = read_single_sequence(design_fasta)
    rows = []
    for event in read_table(differences_tsv):
        pos, ref, alt, edge = vcf_record(design, event)
        if ref == alt or not ref or not alt:
            raise ValueError(f'Invalid VCF alleles for {event}')
        info = f"SOURCE=ASSEMBLY_VS_DESIGN;KIND={event['type']}" + (';CIRCULAR_EDGE' if edge else '')
        rows.append((pos, ref, alt, info))
    rows.sort(key=lambda r: r[0])

    path = Path(vcf_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w') as out:
        out.write('##fileformat=VCFv4.2\n')
        out.write('##source=plasmid_pipeline_assembly_vs_design\n')
        out.write('##INFO=<ID=SOURCE,Number=1,Type=String,Description="Assembly compared with design, not a read-based call">\n')
        out.write('##INFO=<ID=KIND,Number=1,Type=String,Description="Type of difference">\n')
        out.write('##INFO=<ID=CIRCULAR_EDGE,Number=0,Type=Flag,Description="Touches the artificial start/end of the circular FASTA; check by eye">\n')
        out.write(f'##contig=<ID={record.id},length={len(design)}>\n')
        out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for pos, ref, alt, info in rows:
            out.write(f'{record.id}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\n')
    return len(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    for name in ('design', 'differences', 'output'):
        parser.add_argument('--' + name, required=True)
    args = parser.parse_args()
    write_vcf(args.design, args.differences, args.output)


if __name__ == '__main__':
    main()
