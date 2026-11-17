#!/usr/bin/env python3
"""Three figures per sample: read lengths, read depth along the design, and a circular
map of where the assembly differs from the design. Each is saved as PNG and PDF."""
import argparse
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from plasmid_utils import read_table


def save(fig, stem):
    fig.tight_layout()
    Path(stem).parent.mkdir(parents=True, exist_ok=True)
    for extension in ('png', 'pdf'):
        fig.savefig(f'{stem}.{extension}', dpi=220, bbox_inches='tight')
    plt.close(fig)


def plot_read_lengths(args, stem):
    rows = read_table(args.length_histogram)
    centers = [int(r['length_bin_start_bp']) + 50 for r in rows]
    fig, ax = plt.subplots(figsize=(8, 3.5))
    ax.bar(centers, [int(r['all_reads']) for r in rows], width=95, alpha=.45, label='All reads')
    ax.step(centers, [int(r['reads_kept']) for r in rows], where='mid', linewidth=1.5, label='Kept after filtering')
    ax.set(xlabel='Read length (bp, 100 bp bins)', ylabel='Number of reads', title=f'{args.sample}: read lengths')
    ax.legend(frameon=False)
    save(fig, stem)


def plot_depth(args, stem, design_length):
    positions, depths = [], []
    with open(args.depth) as handle:
        for line in handle:
            columns = line.split('\t')
            if len(columns) >= 3:
                positions.append(int(columns[1]))
                depths.append(int(columns[2]))
    fig, ax = plt.subplots(figsize=(9, 3.5))
    ax.plot(positions, depths, linewidth=.8)
    ax.set(xlabel='Position in design (bp)', ylabel='Read depth (x)', xlim=(1, design_length),
           title=f'{args.sample}: read depth')
    save(fig, stem)


def plot_differences(args, stem, design_length):
    differences = read_table(args.differences)
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection='polar')
    ax.plot([i * 2 * math.pi / design_length for i in range(design_length)], [1] * design_length,
            linewidth=7, alpha=.25)
    by_type = {}
    for row in differences:
        by_type.setdefault(row['type'], []).append(int(row['position_1based']))
    for kind, positions in by_type.items():
        angles = [(p - 1) % design_length * 2 * math.pi / design_length for p in positions]
        ax.scatter(angles, [1.0] * len(positions), s=28, label=kind, zorder=3)
    ax.set_ylim(.7, 1.3)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.grid(False)
    ax.spines['polar'].set_visible(False)
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location('N')
    ax.set_title(f'{args.sample}: assembly differences\nfrom the design', pad=20)
    if by_type:
        ax.legend(loc='lower center', bbox_to_anchor=(.5, -.13), ncol=4, frameon=False)
    else:
        ax.text(.5, .5, 'No differences', transform=ax.transAxes, ha='center')
    save(fig, stem)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    for name in ('sample', 'read-stats', 'length-histogram', 'depth', 'differences', 'output-folder'):
        parser.add_argument('--' + name, required=True)
    args = parser.parse_args()
    design_length = int(read_table(args.read_stats)[0]['design_length_bp'])
    plot_read_lengths(args, f'{args.output_folder}/read_lengths')
    plot_depth(args, f'{args.output_folder}/depth', design_length)
    plot_differences(args, f'{args.output_folder}/differences', design_length)


if __name__ == '__main__':
    main()
