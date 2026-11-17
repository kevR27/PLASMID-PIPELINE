#!/usr/bin/env python3
"""One overview figure for all samples: identity to the design, and read depth."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from plasmid_utils import read_table


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--summary', required=True)
    parser.add_argument('--prefix', required=True)
    args = parser.parse_args()

    rows = read_table(args.summary)
    x = list(range(len(rows)))
    identities = [float(r['assembly_identity_percent']) for r in rows]
    depths = [float(r['mapped_mean_depth']) if r['mapped_mean_depth'] != 'NA' else 0 for r in rows]

    fig, (top, bottom) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    top.bar(x, identities)
    top.set_ylim(max(0, min(identities) - 1), 100.2)
    top.set_ylabel('Identity to design (%)')
    top.set_title('Assembled plasmids compared with their designs')
    bottom.bar(x, depths)
    bottom.set_ylabel('Mean read depth (x)')
    bottom.set_xticks(x, [r['barcode'] for r in rows], rotation=35, ha='right')
    # dotted line wherever the construct changes
    for i in range(1, len(rows)):
        if rows[i]['construct'] != rows[i - 1]['construct']:
            for ax in (top, bottom):
                ax.axvline(i - .5, linestyle=':', alpha=.45)
    fig.tight_layout()
    Path(args.prefix).parent.mkdir(parents=True, exist_ok=True)
    for extension in ('png', 'pdf'):
        fig.savefig(f'{args.prefix}.{extension}', dpi=220, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()
