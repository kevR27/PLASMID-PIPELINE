#!/usr/bin/env python3
"""Patch ONLY the ONT forward-query extraction in NANOPlasty0.4.pl.

The original Perl source has the same extraction expression in the PacBio
and Nanopore branches. Scope the replacement to LONGER_LAST_600_NP.
This is a local defensive workaround, not an upstream-tested fix.
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

OLD = 'my $last_600 = substr $read, -($first_length_back+$first_length_back_extra), $first_length_back;'
MARKER = '# ONT local guard: prevent empty query_NP BLAST inputs'

# Keep the original source indentation at the replacement site.
NEW = '''# ONT local guard: prevent empty query_NP BLAST inputs
my $np_available = length($read);
die "NANOPlasty ONT: current assembly extension sequence has only $np_available bp at iteration $y\\n"
    if $np_available < 100;
if ($first_length_back < 100 || $first_length_back > $np_available)
{
    warn "NANOPlasty ONT: limiting query length from $first_length_back to $np_available at iteration $y\\n";
    $first_length_back = $np_available;
}
if ($first_length_back_extra < 0)
{
    warn "NANOPlasty ONT: resetting negative query offset ($first_length_back_extra) at iteration $y\\n";
    $first_length_back_extra = 0;
}
my $np_max_extra = $np_available - $first_length_back;
if ($first_length_back_extra > $np_max_extra)
{
    warn "NANOPlasty ONT: limiting extra query offset from $first_length_back_extra to $np_max_extra at iteration $y\\n";
    $first_length_back_extra = $np_max_extra;
}
my $last_600 = substr $read,
    ($np_available - $first_length_back - $first_length_back_extra),
    $first_length_back;
die "NANOPlasty ONT: generated query has only ".length($last_600)." bp at iteration $y\\n"
    if length($last_600) < 100;'''


def patch(source: str) -> str:
    label = 'LONGER_LAST_600_NP:'
    end_label = 'SKIP_TO_CONFIRMED_NP:'
    if source.count(label) != 1 or source.count(end_label) != 1:
        raise ValueError('Expected exactly one Nanopore query block and one end label; script version may differ')
    start = source.index(label)
    end = source.index(end_label, start)
    block = source[start:end]
    if MARKER in block:
        return source
    count = block.count(OLD)
    if count != 1:
        raise ValueError(f'Expected exactly one unpatched expression within Nanopore query block; found {count}')
    # Replace only inside the Nanopore block, leaving the analogous PacBio line unchanged.
    new_block = block.replace(OLD, NEW, 1)
    return source[:start] + new_block + source[end:]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('script', type=Path, help='Installed NANOPlasty Perl script')
    args = ap.parse_args()
    file = args.script.expanduser().resolve()
    if not file.is_file():
        ap.error(f'Script not found: {file}')
    content = file.read_text()
    try:
        patched = patch(content)
    except ValueError as e:
        ap.error(f'{e}. No changes made to {file}')
    if patched == content:
        print(f'Already patched: {file}')
        return
    backup = Path(str(file) + '.bak')
    if backup.exists():
        ap.error(f'Backup already exists: {backup}; no changes made')
    shutil.copy2(file, backup)
    try:
        file.write_text(patched)
    except Exception:
        shutil.copy2(backup, file)
        raise
    print(f'Patched Nanopore query extraction: {file}')
    print(f'Original backup: {backup}')
    print('Next: perl -c "${CONDA_PREFIX}/bin/NANOPlasty0.4.pl"')


if __name__ == '__main__':
    main()
