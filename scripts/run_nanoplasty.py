#!/usr/bin/env python3
"""Run NANOPlasty without destroying earlier native assemblies.

Every invocation assembles in a uniquely named attempts/run_*/native directory
with a copy of the NANOPlasty config whose Output path is redirected there.
Only a validated candidate is promoted to native/ and to the Snakemake outputs.
Previous native folders and failed attempts are retained for review.
"""
import argparse
import csv
from datetime import datetime, timezone
import fcntl
import os
from pathlib import Path
import re
import shutil
import subprocess
import uuid

from Bio import SeqIO
from Bio.Seq import Seq

EXCLUDE = ('assembled_read', 'read_', 'tmp', 'seed', 'query', 'reference',
           'quality', 'blast_db', 'database', 'partial')


def find_candidates(native, expected):
    possible = []
    for path in sorted(native.rglob('*')):
        if not path.is_file() or path.suffix.lower() not in ('.fa', '.fna', '.fasta'):
            continue
        name = path.name.lower()
        if any(word in name for word in EXCLUDE):
            continue
        if any(word in part.lower() for part in path.relative_to(native).parts[:-1]
               for word in ('blast', 'tmp')):
            continue
        if not any(word in name for word in ('assembly', 'contig', 'circular')):
            continue
        try:
            with path.open() as handle:
                records = list(SeqIO.parse(handle, 'fasta'))
        except (ValueError, OSError):
            continue
        if len(records) != 1:
            continue
        seq = str(records[0].seq).upper()
        if not seq or set(seq) - set('ACGTN'):
            continue
        if not (0.70 * expected <= len(seq) <= 1.30 * expected):
            continue
        if name.startswith('assembly_final_'):
            score = 0
        elif name.startswith('assembly_') and re.search(r'_option_\d+\.(fa|fna|fasta)$', name):
            score = 1
        elif name.startswith('assembly_raw_'):
            score = 2
        elif 'circular' in name:
            score = 3
        else:
            score = 4
        possible.append((score, path, records[0], len(seq)))
    return sorted(possible, key=lambda q: (q[0], str(q[1])))


def circularly_equivalent(seq_a, seq_b):
    if len(seq_a) != len(seq_b) or not seq_a:
        return False
    if seq_a in (seq_b + seq_b):
        return True
    rc = str(Seq(seq_b).reverse_complement())
    return seq_a in (rc + rc)


def capture(native, expected, sample, output_fasta, source_report):
    """Select and write a valid native assembly; no destructive operations."""
    native = Path(native)
    choices = find_candidates(native, expected)
    if not choices:
        files = [str(x.relative_to(native)) for x in native.rglob('*') if x.is_file()]
        raise RuntimeError('No near-expected-length NANOPlasty assembly FASTA was found. '
                           'Inspect native output and log. Files: ' + repr(files[:30]))
    best_score = choices[0][0]
    best = [choice for choice in choices if choice[0] == best_score]
    if len(best) != 1:
        raise RuntimeError('Ambiguous NANOPlasty final assemblies; manual review required: ' +
                           repr([(str(x[1]), x[3]) for x in best]))
    _, source, record, length = best[0]
    selected_seq = str(record.seq).upper()
    alternative_status = []
    for _, path, alt_record, alt_length in choices:
        if path == source:
            continue
        alt_seq = str(alt_record.seq).upper()
        status = ('circularly_equivalent'
                  if circularly_equivalent(selected_seq, alt_seq) else 'different')
        alternative_status.append(f'{path.name}:{status}:{alt_length}bp')
        if status == 'different':
            print(f'WARNING: {path.name} differs from selected final assembly; '
                  'review native outputs before interpreting variants.', flush=True)
    output = Path(output_fasta)
    report = Path(source_report)
    output.parent.mkdir(parents=True, exist_ok=True)
    report.parent.mkdir(parents=True, exist_ok=True)
    with output.open('w') as f:
        f.write(f'>{sample}_nanoplasty_from_{record.id}\n')
        for i in range(0, len(selected_seq), 80):
            f.write(selected_seq[i:i + 80] + '\n')
    with report.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=('sample', 'assembler', 'source_fasta',
                                         'source_record', 'length_bp',
                                         'expected_length_bp', 'alternative_native_stages'),
                           delimiter='\t')
        w.writeheader()
        w.writerow({'sample': sample, 'assembler': 'NANOPlasty',
                    'source_fasta': str(source.resolve()), 'source_record': record.id,
                    'length_bp': length, 'expected_length_bp': expected,
                    'alternative_native_stages': ';'.join(alternative_status)})
    print(f'Selected NANOPlasty contig: {source} ({length} bp)', flush=True)
    return source, length


def redirect_output_config(config, expected_native, target_native, destination):
    """Rewrite only the Output path of a temporary config, never the original."""
    text = config.read_text()
    pattern = re.compile(r'(?m)^(\s*Output path\s*=\s*)([^\r\n]*)$')
    matches = list(pattern.finditer(text))
    if len(matches) != 1:
        raise ValueError('NANOPlasty config must contain exactly one Output path = line')
    existing = Path(matches[0].group(2).strip()).expanduser().resolve()
    if existing != expected_native:
        raise ValueError(f'Config Output path {existing} does not match '
                         f'--native-folder {expected_native}. Refusing cross-sample write.')
    rewritten = pattern.sub(lambda m: m.group(1) + str(target_native) + '/', text, count=1)
    destination.write_text(rewritten)


def rewrite_report_source(report, old_native, new_native):
    """Update source_fasta after moving a validated attempt into native/."""
    with report.open(newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        columns = reader.fieldnames
        rows = list(reader)
    if columns is None or len(rows) != 1:
        raise ValueError(f'Unexpected NANOPlasty source report: {report}')
    old_source = Path(rows[0]['source_fasta'])
    rel = old_source.relative_to(old_native)
    rows[0]['source_fasta'] = str(new_native / rel)
    with report.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def preserve_previous_outputs(output_fasta, source_report, attempt):
    """Keep a recoverable copy outside Snakemake's declared output paths."""
    previous = attempt / 'previous_outputs'
    for path in (Path(output_fasta), Path(source_report)):
        if path.is_file():
            previous.mkdir(exist_ok=True)
            shutil.copy2(path, previous / path.name)
            print(f'Preserved prior pipeline output: {previous / path.name}', flush=True)


def run(a):
    script = Path(a.script).expanduser().resolve(strict=True)
    config = Path(a.config_file).expanduser().resolve(strict=True)
    native = Path(a.native_folder).expanduser().resolve()
    if native.name != 'native' or native.parent.name != 'nanoplasty':
        raise ValueError(f'Refusing unexpected native output directory: {native}')
    if native.is_symlink():
        raise ValueError('Refusing symlinked native output directory')
    root = native.parent
    root.mkdir(parents=True, exist_ok=True)
    attempts_root = root / 'attempts'
    attempts_root.mkdir(exist_ok=True)

    # A second Snakemake session must not promote a competing result.
    with (root / '.nanoplasty_run.lock').open('a+') as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f'Another NANOPlasty attempt is active for {a.sample}') from exc

        run_id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ') + '_' + uuid.uuid4().hex[:8]
        attempt = attempts_root / ('run_' + run_id)
        attempt.mkdir(exist_ok=False)
        attempt_native = attempt / 'native'
        attempt_native.mkdir()
        staged_config = attempt / 'nanoplasty_config.txt'
        saved_log = attempt / 'nanoplasty_stdout.log'
        status = attempt / 'status.txt'
        status.write_text('RUNNING\n')
        preserve_previous_outputs(a.output_fasta, a.source_report, attempt)
        print(f'NANOPlasty attempt directory: {attempt}', flush=True)

        try:
            redirect_output_config(config, native, attempt_native, staged_config)
            cmd = ['perl', str(script), '-c', str(staged_config)]
            print('Running:', ' '.join(cmd), flush=True)
            # Keep the entire assembler output as a permanent per-attempt log.
            with saved_log.open('w') as log, subprocess.Popen(
                    cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    text=True, bufsize=1) as proc:
                for line in proc.stdout:
                    print(line, end='', flush=True)
                    log.write(line)
                returncode = proc.wait()
            if returncode:
                raise RuntimeError(f'NANOPlasty exited with code {returncode}; '
                                   f'inspect {saved_log}')

            # Validate into private staged files. The public outputs are unchanged
            # until a real, near-expected-length assembly has been selected.
            staged_fasta = attempt / 'selected_assembly.fasta'
            staged_report = attempt / 'selected_source.tsv'
            selected, length = capture(attempt_native, a.design_length, a.sample,
                                       staged_fasta, staged_report)
            rewrite_report_source(staged_report, attempt_native, native)

            # Promote only after validation. Archive the previous native folder.
            previous_native = None
            if native.exists():
                if not native.is_dir():
                    raise ValueError(f'Expected directory at {native}')
                previous_native = attempts_root / ('previous_native_' + run_id)
                native.rename(previous_native)
                print(f'Archived previous native output: {previous_native}', flush=True)
            try:
                attempt_native.rename(native)
            except Exception:
                if previous_native is not None:
                    previous_native.rename(native)
                raise

            output = Path(a.output_fasta)
            report = Path(a.source_report)
            output.parent.mkdir(parents=True, exist_ok=True)
            report.parent.mkdir(parents=True, exist_ok=True)
            os.replace(staged_fasta, output)
            os.replace(staged_report, report)
            status.write_text(f'SUCCESS\nlength_bp={length}\ncurrent_native={native}\n')
            print(f'Assembly promoted to {native}; outputs: {output}, {report}', flush=True)
        except Exception as exc:
            status.write_text(f'FAILED\n{type(exc).__name__}: {exc}\n')
            print(f'FAILED NANOPlasty attempt retained for review: {attempt}', flush=True)
            raise


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--script', help='Path to NANOPlasty Perl script')
    p.add_argument('--config-file', help='NANOPlasty configuration file')
    p.add_argument('--native-folder', required=True)
    p.add_argument('--sample', required=True)
    p.add_argument('--design-length', type=int, required=True)
    p.add_argument('--output-fasta', required=True)
    p.add_argument('--source-report', required=True)
    p.add_argument('--capture-only', action='store_true',
                   help='Capture existing native/ without running NANOPlasty')
    a = p.parse_args()
    if a.capture_only:
        native = Path(a.native_folder).expanduser().resolve(strict=True)
        capture(native, a.design_length, a.sample, a.output_fasta, a.source_report)
    else:
        if not a.script or not a.config_file:
            p.error('--script and --config-file are required unless --capture-only is used')
        run(a)


if __name__ == '__main__':
    main()
