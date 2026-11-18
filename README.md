# Plasmid assembly and validation (ONT R10.4.1)

Nine demultiplexed Nanopore samples, three designed constructs (6,655 / 6,164 / 8,280 bp), one sample per barcode. The primary assembler is **NANOPlasty0.4.pl** (seed-assisted long-read assembly), with **Flye** as an independent cross-check. The original Flye contig is preserved and potential tandem-copy assemblies are evaluated before a one-unit sequence is used in comparisons.

> **Validation boundary:** NANOPlasty's documented `Type` values are `mito`, `chloro`, and `mito_plant`, not `plasmid`. This workflow uses `Type = mito` with an expected *plasmid-size* `Genome Range`, because the program supports small circular ONT genomes. That is a pragmatic, **untested use outside its explicitly named target types**. NANOPlasty has produced both successful and failed circularization attempts for barcode 06 on the user's server; the full nine-barcode workflow has **not been confirmed end-to-end**. Confirm the first sample's output before scaling up.

## Standalone GitHub repository

This workflow is maintained in its independent repository: **https://github.com/kevR27/PLASMID-PIPELINE**. It is separate from `kevR27/PRIMER`. The repository is **public**; never commit patient data, unpublished confidential sequences, original FASTQs or server credentials. Generated `results/`, `work/`, FASTQs, BAMs and `.snakemake/` remain untracked; see `.gitignore`. The included construct reference and seed FASTAs are project data and must be reviewed by the repository owner before public distribution.

The repository snapshot includes the preservation-aware `run_nanoplasty.py`. It does not include the third-party NANOPlasty Perl source. The optional local patch in `patches/` guards against an empty ONT BLAST query but does **not** resolve every circularization/finalization failure and is not an upstream patch.

**Configuration note:** The included `config/config.yaml` contains an example FASTQ input directory. Change `reads_folder` and NANOPlasty options for your own environment before running. Avoid committing private sequencing data or personal credentials.

## Installation on your Linux server

Clone the repository to your preferred project directory (back up existing code first; no need to delete the original demultiplexed FASTQs). From the project folder:

```bash
conda env update -n plasmid_pipeline -f environment.yaml --prune
conda activate plasmid_pipeline
mkdir -p tools
git clone --depth 1 https://github.com/ndierckx/NANOPlasty.git tools/NANOPlasty
```

NANOPlasty 0.4 is maintained as a separate upstream Perl script; it is **not** bundled in this archive. Dependencies (BLAST, Kalign, Perl, MCE and Parallel::ForkManager) are declared in `environment.yaml`. Verify:

```bash
command -v kalign
perl -MMCE::Child -e 1
perl -MMCE::Channel -e 1
perl -MParallel::ForkManager -e 1
test -s tools/NANOPlasty/NANOPlasty0.4.pl
python3 scripts/check_setup.py
```

Set `reads_folder` in `config/config.yaml` to the absolute location of your demultiplexed FASTQ files on your system. Inputs matching `{barcode}/*.fastq.gz` are read by the filter script and written as a single `work/<sample>/filtered_reads/reads.fastq`. Thus **you do not need to convert the `.fastq.gz` to FASTA**. `samples.csv` refers to bundled designs and seeds by relative paths, which work in the server project folder.

## NANOPlasty parameters

The generated file `work/<sample>/nanoplasty_config.txt` follows the official `config_nanoplasty.txt` layout:

| Parameter | Value | Reason |
| --- | --- | --- |
| `Type` | `mito` | Only documented small-genome preset; experimental for plasmids |
| `Genome Range` | expected size ±5% | Restricts circular assembly around design length |
| `Seed Input` | respective clean 800-bp seed FASTA | Start from a unique sequence |
| `Extend seed directly` | `no` | Seed comes from the in-silico design, not a same-sample read |
| `R10` | `yes` | ONT R10.4.1 chemistry |
| `Min read length NP` | 2000 | Prefer longer molecules for assembly; separately configurable |
| `Use full dataset` | `yes` | Retain full dataset; change to `no` if reducing depth is appropriate |
| `Threads` | Snakemake allocated threads | Respect resource limits |

NANOPlasty has no NOVOLoci-style `Sequencing depth NP` option. The deprecated `max_depth_for_novoloci` option is removed. Settings can be modified in `config/config.yaml`. To change NANOPlasty settings after running, force the config rule or remove the generated `work/<sample>/nanoplasty_config.txt` and derived outputs before rerunning. Do not set `Extend seed directly = yes` merely to make assembly finish: that assumes a same-sample, exact seed and can conceal genuine mismatches.

### First sample: inspect before scaling up

```bash
cd /path/to/PLASMID-PIPELINE
python3 tests/test_scripts.py
python3 -m unittest discover -s tests -v

snakemake -s Snakefile \
  results/tal_merrf_right_bc06/nanoplasty/assembly_raw_tal_merrf_right_bc06.fasta \
  --cores 16 --use-conda \
  --conda-prefix /path/to/snakemake_envs \
  --rerun-incomplete -p
```

Inspect:

```bash
cat results/tal_merrf_right_bc06/logs/nanoplasty.log | tail -n 80
cat results/tal_merrf_right_bc06/nanoplasty/assembly_source.tsv
seqkit stats results/tal_merrf_right_bc06/nanoplasty/assembly_raw_tal_merrf_right_bc06.fasta
```

`seqkit` is optional; if unavailable, use `grep -v '^>' ... | tr -d '\n' | wc -c`. The separate output-adapter script searches NANOPlasty's *native* output directory for an unambiguous near-expected-length single-contig assembly. It **rejects** seed/temporary/read FASTAs or multiple equally plausible completed assemblies, rather than silently substituting one. If NANOPlasty uses a different final filename on your installation, inspect the native output and update `scripts/run_nanoplasty.py` conservatively.

**Important:** NANOPlasty has returned a 6,163-bp circularized barcode-06 assembly on the user's server, but other attempts completed without a final contig; full validation across all barcodes is incomplete. A nonzero return code, missing assembly, or ambiguous output correctly stops Snakemake, with native outputs preserved for diagnosis.

### Whole pipeline after first-sample success

```bash
snakemake -s Snakefile -n --cores 32
snakemake -s Snakefile --cores 32 --use-conda \
  --conda-prefix /path/to/snakemake_envs \
  --rerun-incomplete -p
```

The workflow is already a clean-run workflow: missing `results/` triggers the necessary assembly/validation steps. There is **no need to delete the original FASTQ data**, and avoid deleting existing useful results in future restarts. On NANOPlasty reruns, `scripts/run_nanoplasty.py` assembles in an isolated `attempts/run_*/native/` directory. It promotes successful results only after validation, archives a previous `native/`, and retains failed attempts and their logs. It does not erase a previously successful native assembly at the start of a retry.

## Pipeline outputs (one sample)

| Directory | Output and interpretation |
| --- | --- |
| `reads/` | Filtering statistics and length histogram |
| `nanoplasty/native/` | Most recently validated NANOPlasty native output |
| `nanoplasty/attempts/` | Retained failed attempts, run logs, and previously successful native directories |
| `nanoplasty/assembly_raw_<sample>.fasta` | Validated, single candidate extracted from native assembly |
| `nanoplasty/assembly_source.tsv` | Source filename, record ID and length provenance |
| `assembly/` | Primary assembly rotated/oriented to design (raw primary unchanged) |
| `assembly_vs_design/` | Sequence differences, assembly-derived VCF, read support and circular-origin read evidence |
| `alignments/` | BAMs, depth, coverage, IGV-ready alignments |
| `blast/` | Whole-plasmid and insert sequence checks |
| `sniffles/` | Raw-read SV candidates aligned to design and rotated design |
| `flye/assembly_raw/` | Unchanged Flye output (even when ~2× expected size) |
| `flye/assembly_normalized.fasta` | One-copy sequence **only if** tandem copies pass validation |
| `flye/copy_boundaries.tsv`, `copy_pairwise.tsv` | Boundary evidence and pairwise agreement across repeated units |
| `flye/vs_design.tsv`, `vs_nanoplasty.tsv` | Independent Flye comparisons |
| `plots/` | Per-sample read length/depth/difference visualizations |
| `logs/` | Native stderr/stdout for tools and normalization |

Overall outputs: `results/summary/all_samples.tsv`, `replicates.tsv`, `overview.png` and `overview.pdf`. In the summary, comparisons are named `flye_identity_to_nanoplasty_percent` and `flye_agrees_with_nanoplasty`. 

### Flye tandem normalization

The existing `scripts/normalize_flye_tandem.py` is retained. It identifies a reference-unique anchor repeated at approximately plasmid-length intervals, checks repeated units' global sequence identity, and only retains one observed unit when sufficiently concordant. It will fail rather than trim divergent copies. It does **not** replace observed bases with design bases and cannot distinguish assembly duplication from true biological tandem products. Keep the unnormalized Flye FASTA and review read alignments before concluding monomeric plasmid topology.

### Interpretation cautions

- NANOPlasty uses a seed; Flye does not. Agreement is independent corroboration, not proof of error-free assembly.
- The primary assembly-oriented sequence is compared to the designed reference. A differences VCF from an assembler is **not** equivalent to a raw-read variant call.
- Sniffles calls SVs above its configured length threshold; small SNVs/indels require assembly comparison and read-level inspection.
- `reads_across_origin.tsv` reports linearized-reference origin support; it is supporting sequencing evidence, not an independent biochemical demonstration of circularity.
- Final validation of each barcode should include the expected length, insert presence, agreement of independent assemblies, read support, and possible contamination or mixed plasmid populations.

## Testing scope

`tests/test_scripts.py` checks read filtering, config generation, alignment-independent comparison logic, VCF, insert/BLAST parsing, plotting and summary tables using synthetic data. `python3 -m unittest discover -s tests -v` additionally checks NANOPlasty output selection with a stand-in Perl assembler and Flye tandem normalization on synthetic one-, two- and three-copy contigs. **No real NANOPlasty assembly or full server Snakemake run was executed here**. The package does not contain raw patient/experiment reads or the upstream NANOPlasty source.

Source: https://github.com/ndierckx/NANOPlasty (v0.4 README, configuration instructions).


## Main validation question: insertion and repeat retention

**Did each sequenced plasmid retain the intended insert and the correct TALEN repeat architecture after cloning?**

The three complete designed FASTA sequences are sufficient for an *end-to-end
observed-versus-designed construct comparison*: if the entire experimental
assembly matches its complete design, the sequence is consistent with retention
of all designed regions, including the intended insert and repetitive arrays.
An exact design match is **not** proof that the insert arose via the intended
cloning event. With no pre-insertion backbone, insert FASTA, or known insert
coordinates, the program cannot uniquely identify *which portion was inserted*
or independently verify insertion-specific junctions. The insert-specific table
is therefore marked `ANNOTATIONS_MISSING` until this information is provided,
rather than reporting a false confirmation.

### Repeat integrity screen (automatic; no annotation required)

Every sample now produces `results/<sample>/repeat_integrity/summary.tsv`.
The screen finds candidate repetitive design intervals using recurrent 30-mers
(at least three occurrences), anchors each candidate to unique exact flanks in
the design-oriented assembly, and compares the extracted assembly sequence
with the designed region. It reports region length, substitutions, gapped
bases, representative 30-mer counts, and primary reads aligned across both
flanks (mapping quality at least 20). A minimum of three spanning alignments
is used as a *review threshold*; this threshold is not clinical or biochemical
proof of sequence correctness.

These intervals are computational repeat candidates, **not verified TALEN
feature annotations**. Similar TALE modules are often around 102 bp; this
screen does not assign every 102-bp module, independently confirm the RVD
motif order, or eliminate mapping ambiguity. A `SEQUENCE_MATCH` result
indicates that the selected assembly retained the corresponding design
subsequence at these anchors. A mismatch, unresolved flank, no discovered
candidate, or insufficient spanning evidence is flagged for manual review.
Use raw read alignments and independent Flye and NANOPlasty assemblies before
claiming repeat architecture is correct.

The integrated `results/summary/all_samples.tsv` now contains
`insert_validation_status`, `repeat_candidate_regions`,
`repeat_candidate_regions_matching`, `repeat_min_spanning_alignments`,
`repeat_integrity_status`, and `needs_review`.

Run `python3 -m unittest discover -s tests -v` to check the software.
Synthetic tests do not replace real-data validation across the nine barcodes.
