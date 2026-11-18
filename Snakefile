"""Plasmid pipeline: nanopore reads -> assembly -> comparison with the designed plasmid.

Every sample gets its own folder:  results/<sample>/
Overall tables and figure:         results/summary/

Steps, in order (each is a "rule" below):
   1. filter_reads                 keep good reads, write read statistics
   2. write_nanoplasty_config       write NANOPlasty settings file
   3. assemble_with_nanoplasty      de novo circular assembly (seed-assisted)
   4. orient_assembly              flip/rotate the assembly to match the design
   5. map_reads / map_assembly     minimap2 alignments (BAM files for IGV)
   6. compare_assembly_to_design   base-by-base differences, VCF, read support
   7. blast_assembly_vs_design     BLAST of the assembly against the design
   8. check_inserts                are your inserts in the assembly?
   9. run_sniffles                 structural-variant candidates from the raw reads
  10. (optional) assemble_with_flye / compare_flye   independent cross-check
  11. plot_sample / make_summary_table / plot_overview
"""
import csv
import glob
import os
import shlex
from pathlib import Path

configfile: "config/config.yaml"

BASE = Path(workflow.basedir).resolve()


# --------------------------------------------------------------------------------------
# Settings and sample sheet
# --------------------------------------------------------------------------------------
with open(BASE / "samples.csv", newline="") as handle:
    SAMPLE_SHEET = {row["sample"]: row for row in csv.DictReader(handle)}
SAMPLES = list(SAMPLE_SHEET)
if not SAMPLES:
    raise ValueError("samples.csv has no rows")

THREADS = int(config.get("threads", 4))
RUN_SNIFFLES = bool(config.get("run_sniffles", True))
RUN_FLYE = bool(config.get("run_flye", False))
INSERT_SEQUENCES = str(config.get("insert_sequences_fasta", "") or "")
SNIFFLES_MIN_SUPPORT = str(config.get("sniffles_min_support", "auto")).lower()

wildcard_constraints:
    sample="|".join(SAMPLES),
    version="design|rotated_design"      # which copy of the design the reads are aligned to


def design_fasta(w):
    return str(BASE / SAMPLE_SHEET[w.sample]["design_fasta"])

def seed_fasta(w):
    return str(BASE / SAMPLE_SHEET[w.sample]["seed_fasta"])

def design_length(w):
    return int(SAMPLE_SHEET[w.sample]["expected_length_bp"])

def construct_name(w):
    return SAMPLE_SHEET[w.sample]["construct"]

def reads_pattern(w):
    """File name or glob pattern of this sample's FASTQ files, e.g. .../barcode01/*.fastq.gz"""
    folder = Path(os.path.expanduser(str(config["reads_folder"]))).resolve()
    return str(folder / config["reads_filename"].format(
        barcode=SAMPLE_SHEET[w.sample]["barcode"], sample=w.sample))

def reads_files(w):
    """Every FASTQ file of this sample (several files per barcode are joined automatically)."""
    files = sorted(f for f in glob.glob(reads_pattern(w)) if os.path.isfile(f))
    if not files:
        raise ValueError(f"No FASTQ files found for {w.sample}. Pattern searched: {reads_pattern(w)}\n"
                         "Check reads_folder and reads_filename in config/config.yaml, "
                         "and the barcode name in samples.csv.")
    return files

def reference_for_version(w):
    """The design itself, or its half-rotated copy."""
    return design_fasta(w) if w.version == "design" else f"work/{w.sample}/rotated_design.fasta"

def sniffles_min_support_option():
    return "" if SNIFFLES_MIN_SUPPORT == "auto" else f"--minsupport {int(SNIFFLES_MIN_SUPPORT)}"

def insert_sequences_input(w):
    return [str(BASE / INSERT_SEQUENCES)] if INSERT_SEQUENCES else []

def insert_sequences_option(w):
    return f"--insert-sequences {shlex.quote(str(BASE / INSERT_SEQUENCES))}" if INSERT_SEQUENCES else ""


# Files per sample that the summary table needs (names relative to results/<sample>/)
SUMMARY_INPUTS = [
    "reads/read_stats.tsv",
    "nanoplasty/assembly_source.tsv",
    "assembly/assembly_oriented.fasta",
    "assembly/orientation.tsv",
    "assembly_vs_design/summary.tsv",
    "assembly_vs_design/differences.vcf",
    "assembly_vs_design/reads_across_origin.tsv",
    "alignments/read_coverage.tsv",
    "blast/whole_plasmid_summary.tsv",
    "blast/inserts_summary.tsv",
    "repeat_integrity/summary.tsv",
]
if RUN_SNIFFLES:
    SUMMARY_INPUTS.append("sniffles/summary.tsv")
if RUN_FLYE:
    SUMMARY_INPUTS += ["flye/normalization.tsv", "flye/orientation.tsv", "flye/vs_design.tsv", "flye/vs_nanoplasty.tsv"]


rule all:
    input:
        "results/summary/all_samples.tsv",
        "results/summary/replicates.tsv",
        "results/summary/overview.png",
        expand("results/{sample}/plots/{plot}.png", sample=SAMPLES,
               plot=["read_lengths", "depth", "differences"]),
        expand("results/{sample}/alignments/assembly_to_design.bam.bai", sample=SAMPLES),
        expand("results/{sample}/assembly_vs_design/read_support.tsv", sample=SAMPLES),


# ======================================================================================
# 1-3. Reads and de novo assembly
# ======================================================================================
rule filter_reads:
    input:
        fastq=reads_files
    output:
        fastq="work/{sample}/filtered_reads/reads.fastq",
        stats="results/{sample}/reads/read_stats.tsv",
        histogram="results/{sample}/reads/read_length_histogram.tsv"
    params:
        pattern=reads_pattern,
        barcode=lambda w: SAMPLE_SHEET[w.sample]["barcode"],
        design_length=design_length,
        min_length=int(config.get("min_read_length", 1000)),
        min_quality=float(config.get("min_mean_read_quality", 0))
    shell:
        "python3 scripts/filter_reads.py --fastq-pattern {params.pattern:q} "
        "--output-fastq {output.fastq:q} --stats {output.stats:q} "
        "--length-histogram {output.histogram:q} --sample {wildcards.sample:q} "
        "--barcode {params.barcode:q} --design-length {params.design_length} "
        "--min-length {params.min_length} --min-quality {params.min_quality}"


rule write_nanoplasty_config:
    input:
        fastq="work/{sample}/filtered_reads/reads.fastq",
        seed=seed_fasta,
        design=design_fasta
    output:
        "work/{sample}/nanoplasty_config.txt"
    params:
        design_length=design_length,
        min_length=int(config.get("nanoplasty_min_read_length", 3000)),
        range_fraction=float(config.get("nanoplasty_genome_range_fraction", 0.15)),
        use_full="yes" if config.get("nanoplasty_use_full_dataset", False) else "no",
        extend_seed="yes" if config.get("nanoplasty_extend_seed_directly", False) else "no",
        save_reads="yes" if config.get("nanoplasty_save_assembled_reads", False) else "no",
        # All paths are absolute for NANOPlasty, regardless of the run directory.
        reads_folder=lambda w: str((BASE / "work" / w.sample / "filtered_reads").resolve()),
        output_folder=lambda w: str((BASE / "results" / w.sample / "nanoplasty" / "native").resolve())
    threads: THREADS
    shell:
        "python3 scripts/write_nanoplasty_config.py --sample {wildcards.sample:q} "
        "--seed {input.seed:q} --design-length {params.design_length} "
        "--reads-folder {params.reads_folder:q} --output-folder {params.output_folder:q} "
        "--config-file {output:q} --threads {threads} "
        "--min-read-length {params.min_length} --range-fraction {params.range_fraction} "
        "--use-full-dataset {params.use_full:q} --extend-seed-directly {params.extend_seed:q} "
        "--save-assembled-reads {params.save_reads:q}"


rule assemble_with_nanoplasty:
    input:
        config_file="work/{sample}/nanoplasty_config.txt",
        fastq="work/{sample}/filtered_reads/reads.fastq"
    output:
        assembly="results/{sample}/nanoplasty/assembly_raw_{sample}.fasta",
        source="results/{sample}/nanoplasty/assembly_source.tsv"
    params:
        script=str(config.get("nanoplasty_script", "tools/NANOPlasty/NANOPlasty0.4.pl")),
        native=lambda w: str((BASE / "results" / w.sample / "nanoplasty" / "native").resolve()),
        design_length=design_length
    threads: THREADS
    log:
        "results/{sample}/logs/nanoplasty.log"
    shell:
        "python3 scripts/run_nanoplasty.py --script {params.script:q} "
        "--config-file {input.config_file:q} --native-folder {params.native:q} "
        "--sample {wildcards.sample:q} --design-length {params.design_length} "
        "--output-fasta {output.assembly:q} --source-report {output.source:q} "
        "> {log:q} 2>&1"


# ======================================================================================
# 4. Make the assembly comparable with the design
# ======================================================================================
rule orient_assembly:
    input:
        assembly="results/{sample}/nanoplasty/assembly_raw_{sample}.fasta",
        design=design_fasta
    output:
        fasta="results/{sample}/assembly/assembly_oriented.fasta",
        report="results/{sample}/assembly/orientation.tsv"
    shell:
        "python3 scripts/orient_assembly.py --assembly {input.assembly:q} "
        "--design {input.design:q} --name {wildcards.sample:q} "
        "--output-fasta {output.fasta:q} --output-report {output.report:q}"


rule rotate_design:
    # A copy of the design rotated by half, so its artificial start/end sits in the middle.
    input:
        design=design_fasta
    output:
        fasta="work/{sample}/rotated_design.fasta",
        origin_file="work/{sample}/rotated_design_origin.txt"
    shell:
        "python3 scripts/rotate_design.py --design {input.design:q} "
        "--output-fasta {output.fasta:q} --origin-file {output.origin_file:q}"


rule copy_and_index_design:
    # bcftools needs an indexed FASTA. Each sample indexes its own copy, so this works for any
    # design file name or location and parallel samples never write the same index file.
    input:
        design=design_fasta
    output:
        fasta="work/{sample}/design.fasta",
        index="work/{sample}/design.fasta.fai"
    shell:
        "cp {input.design:q} {output.fasta:q}; "
        "samtools faidx {output.fasta:q}"


# ======================================================================================
# 5. Alignments (BAM files you can open in IGV)
# ======================================================================================
rule map_reads:
    # Raw reads against the design, and against the rotated design.
    input:
        fastq="work/{sample}/filtered_reads/reads.fastq",
        reference=reference_for_version
    output:
        bam="results/{sample}/alignments/reads_to_{version}.bam",
        bai="results/{sample}/alignments/reads_to_{version}.bam.bai"
    threads: THREADS
    shell:
        "minimap2 -t {threads} -ax map-ont {input.reference:q} {input.fastq:q} "
        "| samtools sort -@ {threads} -o {output.bam:q} -; "
        "samtools index {output.bam:q} {output.bai:q}"


rule map_assembly_to_design:
    input:
        assembly="results/{sample}/assembly/assembly_oriented.fasta",
        design=design_fasta
    output:
        bam="results/{sample}/alignments/assembly_to_design.bam",
        bai="results/{sample}/alignments/assembly_to_design.bam.bai"
    threads: THREADS
    shell:
        "minimap2 -t {threads} -ax asm5 --cs {input.design:q} {input.assembly:q} "
        "| samtools sort -@ {threads} -o {output.bam:q} -; "
        "samtools index {output.bam:q} {output.bai:q}"


rule read_alignment_stats:
    input:
        bam="results/{sample}/alignments/reads_to_design.bam",
        bai="results/{sample}/alignments/reads_to_design.bam.bai"
    output:
        flagstat="results/{sample}/alignments/read_alignment_flagstat.txt",
        coverage="results/{sample}/alignments/read_coverage.tsv",
        depth="results/{sample}/alignments/read_depth.tsv"
    threads: THREADS
    shell:
        "samtools flagstat -@ {threads} {input.bam:q} > {output.flagstat:q}; "
        "samtools coverage {input.bam:q} > {output.coverage:q}; "
        "samtools depth -aa {input.bam:q} > {output.depth:q}"


# ======================================================================================
# 6. Assembly versus design
# ======================================================================================
rule compare_assembly_to_design:
    input:
        design=design_fasta,
        assembly="results/{sample}/assembly/assembly_oriented.fasta"
    output:
        summary="results/{sample}/assembly_vs_design/summary.tsv",
        differences="results/{sample}/assembly_vs_design/differences.tsv"
    shell:
        "python3 scripts/compare_sequences.py --reference {input.design:q} "
        "--other {input.assembly:q} --reference-label design --other-label assembly "
        "--summary {output.summary:q} --differences {output.differences:q}"


rule write_vcf:
    # Differences -> VCF, then left-aligned with bcftools norm.
    input:
        design="work/{sample}/design.fasta",
        design_index="work/{sample}/design.fasta.fai",
        differences="results/{sample}/assembly_vs_design/differences.tsv"
    output:
        "results/{sample}/assembly_vs_design/differences.vcf"
    shell:
        "python3 scripts/differences_to_vcf.py --design {input.design:q} "
        "--differences {input.differences:q} --output {output:q}.unnormalized; "
        "bcftools norm -f {input.design:q} -c e -O v -o {output:q} {output:q}.unnormalized; "
        "rm {output:q}.unnormalized"


rule check_read_support:
    input:
        bam="results/{sample}/alignments/reads_to_design.bam",
        bai="results/{sample}/alignments/reads_to_design.bam.bai",
        differences="results/{sample}/assembly_vs_design/differences.tsv"
    output:
        "results/{sample}/assembly_vs_design/read_support.tsv"
    shell:
        "python3 scripts/check_read_support.py --bam {input.bam:q} "
        "--differences {input.differences:q} --output {output:q}"


rule check_reads_across_origin:
    input:
        bam="results/{sample}/alignments/reads_to_rotated_design.bam",
        bai="results/{sample}/alignments/reads_to_rotated_design.bam.bai",
        origin_file="work/{sample}/rotated_design_origin.txt"
    output:
        "results/{sample}/assembly_vs_design/reads_across_origin.tsv"
    shell:
        "python3 scripts/check_reads_across_origin.py --bam {input.bam:q} "
        "--origin-file {input.origin_file:q} --output {output:q}"


# ======================================================================================
# 7-8. BLAST: whole plasmid, and your inserts
# ======================================================================================
rule blast_assembly_vs_design:
    input:
        assembly="results/{sample}/assembly/assembly_oriented.fasta",
        design=design_fasta
    output:
        hits="results/{sample}/blast/whole_plasmid.tsv",
        summary="results/{sample}/blast/whole_plasmid_summary.tsv"
    shell:
        "blastn -task blastn -dust no -query {input.assembly:q} -subject {input.design:q} "
        "-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send "
        "evalue bitscore qlen slen' > {output.hits:q}; "
        "python3 scripts/summarize_blast.py --blast {output.hits:q} --output {output.summary:q}"



rule check_repeat_regions:
    # Candidate repeat-rich regions inferred from the design, not known annotations.
    input:
        design=design_fasta,
        assembly="results/{sample}/assembly/assembly_oriented.fasta",
        bam="results/{sample}/alignments/reads_to_design.bam",
        bai="results/{sample}/alignments/reads_to_design.bam.bai"
    output:
        "results/{sample}/repeat_integrity/summary.tsv"
    shell:
        "python3 scripts/check_repeat_regions.py --design {input.design:q} "
        "--assembly {input.assembly:q} --bam {input.bam:q} "
        "--sample {wildcards.sample:q} --output {output:q}"


rule check_inserts:
    input:
        design=design_fasta,
        assembly="results/{sample}/assembly/assembly_oriented.fasta",
        coordinates=str(BASE / "config" / "insert_coordinates.csv"),
        sequences=insert_sequences_input
    output:
        summary="results/{sample}/blast/inserts_summary.tsv",
        queries="results/{sample}/blast/insert_sequences.fasta"
    params:
        construct=construct_name,
        sequences_option=insert_sequences_option
    shell:
        "python3 scripts/check_inserts.py --design {input.design:q} "
        "--assembly {input.assembly:q} --coordinates {input.coordinates:q} "
        "--construct {params.construct:q} --output {output.summary:q} "
        "--queries-out {output.queries:q} {params.sequences_option}"


# ======================================================================================
# 9. Sniffles2: structural-variant candidates from the raw reads
#    Run on the design AND on the rotated design, so the artificial start/end
#    of the circular FASTA is not a blind spot.
# ======================================================================================
rule run_sniffles:
    input:
        bam="results/{sample}/alignments/reads_to_{version}.bam",
        bai="results/{sample}/alignments/reads_to_{version}.bam.bai",
        reference=reference_for_version
    output:
        "results/{sample}/sniffles/reads_to_{version}.vcf"
    threads: THREADS
    params:
        min_sv_length=int(config.get("sniffles_min_sv_length", 50)),
        min_support=sniffles_min_support_option()
    log:
        "results/{sample}/logs/sniffles_{version}.log"
    shell:
        "sniffles --input {input.bam:q} --vcf {output:q} --reference {input.reference:q} "
        "--threads {threads} --minsvlen {params.min_sv_length} {params.min_support} "
        "> {log:q} 2>&1"


rule summarize_sniffles:
    input:
        design_vcf="results/{sample}/sniffles/reads_to_design.vcf",
        rotated_vcf="results/{sample}/sniffles/reads_to_rotated_design.vcf"
    output:
        "results/{sample}/sniffles/summary.tsv"
    params:
        length=design_length
    shell:
        "python3 scripts/summarize_sniffles.py --design-vcf {input.design_vcf:q} "
        "--rotated-vcf {input.rotated_vcf:q} --length {params.length} --output {output:q}"


# ======================================================================================
# 10. Optional: Flye, an assembly that never sees the design or the seed
# ======================================================================================
rule assemble_with_flye:
    input:
        fastq="work/{sample}/filtered_reads/reads.fastq"
    output:
        "results/{sample}/flye/assembly_raw/assembly.fasta"
    threads: THREADS
    params:
        folder="results/{sample}/flye/assembly_raw",
        # --genome-size only lets Flye subsample very deep data; 0 turns subsampling off
        size_options=lambda w: (
            f"--genome-size {design_length(w)} --asm-coverage {int(config.get('flye_max_coverage', 100))}"
            if int(config.get("flye_max_coverage", 100)) > 0 else "")
    log:
        "results/{sample}/logs/flye.log"
    shell:
        "flye --nano-hq {input.fastq:q} --out-dir {params.folder:q} --threads {threads} "
        "{params.size_options} > {log:q} 2>&1; "
        "test -s {output:q}"


rule normalize_flye_assembly:
    # Preserve the original Flye assembly. A single observed copy is retained only
    # when a unique design anchor identifies near-length tandem units that agree.
    # This is sequence normalization, NOT proof that the biological plasmid is monomeric.
    input:
        raw="results/{sample}/flye/assembly_raw/assembly.fasta",
        graph="results/{sample}/flye/assembly_raw/assembly_graph.gfa",
        nanoplasty="results/{sample}/nanoplasty/assembly_raw_{sample}.fasta",
        design=design_fasta
    output:
        normalized="results/{sample}/flye/assembly_normalized.fasta",
        report="results/{sample}/flye/normalization.tsv",
        boundaries="results/{sample}/flye/copy_boundaries.tsv",
        pairwise="results/{sample}/flye/copy_pairwise.tsv",
        copies="results/{sample}/flye/copies.fasta"
    log:
        "results/{sample}/logs/flye_normalization.log"
    shell:
        "python3 scripts/normalize_flye_tandem.py --assembly {input.raw:q} "
        "--graph {input.graph:q} --nanoplasty {input.nanoplasty:q} "
        "--design {input.design:q} --sample {wildcards.sample:q} "
        "--output-fasta {output.normalized:q} --report {output.report:q} "
        "--boundaries {output.boundaries:q} --pairwise {output.pairwise:q} "
        "--copies-fasta {output.copies:q} > {log:q} 2>&1"


rule compare_flye:
    # Orient the *normalized* Flye unit, then compare it against the design
    # and independently generated NANOPlasty assembly. Keep Flye's raw contig.
    input:
        flye="results/{sample}/flye/assembly_normalized.fasta",
        nanoplasty="results/{sample}/assembly/assembly_oriented.fasta",
        design=design_fasta
    output:
        oriented="results/{sample}/flye/assembly_oriented.fasta",
        orientation="results/{sample}/flye/orientation.tsv",
        vs_design="results/{sample}/flye/vs_design.tsv",
        vs_design_differences="results/{sample}/flye/vs_design_differences.tsv",
        vs_nanoplasty="results/{sample}/flye/vs_nanoplasty.tsv",
        vs_nanoplasty_differences="results/{sample}/flye/vs_nanoplasty_differences.tsv"
    shell:
        "python3 scripts/orient_assembly.py --assembly {input.flye:q} --design {input.design:q} "
        "--name {wildcards.sample:q} --output-fasta {output.oriented:q} "
        "--output-report {output.orientation:q}; "
        "python3 scripts/compare_sequences.py --reference {input.design:q} "
        "--other {output.oriented:q} --reference-label design --other-label flye "
        "--summary {output.vs_design:q} --differences {output.vs_design_differences:q}; "
        "python3 scripts/compare_sequences.py --reference {input.nanoplasty:q} "
        "--other {output.oriented:q} --reference-label nanoplasty --other-label flye "
        "--summary {output.vs_nanoplasty:q} --differences {output.vs_nanoplasty_differences:q}"


# ======================================================================================
# 11. Figures and summary tables
# ======================================================================================
rule plot_sample:
    input:
        read_stats="results/{sample}/reads/read_stats.tsv",
        histogram="results/{sample}/reads/read_length_histogram.tsv",
        depth="results/{sample}/alignments/read_depth.tsv",
        differences="results/{sample}/assembly_vs_design/differences.tsv"
    output:
        expand("results/{{sample}}/plots/{plot}.{extension}",
               plot=["read_lengths", "depth", "differences"], extension=["png", "pdf"])
    shell:
        "python3 scripts/plot_sample.py --sample {wildcards.sample:q} "
        "--read-stats {input.read_stats:q} --length-histogram {input.histogram:q} "
        "--depth {input.depth:q} --differences {input.differences:q} "
        "--output-folder results/{wildcards.sample}/plots"


rule make_summary_table:
    input:
        sample_sheet=str(BASE / "samples.csv"),
        per_sample=expand("results/{sample}/{file}", sample=SAMPLES, file=SUMMARY_INPUTS)
    output:
        summary="results/summary/all_samples.tsv",
        replicates="results/summary/replicates.tsv"
    params:
        sniffles="yes" if RUN_SNIFFLES else "no",
        flye="yes" if RUN_FLYE else "no"
    shell:
        "python3 scripts/make_summary_table.py --sample-sheet {input.sample_sheet:q} "
        "--results-folder results --summary-out {output.summary:q} "
        "--replicates-out {output.replicates:q} --sniffles {params.sniffles} --flye {params.flye}"


rule plot_overview:
    input:
        "results/summary/all_samples.tsv"
    output:
        png="results/summary/overview.png",
        pdf="results/summary/overview.pdf"
    shell:
        "python3 scripts/plot_overview.py --summary {input:q} --prefix results/summary/overview"
