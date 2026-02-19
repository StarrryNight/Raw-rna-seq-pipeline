import os

configfile: "config/config.yaml"

SAMPLES = list(config["samples"].keys())
PROCESS = config["process_dir"]
RESULTS = config["results_dir"]

def get_yac_bed_values(wildcards):
    bed_file = checkpoints.find_yac_region.get(sample=wildcards.sample).output.bed
    with open(bed_file) as f:
        parts = f.read().strip().split("\t")
        return {"chrom": parts[0], "start": parts[1], "end": parts[2]}

include: "workflow/rules/qc.smk"
include: "workflow/rules/trimmomatic.smk"
include: "workflow/rules/genome_generation1.smk"
include: "workflow/rules/concatenate1.smk"
include: "workflow/rules/alignment1.smk"
include: "workflow/rules/find_human_yac.smk"
include: "workflow/rules/create_bw.smk"
include: "workflow/rules/concatenate2.smk"
include: "workflow/rules/genome_generation2.smk"
include: "workflow/rules/alignment2.smk"
include: "workflow/rules/plot_yac.smk"

rule all:
    input:
        expand(f"{RESULTS}/{{sample}}/trimmed/R1_paired.fq.gz", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/trimmed/R2_paired.fq.gz", sample=SAMPLES),
        f"{RESULTS}/full_hybrid_fasta/hybrid_fasta.fa",
        f"{RESULTS}/full_hybrid_genome_index",
        expand(f"{RESULTS}/{{sample}}/aligned/Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/aligned/Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/yac/yac_region.bed", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/yac_alignment/{{sample}}_forward.bw", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/yac_alignment/{{sample}}_reverse.bw", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/yac/hybrid_fasta.fa", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/yac/yac_genome_index", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/yac_alignment/Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/yac_alignment/Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/yac_alignment/{{sample}}_coverage.png", sample=SAMPLES),
