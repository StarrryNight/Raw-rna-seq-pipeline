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

include: "workflow/rules/0_qc.smk"
include: "workflow/rules/1_preprocess.smk"
include: "workflow/rules/2_trimmomatic.smk"
include: "workflow/rules/3_genome_generation1.smk"
include: "workflow/rules/4_concatenate1.smk"
include: "workflow/rules/5_alignment1.smk"
include: "workflow/rules/6_find_human_yac.smk"
include: "workflow/rules/7_create_bw.smk"
include: "workflow/rules/8_concatenate2.smk"
include: "workflow/rules/9_genome_generation2.smk"
include: "workflow/rules/10_alignment2.smk"
include: "workflow/rules/11_plot_yac.smk"
include: "workflow/rules/12_stats.smk"

rule all:
    input:
        expand(f"{PROCESS}/{{sample}}/trimmed/R1_paired.fq.gz", sample=SAMPLES),
        expand(f"{PROCESS}/{{sample}}/trimmed/R2_paired.fq.gz", sample=SAMPLES),
        "reference_fastas/hybrid_fasta.fa",
        "reference_fastas/full_hybrid_genome_index",
        expand(f"{PROCESS}/{{sample}}/aligned/Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES),
        expand(f"{PROCESS}/{{sample}}/aligned/Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        expand(f"{PROCESS}/{{sample}}/yac/yac_region.bed", sample=SAMPLES),
        expand(f"{PROCESS}/{{sample}}/yac_alignment/{{sample}}_forward.bw", sample=SAMPLES),
        expand(f"{PROCESS}/{{sample}}/yac_alignment/{{sample}}_reverse.bw", sample=SAMPLES),
        expand(f"{PROCESS}/{{sample}}/yac/hybrid_fasta.fa", sample=SAMPLES),
        expand(f"{PROCESS}/{{sample}}/yac/yac_genome_index", sample=SAMPLES),
        expand(f"{PROCESS}/{{sample}}/yac_alignment/Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES),
        expand(f"{PROCESS}/{{sample}}/yac_alignment/Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/{{sample}}_fwd.npz", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/{{sample}}_rev.npz", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/{{sample}}_yac_coverage.png", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/{{sample}}_yeast_fwd.png", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/{{sample}}_yeast_rev.png", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/stats.json", sample=SAMPLES),
