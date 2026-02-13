configfile: "config/config.yaml"

SAMPLES = list(configp"samples".keys())
RESULTS = config["results_dir"]

def get_yac_bed_values(wildcards):
      bed_file = checkpoints.find_yac_region.get(sam
  ple=wildcards.sample).output.bed
      with open(bed_file) as f:
          parts = f.read().strip().split("\t")
          return {"chrom": parts[0], "start":
  parts[1], "end": parts[2]}

  include: "workflow/rules/qc.smk"
  include: "workflow/rules/genome_prep.smk"
  include: "workflow/rules/alignment_pass1.smk"
  include: "workflow/rules/yac_discovery.smk"
  include: "workflow/rules/alignment_pass2.smk"
  include: "workflow/rules/coverage.smk"
  include: "workflow/rules/final_outputs.smk"

 rule all:
      input:
          expand(f"{RESULTS}/{{sample}}/final/{{samp
  le}}_fwd.npy", sample=SAMPLES),
          expand(f"{RESULTS}/{{sample}}/final/{{samp
  le}}_rev.npy", sample=SAMPLES),
          expand(f"{RESULTS}/{{sample}}/final/human_
  yac_insert.fa", sample=SAMPLES),
          expand(f"{RESULTS}/{{sample}}/forward.bw", sample=SAMPLES),
          expand(f"{RESULTS}/{{sample}}/reverse.bw", sample=SAMPLES),

