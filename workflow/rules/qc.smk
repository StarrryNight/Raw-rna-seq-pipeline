  rule fastqc_pretrim:                              
      input:                                        
          fq = lambda wc:                           
  config["samples"][wc.sample][wc.read]             
      output:
          html = f"{RESULTS}/{{sample}}/qc/pretrim/{
  {read}}_fastqc.html",
          zip = f"{RESULTS}/{{sample}}/qc/pretrim/{{
  read}}_fastqc.zip"
      params:
          outdir = lambda wc:
  f"{RESULTS}/{wc.sample}/qc/pretrim"
      threads: config["fastqc"]["threads"]
      resources:
          mem_mb = 16000,
          time = "2:00:00"
      envmodules:
          config["modules"]["fastqc"]
      shell:
          "fastqc -o {params.outdir} -t {threads}
  {input.fq}"