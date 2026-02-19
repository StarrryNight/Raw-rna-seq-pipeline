rule fastqc_pretrim:
    input:
        fq = lambda wc: config["samples"][wc.sample][wc.read]
    output:
        html = f"{PROCESS}/{{sample}}/qc/pretrim/{{read}}_fastqc.html",
        zip = f"{PROCESS}/{{sample}}/qc/pretrim/{{read}}_fastqc.zip"
    params:
        outdir = lambda wc: f"{PROCESS}/{wc.sample}/qc/pretrim"
    threads: config["fastqc"]["threads"]
    resources:
        mem_mb = 16000,
        time = "2:00:00"
    envmodules:
        config["modules"]["fastqc"]
    shell:
        "fastqc -o {params.outdir} -t {threads} {input.fq}"
