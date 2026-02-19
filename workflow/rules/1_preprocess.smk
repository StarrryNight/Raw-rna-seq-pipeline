rule preprocess:
    input:
        r1 = lambda wc: config["samples"][wc.sample]["R1"],
        r2 = lambda wc: config["samples"][wc.sample]["R2"]
    output:
        r1 = f"{PROCESS}/{{sample}}/raw/R1_cat.fq",
        r2 = f"{PROCESS}/{{sample}}/raw/R2_cat.fq"
    resources:
        mem_mb = 4000,
        time = "00:15:00"
    shell:
        '''
        zcat {input.r1} > {output.r1}
        zcat {input.r2} > {output.r2}
        '''
