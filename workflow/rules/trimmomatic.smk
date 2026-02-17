rule trimmomatic:
    input:
        r1 = lambda wc: config["samples"][wc.sample]["R1"],
        r2 = lambda wc: config["samples"][wc.sample]["R2"]
    output:
        r1_paired = f"{RESULTS}/{{sample}}/trimmed/R1_paired.fq.gz",
        r1_unpaired = f"{RESULTS}/{{sample}}/trimmed/R1_unpaired.fq.gz",
        r2_paired = f"{RESULTS}/{{sample}}/trimmed/R2_paired.fq.gz",
        r2_unpaired = f"{RESULTS}/{{sample}}/trimmed/R2_unpaired.fq.gz"
    threads: config["trimmomatic"]["threads"]
    resources:
        mem_mb = 16000,
        time = "4:00:00"
    envmodules:
        config["modules"]["trimmomatic"]
    shell:
        '''
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
        -threads {threads} \
        -phred33 \
        {input.r1} {input.r2} \
        {output.r1_paired} {output.r1_unpaired} \
        {output.r2_paired} {output.r2_unpaired} \
        ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        '''