rule trimmomatic:
    input:
        r1 = f"{PROCESS}/{{sample}}/raw/R1_cat.fq",
        r2 = f"{PROCESS}/{{sample}}/raw/R2_cat.fq"
    output:
        r1_paired = f"{PROCESS}/{{sample}}/trimmed/R1_paired.fq.gz",
        r1_unpaired = f"{PROCESS}/{{sample}}/trimmed/R1_unpaired.fq.gz",
        r2_paired = f"{PROCESS}/{{sample}}/trimmed/R2_paired.fq.gz",
        r2_unpaired = f"{PROCESS}/{{sample}}/trimmed/R2_unpaired.fq.gz"
    threads: config["trimmomatic"]["threads"]
    resources:
        mem_mb = 16000,
        time = "00:30:00"
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