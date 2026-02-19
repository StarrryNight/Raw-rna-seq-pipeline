rule alignment:
    input:
        r1_paired = f"{PROCESS}/{{sample}}/trimmed/R1_paired.fq.gz",
        r2_paired = f"{PROCESS}/{{sample}}/trimmed/R2_paired.fq.gz",
        hybrid_index = "reference_fastas/full_hybrid_genome_index"
    output:
        bam = f"{PROCESS}/{{sample}}/aligned/Aligned.sortedByCoord.out.bam",
        bai = f"{PROCESS}/{{sample}}/aligned/Aligned.sortedByCoord.out.bam.bai"
    params:
        outprefix = lambda wc, output: os.path.dirname(output.bam) + "/",
        align_intron_max = config["star"]["align_intron_max"]
    threads: config["star"]["threads"]
    resources:
        mem_mb = 32000,
        time = "3:00:00"
    envmodules:
        config["modules"]["star"],
        config["modules"]["samtools"]
    shell:
        '''
        STAR --runMode alignReads \
             --runThreadN {threads} \
             --genomeDir {input.hybrid_index} \
             --readFilesIn {input.r1_paired} {input.r2_paired} \
             --readFilesCommand zcat \
             --alignIntronMax {params.align_intron_max} \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix {params.outprefix}

        samtools index {output.bam}
        '''
