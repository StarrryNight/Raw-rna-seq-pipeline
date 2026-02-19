rule create_bigwig:
    input:
        bam = f"{PROCESS}/{{sample}}/yac_alignment/Aligned.sortedByCoord.out.bam",
        bai = f"{PROCESS}/{{sample}}/yac_alignment/Aligned.sortedByCoord.out.bam.bai",
        bed = lambda wc: checkpoints.find_yac_region.get(sample=wc.sample).output.bed
    output:
        bw_forward = f"{PROCESS}/{{sample}}/yac_alignment/{{sample}}_forward.bw",
        bw_reverse = f"{PROCESS}/{{sample}}/yac_alignment/{{sample}}_reverse.bw"
    threads: config["bamcoverage"]["threads"]
    resources:
        mem_mb = 16000,
        time = "4:00:00"
    params:
        bin_size = config["bamcoverage"]["bin_size"],
        normalize_using = config["bamcoverage"]["normalize_using"]
    envmodules:
        config["modules"]["python"]
    shell:
        '''
        bamCoverage -b {input.bam} -o {output.bw_reverse} \
                    --binSize {params.bin_size} \
                    --normalizeUsing {params.normalize_using} \
                    --filterRNAstrand reverse \
                    -p {threads}
        bamCoverage -b {input.bam} -o {output.bw_forward} \
                    --binSize {params.bin_size} \
                    --normalizeUsing {params.normalize_using} \
                    --filterRNAstrand forward \
                    -p {threads}
        '''
