rule concat2:
    input:
        human_fa = config["ref"]["human_fasta"],
        yeast_fa = config["ref"]["yeast_fasta"],
        bed_file = lambda wc: checkpoints.find_yac_region.get(sample=wc.sample).output.bed
    output:
        hybrid_fa = f"{PROCESS}/{{sample}}/yac/hybrid_fasta.fa"
    threads: 1
    envmodules:
        config["modules"]["samtools"]
    resources:
        mem_mb = 16000,
        time = "1:00:00"
    shell:
        '''
        TOP_CHR=$(awk '{{print $1}}' {input.bed_file})
        START=$(awk '{{print $2}}' {input.bed_file})
        END=$(awk '{{print $3}}' {input.bed_file})

        samtools faidx {input.human_fa} "$TOP_CHR:$START-$END" | \
        sed "s/^>.*/>$TOP_CHR/" > {output.hybrid_fa}.tmp

        cat {input.yeast_fa} {output.hybrid_fa}.tmp > {output.hybrid_fa}
        rm {output.hybrid_fa}.tmp

        samtools faidx {output.hybrid_fa}
        '''
