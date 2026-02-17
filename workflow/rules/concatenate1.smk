rule concat1:
    input:
        human_fa = config["ref"]["human_fasta"],
        yeast_fa = config["ref"]["yeast_fasta"]
    output:
        hybrid_fa = f"{RESULTS}/full_hybrid_fasta/hybrid_fasta.fa"
    threads: 1
    resources:
        mem_mb = 16000,
        time = "1:00:00"
    shell:
        '''
        sed 's/>/>sacCer3_/g' {input.yeast_fa} | cat {input.human_fa} - > {output.hybrid_fa}
        '''