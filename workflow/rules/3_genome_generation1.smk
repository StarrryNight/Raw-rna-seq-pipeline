rule genome_generation1:
    input:
        hybrid_fa = "reference_fastas/hybrid_fasta.fa"
    output:
        hybrid_index = directory(f"reference_fastas/full_hybrid_genome_index")
    threads: config["star"]["threads"]
    resources:
        mem_mb = 50000,
        time = "8:00:00"
    envmodules:
        config["modules"]["star"]
    shell:
        '''
        mkdir -p {output.hybrid_index}
        STAR --runMode genomeGenerate \
             --runThreadN {threads} \
             --genomeDir {output.hybrid_index} \
             --genomeFastaFiles {input.hybrid_fa} \
             --genomeSAindexNbases {config[star][genome_sa_index_nbases]}
        '''