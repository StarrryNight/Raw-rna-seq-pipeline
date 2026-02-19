rule genome_generation2:
    input:
        hybrid_fa = f"{RESULTS}/{{sample}}/yac/hybrid_fasta.fa"
    output:
        hybrid_index = directory(f"{RESULTS}/{{sample}}/yac/yac_genome_index")
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