rule record_stats:
    input:
        log    = f"{PROCESS}/{{sample}}/yac_alignment/Log.final.out",
        bam    = f"{PROCESS}/{{sample}}/yac_alignment/Aligned.sortedByCoord.out.bam",
        bai    = f"{PROCESS}/{{sample}}/yac_alignment/Aligned.sortedByCoord.out.bam.bai",
        fwd_bw = f"{PROCESS}/{{sample}}/yac_alignment/{{sample}}_forward.bw",
        rev_bw = f"{PROCESS}/{{sample}}/yac_alignment/{{sample}}_reverse.bw",
        bed    = lambda wc: checkpoints.find_yac_region.get(sample=wc.sample).output.bed
    output:
        stats = f"{RESULTS}/{{sample}}/final_result/stats.json"
    envmodules:
        config["modules"]["python"],
        config["modules"]["samtools"]
    resources:
        mem_mb = 2000,
        time   = "00:05:00"
    shell:
        '''
        source {config[python_venv]}/bin/activate
        python workflow/scripts/record_stats.py \
            --sample  {wildcards.sample} \
            --log     {input.log} \
            --bam     {input.bam} \
            --bed     {input.bed} \
            --fwd_bw  {input.fwd_bw} \
            --rev_bw  {input.rev_bw} \
            --output  {output.stats}
        '''
