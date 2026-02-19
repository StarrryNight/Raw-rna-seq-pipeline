rule yac_numpy:
    input:
        fwd_bw = f"{PROCESS}/{{sample}}/yac_alignment/{{sample}}_forward.bw",
        rev_bw = f"{PROCESS}/{{sample}}/yac_alignment/{{sample}}_reverse.bw",
        bed = lambda wc: checkpoints.find_yac_region.get(sample=wc.sample).output.bed
    output:
        fwd_npy = f"{PROCESS}/{{sample}}/final_result/yac_sample_1_fwd.npy",
        rev_npy = f"{PROCESS}/{{sample}}/final_result/yac_sample_1_rev.npy"
    params:
        outdir = lambda wc, output: os.path.dirname(output.fwd_npy),
        bin_size = config["bigwig_to_numpy"]["bin_size"]
    envmodules:
        config["modules"]["python"]
    resources:
        mem_mb = 4000,
        time = "00:30:00"
    shell:
        '''
        source {config[python_venv]}/bin/activate
        CHR=$(awk '{{print $1}}' {input.bed})
        LEN=$(awk '{{print $3-$2}}' {input.bed})
        python workflow/scripts/bw_to_np.py \
            --chrom $CHR \
            --start 0 \
            --end $LEN \
            --fwd_bw {input.fwd_bw} \
            --rev_bw {input.rev_bw} \
            --output_dir {params.outdir}
        '''


rule yac_plot:
    input:
        fwd_npy = f"{RESULTS}/{{sample}}/final_result/yac_sample_1_fwd.npy",
        rev_npy = f"{RESULTS}/{{sample}}/final_result/yac_sample_1_rev.npy",
        bed = lambda wc: checkpoints.find_yac_region.get(sample=wc.sample).output.bed
    output:
        plot = f"{RESULTS}/{{sample}}/final_result/{{sample}}_coverage.png"
    envmodules:
        config["modules"]["python"]
    resources:
        mem_mb = 4000,
        time = "00:15:00"
    shell:
        '''
        source {config[python_venv]}/bin/activate
        python workflow/scripts/plot_yac.py \
            --fwd_npy {input.fwd_npy} \
            --rev_npy {input.rev_npy} \
            --bed {input.bed} \
            --output {output.plot}
        '''
