rule yac_numpy:
    input:
        fwd_bw = f"{PROCESS}/{{sample}}/yac_alignment/{{sample}}_forward.bw",
        rev_bw = f"{PROCESS}/{{sample}}/yac_alignment/{{sample}}_reverse.bw",
        bed = lambda wc: checkpoints.find_yac_region.get(sample=wc.sample).output.bed
    output:
        npz = f"{RESULTS}/{{sample}}/final_result/coverage.npz"
    envmodules:
        config["modules"]["python"]
    resources:
        mem_mb = 4000,
        time = "00:30:00"
    shell:
        '''
        source {config[python_venv]}/bin/activate
        python workflow/scripts/bw_to_np.py \
            --fwd_bw {input.fwd_bw} \
            --rev_bw {input.rev_bw} \
            --bed {input.bed} \
            --output {output.npz}
        '''


rule yac_plot:
    input:
        npz = f"{RESULTS}/{{sample}}/final_result/coverage.npz",
        bed = lambda wc: checkpoints.find_yac_region.get(sample=wc.sample).output.bed
    output:
        yac_plot  = f"{RESULTS}/{{sample}}/final_result/{{sample}}_yac_coverage.png",
        yeast_fwd = f"{RESULTS}/{{sample}}/final_result/{{sample}}_yeast_fwd.png",
        yeast_rev = f"{RESULTS}/{{sample}}/final_result/{{sample}}_yeast_rev.png"
    envmodules:
        config["modules"]["python"]
    resources:
        mem_mb = 4000,
        time = "00:15:00"
    shell:
        '''
        source {config[python_venv]}/bin/activate
        python workflow/scripts/plot_yac.py \
            --npz {input.npz} \
            --bed {input.bed} \
            --out_yac {output.yac_plot} \
            --out_yeast_fwd {output.yeast_fwd} \
            --out_yeast_rev {output.yeast_rev}
        '''
