checkpoint find_yac_region:
    input:
        bam = f"{PROCESS}/{{sample}}/aligned/Aligned.sortedByCoord.out.bam",
        bai = f"{PROCESS}/{{sample}}/aligned/Aligned.sortedByCoord.out.bam.bai"
    output:
        bed = f"{PROCESS}/{{sample}}/yac/yac_region.bed"
    params:
        trim_start_percentile = config["find_yac_region"]["trim_start_percentile"],
        trim_end_percentile = config["find_yac_region"]["trim_end_percentile"]
    resources:
        mem_mb = 4000,
        time = "00:15:00"
    envmodules:
        config["modules"]["samtools"]
    shell:
        '''
TOP_CHR=$(samtools idxstats {input.bam} | \
grep "^chr" | \
awk '$2 > 10000000' | \
sort -k3,3rn | \
head -n 1 | \
cut -f1)

echo "Top chromosome is: $TOP_CHR"

NUM_READS=$(samtools view -c {input.bam} "$TOP_CHR")
START_INDEX=$(( NUM_READS * {params.trim_start_percentile} / 100 ))
END_INDEX=$(( NUM_READS * {params.trim_end_percentile} / 100 ))

#TODO rename duplicated chrom names
read START END <<< $(samtools view {input.bam} "$TOP_CHR" 2>/dev/null | \
    awk '{{print $4}}' | \
    sort -n | \
    sed -n "${{START_INDEX}}p;${{END_INDEX}}p" | \
    tr "\n" " ")

echo "Calculated YAC region: $START to $END"

echo -e "${{TOP_CHR}}\t${{START}}\t${{END}}" > {output.bed}
        '''