checkpoint find_yac_region:
    input:
        bam = f"{RESULTS}/{{sample}}/aligned/Aligned.sortedByCoord.out.bam",
        bai = f"{RESULTS}/{{sample}}/aligned/Aligned.sortedByCoord.out.bam.bai"
    output:
        bed = f"{RESULTS}/{{sample}}/yac/yac_region.bed"
    resources:
        mem_mb = 4000,
        time = "1:00:00"
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
START_INDEX=$(( NUM_READS / 100 ))
END_INDEX=$(( NUM_READS * 99 / 100 ))

#TODO rename duplicated chrom names
read START END <<< $(samtools view {input.bam} "$TOP_CHR" 2>/dev/null | \
    awk '{{print $4}}' | \
    sort -n | \
    sed -n "${START_INDEX}p;${END_INDEX}p" | \
    tr "\n" " ")

echo "Calculated YAC region: $START to $END"

echo -e "${TOP_CHR}\t${START}\t${END}" > {output.bed}
        '''