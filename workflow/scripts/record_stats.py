import argparse, json, subprocess, sys

YEAST_CHROMS = {
    'chrI','chrII','chrIII','chrIV','chrV','chrVI','chrVII','chrVIII',
    'chrIX','chrX','chrXI','chrXII','chrXIII','chrXIV','chrXV','chrXVI'
}

def parse_star_log(log_path):
    total, unmapped = 0, 0
    with open(log_path) as f:
        for line in f:
            if 'Number of input reads' in line:
                total = int(line.split('|')[1].strip())
            elif 'Number of reads unmapped:' in line:
                unmapped += int(line.split('|')[1].strip())
    return total, unmapped

def idxstats(bam_path):
    result = subprocess.run(
        ['samtools', 'idxstats', bam_path],
        capture_output=True, text=True, check=True
    )
    counts = {}
    for line in result.stdout.strip().splitlines():
        chrom, length, mapped, unmapped = line.split('\t')
        counts[chrom] = int(mapped)
    return counts

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample',  required=True)
    parser.add_argument('--log',     required=True, help='alignment2 Log.final.out')
    parser.add_argument('--bam',     required=True, help='alignment2 BAM for idxstats')
    parser.add_argument('--bed',     required=True, help='yac_region.bed')
    parser.add_argument('--fwd_bw',  required=True)
    parser.add_argument('--rev_bw',  required=True)
    parser.add_argument('--output',  required=True, help='output stats.json path')
    args = parser.parse_args()

    with open(args.bed) as f:
        parts = f.read().strip().split('\t')
        yac_chrom, yac_start, yac_end = parts[0], int(parts[1]), int(parts[2])

    total_reads, unmapped_reads = parse_star_log(args.log)
    chrom_counts = idxstats(args.bam)

    yac_reads = chrom_counts.get(yac_chrom, 0)
    yeast_reads = 0
    for chrom, count in chrom_counts.items():
        if chrom == yac_chrom or chrom == 'chrM' or chrom == '*':
            continue
        if chrom in YEAST_CHROMS:
            yeast_reads += count
        else:
            print(f"Warning: unexpected chrom '{chrom}' skipped", file=sys.stderr)

    def pct(n):
        return round(n / total_reads * 100, 2) if total_reads else 0.0

    stats = {
        "name":                      args.sample,
        "chromosome":                yac_chrom,
        "location":                  f"{yac_start}-{yac_end}",
        "sources":                   {"fwd": args.fwd_bw, "rev": args.rev_bw},
        "total_reads":               total_reads,
        "reads_mapped_to_yac":       yac_reads,
        "reads_mapped_to_yac_pct":   pct(yac_reads),
        "reads_mapped_to_yeast":     yeast_reads,
        "reads_mapped_to_yeast_pct": pct(yeast_reads),
        "unmapped_reads":            unmapped_reads,
        "unmapped_reads_pct":        pct(unmapped_reads),
    }

    with open(args.output, 'w') as f:
        json.dump(stats, f, indent=2)
    print(f"Stats written to {args.output}")
