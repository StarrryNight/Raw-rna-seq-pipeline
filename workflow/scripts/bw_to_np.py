import pyBigWig
import numpy as np
import argparse


def extract_all(bw_path, yac_chrom, yac_start, yac_end, suffix):
    bw = pyBigWig.open(bw_path)
    arrays = {}
    for chrom, length in bw.chroms().items():
        if chrom == yac_chrom:
            vals = bw.values(chrom, yac_start, yac_end, numpy=True)
            key = "yac"
        else:
            vals = bw.values(chrom, 0, length, numpy=True)
            key = chrom
        arrays[f"{key}_{suffix}"] = np.nan_to_num(vals)
    bw.close()
    return arrays


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert BigWig coverage to NumPy arrays (all chromosomes)')
    parser.add_argument('--fwd_bw', required=True, help='Forward strand BigWig file')
    parser.add_argument('--rev_bw', required=True, help='Reverse strand BigWig file')
    parser.add_argument('--bed', required=True, help='YAC region BED file')
    parser.add_argument('--output', required=True, help='Output .npz file path')
    args = parser.parse_args()

    with open(args.bed) as f:
        parts = f.read().strip().split('\t')
        yac_chrom, yac_start, yac_end = parts[0], int(parts[1]), int(parts[2])

    arrays = {}
    arrays.update(extract_all(args.fwd_bw, yac_chrom, yac_start, yac_end, "fwd"))
    arrays.update(extract_all(args.rev_bw, yac_chrom, yac_start, yac_end, "rev"))
    np.savez(args.output, **arrays)
    print(f"Saved {len(arrays)} arrays to {args.output}")
