import pyBigWig
import numpy as np
import argparse


def extract_all(bw_path, yac_chrom):
    bw = pyBigWig.open(bw_path)
    arrays = {}
    for chrom, length in bw.chroms().items():
        vals = bw.values(chrom, 0, length, numpy=True)
        key = "yac" if chrom == yac_chrom else chrom
        arrays[key] = np.nan_to_num(vals)
    bw.close()
    return arrays


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert BigWig coverage to NumPy arrays (all chromosomes)')
    parser.add_argument('--fwd_bw', required=True, help='Forward strand BigWig file')
    parser.add_argument('--rev_bw', required=True, help='Reverse strand BigWig file')
    parser.add_argument('--bed', required=True, help='YAC region BED file')
    parser.add_argument('--fwd_output', required=True, help='Output forward .npz file path')
    parser.add_argument('--rev_output', required=True, help='Output reverse .npz file path')
    args = parser.parse_args()

    with open(args.bed) as f:
        parts = f.read().strip().split('\t')
        yac_chrom = parts[0]

    fwd_arrays = extract_all(args.fwd_bw, yac_chrom)
    rev_arrays = extract_all(args.rev_bw, yac_chrom)

    np.savez(args.fwd_output, **fwd_arrays)
    np.savez(args.rev_output, **rev_arrays)
    print(f"Saved {len(fwd_arrays)} arrays to {args.fwd_output}")
    print(f"Saved {len(rev_arrays)} arrays to {args.rev_output}")
