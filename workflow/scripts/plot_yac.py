import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse


def main():
    parser = argparse.ArgumentParser(description='Plot YAC coverage from numpy arrays')
    parser.add_argument('--fwd_npy', required=True, help='Forward strand numpy array (.npy)')
    parser.add_argument('--rev_npy', required=True, help='Reverse strand numpy array (.npy)')
    parser.add_argument('--bed', required=True, help='YAC region BED file')
    parser.add_argument('--output', required=True, help='Output PNG path')
    args = parser.parse_args()

    fwd = np.load(args.fwd_npy)
    rev = np.load(args.rev_npy)

    with open(args.bed) as f:
        parts = f.read().strip().split('\t')
        chrom, start, end = parts[0], int(parts[1]), int(parts[2])

    positions = np.arange(start, start + len(fwd))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 6), sharex=True)

    ax1.fill_between(positions, fwd, alpha=0.8, color='steelblue', label='Forward')
    ax1.set_ylabel('CPM')
    ax1.set_title(f'YAC Coverage ({chrom}:{start}-{end})')
    ax1.legend()

    ax2.fill_between(positions, rev, alpha=0.8, color='tomato', label='Reverse')
    ax2.set_ylabel('CPM')
    ax2.set_xlabel('Genomic Position (bp)')
    ax2.legend()

    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Plot saved to {args.output}")


if __name__ == '__main__':
    main()
