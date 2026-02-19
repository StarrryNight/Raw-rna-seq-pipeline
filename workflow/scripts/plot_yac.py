import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse

YEAST_CHROMS = [
    'chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 'chrVII', 'chrVIII',
    'chrIX', 'chrX', 'chrXI', 'chrXII', 'chrXIII', 'chrXIV', 'chrXV', 'chrXVI'
]


def plot_yac_coverage(fwd_data, rev_data, chrom, start, end, output_path):
    fwd = fwd_data['yac']
    rev = rev_data['yac']
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
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"YAC plot saved to {output_path}")


def plot_yeast_grid(data, color, output_path, strand_label):
    fig, axes = plt.subplots(4, 4, figsize=(20, 16))
    axes = axes.flatten()

    for i, chrom in enumerate(YEAST_CHROMS):
        ax = axes[i]
        if chrom in data:
            vals = data[chrom]
            ax.fill_between(np.arange(len(vals)), vals, alpha=0.8, color=color)
        ax.set_title(chrom, fontsize=10)
        ax.set_xlabel('Position (bp)', fontsize=8)
        ax.set_ylabel('CPM', fontsize=8)
        ax.tick_params(labelsize=7)

    fig.suptitle(f'Yeast Chromosome Coverage — {strand_label} Strand', fontsize=14)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Yeast {strand_label} grid saved to {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Plot YAC and yeast coverage from npz')
    parser.add_argument('--fwd_npz', required=True, help='Forward strand .npz file')
    parser.add_argument('--rev_npz', required=True, help='Reverse strand .npz file')
    parser.add_argument('--bed', required=True, help='YAC region BED file')
    parser.add_argument('--out_yac', required=True, help='Output YAC coverage PNG')
    parser.add_argument('--out_yeast_fwd', required=True, help='Output yeast fwd grid PNG')
    parser.add_argument('--out_yeast_rev', required=True, help='Output yeast rev grid PNG')
    args = parser.parse_args()

    with open(args.bed) as f:
        parts = f.read().strip().split('\t')
        chrom, start, end = parts[0], int(parts[1]), int(parts[2])

    fwd_data = np.load(args.fwd_npz)
    rev_data = np.load(args.rev_npz)

    plot_yac_coverage(fwd_data, rev_data, chrom, start, end, args.out_yac)
    plot_yeast_grid(fwd_data, 'steelblue', args.out_yeast_fwd, 'Forward')
    plot_yeast_grid(rev_data, 'tomato', args.out_yeast_rev, 'Reverse')


if __name__ == '__main__':
    main()
