import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from collections import defaultdict

# PARAMETERS
VCF_GLOB = "breseq_results/raw_files/*/output/*.vcf"  # adjust if needed
BIN_SIZE = 10000  # bp window size for SNP density
PLOT_DIR = "breseq_results/plots"

os.makedirs(PLOT_DIR, exist_ok=True)

# Read the VCF file and extract the genome length and the position of all SNPs
def parse_vcf(vcf_file):
    contigs = {}
    current_contig = None

    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                # Parse contig metadata from the header
                if line.startswith('##contig'):
                    line_parts = line.strip().split(',')
                    contig_id = line_parts[0].split('ID=')[1]  # Extract contig ID
                    contig_length = int(line_parts[1].split('length=')[1][:-1])  # Extract contig length
                    contigs[contig_id] = {'length': contig_length, 'positions': []}
                continue  # Skip other header lines

            # Parse SNP entries
            parts = line.strip().split('\t')
            chrom = parts[0]  # CHROM column
            pos = int(parts[1])  # POS column

            # Add SNP position to the correct contig
            if chrom in contigs:
                contigs[chrom]['positions'].append(pos)

    return contigs

def bin_positions(positions, bin_size):
    bins = defaultdict(int)
    for pos in positions:
        bin_start = (pos // bin_size) * bin_size
        bins[bin_start] += 1
    return bins

# Main
all_samples = glob.glob(VCF_GLOB)
print(f"Found {len(all_samples)} VCF files.")

snp_matrix = {}
gen_lengths = {}

for vcf_file in all_samples:
    # TODO: Use the metafile to get the strain ID
    sample_name = vcf_file.split(os.sep)[2]  # "TKDHS4_23_S55"
    print(f"Processing {sample_name}")

    contigs = parse_vcf(vcf_file)

    for contig in contigs:
        positions = contigs[contig]['positions']
        gen_length = contigs[contig]['length']

        # Bin the SNP positions
        bins = bin_positions(positions, BIN_SIZE)

        # Add the binned SNP counts to the matrix
        snp_matrix[sample_name + "_" + contig] = bins
        gen_lengths[sample_name + "_" + contig] = gen_length

# Create a sorted list of all possible bin starts
# Get the list of all possible bin starts from the length of the longest genome
all_bins = set()
all_bins.update(range(0, max(gen_lengths.values())+1, BIN_SIZE))
all_bins_sorted = sorted(all_bins)

# Create a dataframe with samples as rows, bins as columns
heatmap_df = pd.DataFrame(index=snp_matrix.keys(), columns=all_bins_sorted)

# Fill the dataframe
for sample, bin_counts in snp_matrix.items():
    for bin_start in all_bins_sorted:
        if bin_start in bin_counts:
            heatmap_df.loc[sample, bin_start] = bin_counts[bin_start]
        elif bin_start > gen_lengths[sample]:
            # If the bin start is greater than the genome length, set it to NaN
            heatmap_df.loc[sample, bin_start] = float('nan')
        else:
            # If the bin start is not in the counts, set it to 0
            heatmap_df.loc[sample, bin_start] = 0

# Ensure the DataFrame is of numeric type
heatmap_df = heatmap_df.astype(float)

# Define my own color map, with 0 as gray, and then the viridis color map
viridis = plt.get_cmap("viridis", 256)  # Get the viridis colormap
new_colors = viridis(np.linspace(0, 1, 256))  # Extract the colors
new_colors[0] = np.array([0.5, 0.5, 0.5, 1.0])  # Replace the first color with gray (R=0.5, G=0.5, B=0.5, A=1.0)
custom_cmap = mcolors.ListedColormap(new_colors)

# Plot heatmap
plt.figure(figsize=(20, len(heatmap_df) * 0.5))
ax = sns.heatmap(
    heatmap_df,
    cmap=custom_cmap,
    cbar_kws={'label': 'SNP Count'},
    mask=heatmap_df.isna(),  # Mask NaN values to make them white,
    )

# Add white lines between rows
for y in range(1, heatmap_df.shape[0]):  # Iterate over row indices
    ax.hlines(y=y, xmin=0, xmax=heatmap_df.shape[1], color='white', linewidth=2)

# Only plot X-ticks for the megabase positions
xticks = np.arange(0, len(heatmap_df.columns), 100)
plt.xticks(xticks, [f"{int(x * BIN_SIZE / 1e6)} Mb" for x in xticks], rotation=45)

plt.xlabel("Genomic Bin Start Position (bp)")
plt.ylabel("Sample")
plt.title("SNP Density Heatmap")
plt.tight_layout()

# Save
plot_path = os.path.join(PLOT_DIR, "snp_density.png")
plt.savefig(plot_path)
plt.close()
