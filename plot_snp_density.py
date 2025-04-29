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
    positions = []
    with open(vcf_file, 'r') as f:
        contigs_found = False
        for line in f:
            if line.startswith('#'):
                # Search for the header line to get the length of the genome
                if line.startswith('##contig'):
                    # If there are multiple contigs, throw an error
                    if contigs_found:
                        raise ValueError("Multiple contigs found in VCF file.")
                        # TODO: Handle this case
                    # Mark when a contig line is found, so we can track if there are multiple
                    contigs_found = True
                    # Extract the genome length from the header
                    # We are assuming header line looks like:
                    # ##contig=<ID=NEIAMEAH_1,length=4653851>
                    # So first split on the coma, and then on the equal sign and remove the ">"
                    line_parts = line.strip().split(',')
                    gen_length = int(line_parts[1].split('=')[1][:-1])
                continue  # skip headers
            parts = line.strip().split('\t')
            pos = int(parts[1])  # 2nd column is POS
            positions.append(pos)
    return positions, gen_length

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

    # Harcode skipping samples with multiple contigs, until I figure out how to handle them
    if sample_name not in ["B9NRXP_1_S105",
                          "B9NRXP_2_S106",
                          "TKDHS4_5_S37",
                          "TKDHS4_19_S51",
                          "TKDHS4_22_S54"]:
        print(f"Skipping {sample_name} due to multiple contigs.")
        continue

    positions, gen_length = parse_vcf(vcf_file)
    bins = bin_positions(positions, BIN_SIZE)

    snp_matrix[sample_name] = bins
    gen_lengths[sample_name] = gen_length

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
sns.heatmap(
    heatmap_df,
    cmap=custom_cmap,
    cbar_kws={'label': 'SNP Count'},
    mask=heatmap_df.isna(),  # Mask NaN values to make them white
    )

plt.xlabel("Genomic Bin Start Position (bp)")
plt.ylabel("Sample")
plt.title("SNP Density Heatmap")
plt.tight_layout()

# Save
plot_path = os.path.join(PLOT_DIR, "snp_density.png")
plt.savefig(plot_path)
plt.close()
