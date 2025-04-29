import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

# PARAMETERS
VCF_GLOB = "breseq_results/raw_files/*/output/*.vcf"  # adjust if needed
BIN_SIZE = 1000  # bp window size for SNP density
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
                    # Mark when a contig line is found, so we can track if there are multiple
                    contigs_found = True
                    # Extract the genome length from the header
                    # We are assuming header line looks like:
                    # ##contig=<ID=NEIAMEAH_1,length=4653851>
                    # So first split on the coma, and then on the equal sign and remove the ">"
                    line_parts = line.strip().split(',')
                    gen_length = line_parts[1].split('=')[1][:-1]
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

all_bins = set()
snp_matrix = {}

for vcf_file in all_samples:
    # TODO: Use the metafile to get the strain ID
    sample_name = vcf_file.split(os.sep)[2]  # "TKDHS4_23_S55"
    print(f"Processing {sample_name}")

    positions = parse_vcf(vcf_file)
    bins = bin_positions(positions, BIN_SIZE)

    all_bins.update(bins.keys())
    snp_matrix[sample_name] = bins

# Create a sorted list of all bin starts
all_bins_sorted = sorted(all_bins)

# Create a dataframe with samples as rows, bins as columns
heatmap_df = pd.DataFrame(index=snp_matrix.keys(), columns=all_bins_sorted)

# Fill the dataframe
for sample, bin_counts in snp_matrix.items():
    for bin_start in all_bins_sorted:
        heatmap_df.loc[sample, bin_start] = bin_counts.get(bin_start, 0)

print(heatmap_df)

# Convert to integer (was object due to missing values)
# heatmap_df = heatmap_df.fillna(0).astype(int)

# Plot heatmap
plt.figure(figsize=(20, len(heatmap_df) * 0.5))
sns.heatmap(heatmap_df, cmap="viridis", cbar_kws={'label': 'SNP Count'})

plt.xlabel("Genomic Bin Start Position (bp)")
plt.ylabel("Sample")
plt.title("SNP Density Heatmap")
plt.tight_layout()

# Save
plot_path = os.path.join(PLOT_DIR, "snp_density.png")
plt.savefig(plot_path)
plt.close()
