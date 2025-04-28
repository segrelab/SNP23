import os

import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# Set the output directory
OUT_DIR = os.path.dirname(os.path.realpath(__file__))

# Function to extract unique variants
def extract_variants(vcf_file):
    variants = set()
    with open(vcf_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):  # Ignore header lines
                fields = line.strip().split('\t')
                chrom, pos = fields[0], fields[1]  # Extract CHROM and POS
                variants.add(f"{chrom}:{pos}")
    return variants

# Load the VCF files
read1_variants = extract_variants(os.path.join(OUT_DIR, '01-HOT1A3_genbank_read1/output/output.vcf'))
read2_variants = extract_variants(os.path.join(OUT_DIR, '02-HOT1A3_genbank_read2/output/output.vcf'))
both_variants = extract_variants(os.path.join(OUT_DIR, '03-HOT1A3_genbank_both_reads/output/output.vcf'))
fasta_variants = extract_variants(os.path.join(OUT_DIR, '04-HOT1A3_fast_both_reads/output/output.vcf'))
old_variants_whole = extract_variants('/projectnb/hfsp/SNP23/results/raw_files/D20-160028-4500T/D20-160028.vcf')
old_variants_filtered = extract_variants('/projectnb/hfsp/SNP23/results/raw_files/D20-160028-4500T/filtered_D20-160028.vcf')

# Manually update the chromosome names for the genbank results to match the fasta file
# In GenBank the chromosomes are named GHNKPGAE_1 and GHNKPGAE_2
# In the FASTA file the chromosomes are named 2687454166.2687454166 and 2687454167.2687453488
# We will rename the GenBank chromosomes to match the FASTA file
# This is necessary for the Venn diagram to work correctly
both_variants_chrom_rename = {f"{chrom.replace('GHNKPGAE_1', '2687454166.2687454166')}" for chrom in both_variants}
both_variants_chrom_rename = {f"{chrom.replace('GHNKPGAE_2', '2687454167.2687453488')}" for chrom in both_variants_chrom_rename}

# Create the Venn diagram
fig = plt.figure()
venn = venn2([read1_variants, read2_variants], ('Read 1 (Trimmed)', 'Read 2 (Trimmed)'))
plt.title("Venn Diagram of VCF Variants")
plt.savefig(os.path.join(OUT_DIR, 'read1_vs_read2.png'))

# Combine the variants from both reads
combined_variants = read1_variants.union(read2_variants)
fig = plt.figure()
venn = venn2([combined_variants, both_variants], ('Combined Reads', 'Both Reads'))
plt.title("Venn Diagram of VCF Variants")
plt.savefig(os.path.join(OUT_DIR, 'combined_vs_both.png'))

# Compare the Genbank and FASTA results
fig = plt.figure()
venn = venn2([both_variants_chrom_rename, fasta_variants], ('Genbank', 'FASTA'))
plt.title("Venn Diagram of VCF Variants")
plt.savefig(os.path.join(OUT_DIR, 'genbank_vs_fasta_fixed.png'))

# With the renamed chromosomes
fig = plt.figure()
venn = venn2([both_variants, fasta_variants], ('Genbank', 'FASTA'))
plt.title("Venn Diagram of VCF Variants")
plt.savefig(os.path.join(OUT_DIR, 'genbank_vs_fasta.png'))

# Compare the old and new results
fig = plt.figure()
venn = venn2([fasta_variants, old_variants_filtered], ('breseq (on fasta reference )', 'Old Pipeline (Filtered)'))
plt.title("Venn Diagram of VCF Variants")
plt.savefig(os.path.join(OUT_DIR, 'breseq_vs_old.png'))

fig = plt.figure()
venn = venn2([fasta_variants, old_variants_whole], ('breseq (on fasta reference )', 'Old Pipeline (All)'))
plt.title("Venn Diagram of VCF Variants")
plt.savefig(os.path.join(OUT_DIR, 'breseq_vs_old_unfiltered.png'))
