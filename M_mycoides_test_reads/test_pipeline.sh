#!/bin/bash -l
#$ -l h_rt=75:00:00 # Specify the hard time limit for the job
#$ -j y             # Merge the error and output streams into a single file
#$ -P hfsp          # Specify the SCC project name you want to use
#$ -N M_mycoides_reads_to_vcf  # Give your job a name
#$ -pe omp 5       # Request multiple slots for the Shared Memory application (OpenMP)

module load bwa/0.7.17 
module load samtools
module load trimmomatic/0.36
module load htslib/1.16
module load bcftools/1.16

# Set path to the reads and the genome
ref='/projectnb/hfsp/SNP23/M_mycoides_test_reads/CP014940.1.fasta'
reads='/projectnb/hfsp/SNP23/M_mycoides_test_reads/SRR15032631.fastq'

# Define a file name for the read count spreadsheet
output_csv="/projectnb/hfsp/SNP23/M_mycoides_test_reads/summary.csv"
# Initialize the read count file with headers
echo "sample_name,reference_genome,input_pairs,surviving_pairs,forward_surviving,reverse_surviving,dropped_pairs,total_reads,mapped_reads,unmapped_reads,avg_coverage,min_coverage,max_coverage,SNPs" > $output_csv

# Varaibles for the analysis pipeline
out_directory="/projectnb/hfsp/SNP23/M_mycoides_test_reads/raw_results"
label="M_mycoides"

# Define the files that will be made in the pipeline
bam_file="$out_directory"/"$label".bam
sorted_bam_file="$out_directory"/"$label".sorted.bam
depth_file="$out_directory/${label}_depth.txt"

vcf="$out_directory"/"$label".vcf
f_vcf="$out_directory"/filtered_"$label".vcf

# Run the pipeline
# Create an index for the reference genome
# When you run bwa index, BWA will create several index files, typically with extensions like .amb, .ann, .bwt, .pac, and .sa. These files are saved in the same directory as the reference genome and are used in the alignment step to map sequencing reads to the reference.
bwa index $ref

# Use samtools to generate a FASTA index for the reference genome
# This will create a .fai file. This file contains the byte offsets and lengths of each sequence in the FASTA file, allowing for quick lookup of specific sequences.
# Indexing the reference genome is required for efficient random access to the genome's sequences. Many bioinformatics tools, including Samtools, BWA, and bcftools, use the .fai index file to quickly retrieve specific regions of the reference genome during tasks like read alignment or variant calling. Without the index, every operation would require reloading the entire genome into memory, which would be slow and inefficient for large genomes.
samtools faidx $ref

# Align the reads with bwa and convert to a bam file using samtools
# The BWA MEM algorithm aligns the paired-end reads ($trim_1 and $trim_2) to the reference genome ($ref) and produces alignments in the SAM (Sequence Alignment/Map) format.
bwa mem $ref $reads | samtools view -S -b > 'raw_results/M_mycoides.bam'

# Sort the bam file by genomic coordinates and save the sorted output to a new file
samtools sort $bam_file -o $sorted_bam_file

# Get the coverage depth and save it to a file
samtools depth -a $sorted_bam_file > $depth_file

# Generate a VCF file with bcftools
# bcftools mpileup: Generates a pilup of reads aligned to a reference genome. The pileup format summarizes the base calls at each position of the reference genome, which is then used to detect variants (such as SNPs or indels).
    # -B: which is used to improve the accuracy of indel calling. If you remove this, the pileup will be faster but less accurate in detecting indels.
    # -q: Filters out reads with a mapping quality lower than the given number (e.g. 30). Reads with low mapping quality are more likely to be incorrectly aligned, so this helps to reduce noise in the variant calling process.
    # FIXME: Check what -g really is in the documentation. ChatGPT and copilot disagreed.
    # -g: Specifies the maximum per-sample depth to use. This avoids excessive computation at positions with very high coverage (more than 50 reads at a single position).
    # -Ou: Output the pileup in uncompressed BCF format, which is more efficient for piping into the next command (bcftools call).
    # -f $ref: Specifies the reference genome file (FASTA format) that the reads are aligned to. The $ref variable points to the reference genome.
# bcftools call: Performs the variant calling
    # --ploidy 1: Specifies that the organism is haploid (1 copy of each chromosome).
    # -m: Use the multiallelic caller model. This model can call multiple alleles at the same position, which is useful for handling complex variants.
    # -A (--keep-alts): Output all alternate alleles present in the alignments even if they do not appear in any of the genotypes
bcftools mpileup -B -q 30 -g 50 -Ou -f $ref $sorted_bam_file | bcftools call --ploidy 1 -m -A > $vcf

# Filter the VCF file for all variants with a frequency of at least 50% at a given position
# bcftools view: View, subset and filter VCF or BCF files by position and filtering expression.
    # -i (--include): Include only sites for which the following expression is true.
        # '(DP4[2]+DP4[3])/sum(DP4) >= 0.5 & sum(DP4)>=5': This expression filters the VCF file to include only sites where the frequency of the alternate allele is at least 50% (i.e., the alternate allele is present in at least half of the reads at that position) and the total depth of coverage is at least 5 reads.
            # DP4 is a field in VCF files that contains the high-quality read depths for each of four categories:
                # DP4[0]: Reference reads on the forward strand.
                # DP4[1]: Reference reads on the reverse strand.
                # DP4[2]: Alternate reads on the forward strand.
                # DP4[3]: Alternate reads on the reverse strand.
            # DP4[2] + DP4[3] calculates the total number of reads supporting the alternate allele (both forward and reverse strands).
            # sum(DP4) calculates the total depth of coverage at a given position, including both reference and alternate reads.
            # (DP4[2] + DP4[3]) / sum(DP4) >= 0.5  keeps only variants where at least 50% of the total reads support the alternate allele.
            # sum(DP4) >= 5 ensures that the total depth of coverage is at least 5 reads.
# TODO: Does this keep only SNPs?
bcftools view -i '(DP4[2]+DP4[3])/sum(DP4) >= 0.5 & sum(DP4)>=5' $vcf > $f_vcf

# Function to gather results and write to the CSV file
gather_results_and_write_csv() {
    local sample_name=$1
    local ref=$2
    local trimmomatic_log=$3
    local bam_file=$4
    local depth_file=$5
    local f_vcf=$6

    # Extract key statistics from the Trimmomatic log file
    input_pairs="NA"
    surviving_pairs="NA"
    forward_surviving="NA"
    reverse_surviving="NA"
    dropped_pairs="NA"

    # Get total, mapped, and unmapped reads
    total_reads=$(samtools view -c $bam_file) # -c is the option that tells samtools view to count the total number of reads in the BAM file.
    mapped_reads=$(samtools view -c -F 4 $bam_file) # -F 4 excludes reads that have the "unmapped" flag set (so it only counts mapped reads).
    unmapped_reads=$(samtools view -c -f 4 $bam_file) # -f 4 includes only reads that have the "unmapped" flag set.

    # Calculate coverage statistics
    if [ -s $depth_file ]; then
        avg_coverage=$(awk '{sum+=$3} END {if (NR>0) print sum/NR; else print 0}' $depth_file)
        min_coverage=$(awk 'NR == 1 || $3 < min {min = $3} END {if (NR>0) print min; else print 0}' $depth_file)
        max_coverage=$(awk 'NR == 1 || $3 > max {max = $3} END {if (NR>0) print max; else print 0}' $depth_file)
    else
        avg_coverage=0
        min_coverage=0
        max_coverage=0
    fi

    # Get the number of SNPs in the filtered VCF file
    # -H option is used to suppress the header lines in the output, so only the SNP records are counted.
    snp_count=$(bcftools view -H $f_vcf | wc -l)

    # Append results to CSV
    echo "$sample_name,$ref,$input_pairs,$surviving_pairs,$forward_surviving,$reverse_surviving,$dropped_pairs,$total_reads,$mapped_reads,$unmapped_reads,$avg_coverage,$min_coverage,$max_coverage,$snp_count" >> $output_csv
}

#  Append the sample name and results to the CSV file
gather_results_and_write_csv $sample_name $ref "NONE" $sorted_bam_file $depth_file $f_vcf

