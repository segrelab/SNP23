#!/bin/bash -l
#$ -l h_rt=75:00:00
#$ -j y
#$ -P hfsp
#$ -pe omp 10

module load bwa/0.7.17 
module load samtools
module load trimmomatic/0.36
module load htslib/1.16
module load bcftools/1.16

# Set the path to the metafile
input_fi="mini_metafile.csv"
skip_first_row=true # A flag to skip the first row, set to false if the column names are not included in the metafile

# Define a file name for the read count spreadsheet
read_count_file="read_counts.csv"
# Initialize the read count file with headers
echo "sample_name,total_reads,mapped_reads,unmapped_reads" > $read_count_file

# These variables must match the order/contents of columns in the input file
while IFS=, read -r sample_number sample_name species_name strain_id in_IAMM accession_number source dna_prep ref_genome ref_path;do
    # Skip the first row (not a real sample, just a header)
    if $skip_first_row; then
        skip_first_row=false
        continue
    fi

    # Trim any leading/trailing whitespace including special characters
    ref_path_trimmed=$(echo "$ref_path" | xargs | tr -d '\r')
    ref_genome_trimmed=$(echo "$ref_genome" | xargs | tr -d '\r')

    # Remove any trailing slashes from ref_path_trimmed
    ref_path_trimmed="${ref_path_trimmed%/}"

    # If the reference genome path is empty, skip it
    if [ -z "$ref_path_trimmed" ]; then
        continue
    fi

    # Create a single varaible for the refernece genome path and file name
    ref="${ref_path_trimmed}/${ref_genome_trimmed}"

    # Get the sample label, the sample_name minus the -4500T, eg D20-160027
    label="${sample_name%-4500T}"

    # Setting up variables for inputs, file names and etc
    directory=/projectnb2/hfsp/Challenge21/joseline/data_for_analysis/$sample_name
    read_1="$directory"/200219Seg_"$label"_1_sequence.fastq.gz
    read_2="$directory"/200219Seg_"$label"_2_sequence.fastq.gz

    trim_1="$directory"/trimmed_"$label"_1_sequence.fastq.gz
    trimu_1="$directory"/trimmedu_"$label"_1_sequence.fastq.gz

    trim_2="$directory"/trimmed_"$label"_2_sequence.fastq.gz
    trimu_2="$directory"/trimmedu_"$label"_2_sequence.fastq.gz

    bam_file="$directory"/"$label".bam
    sorted_bam_file="$directory"/"$label".sorted.bam
    vcf="$directory"/"$label".vcf
    f_vcf="$directory"/filtered_"$label".vcf

    ## Actually runnning tools
    # Trimming
    # trimmomatic PE: Invokes Trimmomatic in paired-end mode (PE). Paired-end mode means the tool is processing two files of paired sequencing reads (typically from two ends of a DNA fragment).
    # -threads 10: Specifies the number of CPU threads to use for parallel processing. In this case, 10 threads are being used to speed up the trimming process.
    # $read_1: The file containing the first set of paired-end reads (forward reads).
    # $read_2: The file containing the second set of paired-end reads (reverse reads).
    # $trim_1: The output file for the trimmed forward reads (surviving paired reads from the first file).
    # $trimu_1: The output file for unpaired reads that were originally part of the forward reads but lost their pair during trimming.
    # $trim_2: The output file for the trimmed reverse reads (surviving paired reads from the second file).
    # $trimu_2: The output file for unpaired reads that were originally part of the reverse reads but lost their pair during trimming.
    # ILLUMINACLIP:/usr4/bf527/smit2/.conda/pkgs/trimmomatic-0.39-1/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10:1:TRUE:
        # ILLUMINACLIP: This is a specific step in Trimmomatic to clip (remove) adapter sequences from the reads. Adapters are sequences added to the ends of DNA fragments during library preparation and can interfere with downstream analysis.
        # /usr4/bf527/smit2/.conda/pkgs/trimmomatic-0.39-1/share/trimmomatic/adapters/NexteraPE-PE.fa: Path to the adapter sequence file, which contains the adapter sequences that need to be clipped. In this case, the file is for NexteraPE adapters.
        # 2:30:10:1:TRUE: These are parameters for the adapter clipping step:
            # 2: Minimum number of seed mismatches to allow in adapter matching.
            # 30: Palindrome clip threshold for identifying "adapter dimer" artifacts, where the forward and reverse adapters are ligated together.
            # 10: Simple clip threshold; this is used to remove simple adapter sequences.
            # 1: Minimum length of a match that will be clipped.
            # TRUE: Specifies that the reads should be clipped only if both reads (forward and reverse) contain adapters.
    # FIXME: The file usr4/bf527/smit2/.conda/pkgs/trimmomatic-0.39-1/share/trimmomatic/adapters/NexteraPE-PE.fa is not found
    trimmomatic PE -threads 10 $read_1 $read_2 \ $trim_1 $trimu_1 \ $trim_2 $trimu_2 \ ILLUMINACLIP:/usr4/bf527/smit2/.conda/pkgs/trimmomatic-0.39-1/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10:1:TRUE

    # Create an index for the reference genome
    # When you run bwa index, BWA will create several index files, typically with extensions like .amb, .ann, .bwt, .pac, and .sa. These files are saved in the same directory as the reference genome and are used in the alignment step to map sequencing reads to the reference.
    bwa index $ref

    # Use samtools to generate a FASTA index for the reference genome
    # This will create a .fai file. This file contains the byte offsets and lengths of each sequence in the FASTA file, allowing for quick lookup of specific sequences.
    # Indexing the reference genome is required for efficient random access to the genome's sequences. Many bioinformatics tools, including Samtools, BWA, and bcftools, use the .fai index file to quickly retrieve specific regions of the reference genome during tasks like read alignment or variant calling. Without the index, every operation would require reloading the entire genome into memory, which would be slow and inefficient for large genomes.
    samtools faidx $ref

    # Align the reads with bwa and convert to a bam file using samtools
    # The BWA MEM algorithm aligns the paired-end reads ($trim_1 and $trim_2) to the reference genome ($ref) and produces alignments in the SAM (Sequence Alignment/Map) format.
    bwa mem $ref $trim_1 $trim_2 | samtools view -S -b > $bam_file

    # Get total, mapped, and unmapped reads
    total_reads=$(samtools view -c $bam_file) # -c is the option that tells samtools view to count the total number of reads in the BAM file.
    mapped_reads=$(samtools view -c -F 4 $bam_file) # -F 4 excludes reads that have the "unmapped" flag set (so it only counts mapped reads).
    unmapped_reads=$(samtools view -c -f 4 $bam_file) # -f 4 includes only reads that have the "unmapped" flag set.

    # Append the sample name and counts to the CSV file
    echo "$sample_name,$total_reads,$mapped_reads,$unmapped_reads" >> $output_csv

    # Sort the bam file by genomic coordinates and save the sorted output to a new file
    samtools sort $bam_file -o $sorted_bam_file

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
        # FIXME: Check what -A really is in the documentation. ChatGPT and copilot disagreed.
        # -A: Output all sites, including non-variant sites. This is useful for generating a complete VCF file with all positions in the reference genome.
    bcftools mpileup -B -q 30 -g 50 -Ou -f $ref $sorted_bam_file | bcftools call --ploidy 1 -m -A > $vcf

    # Filter the VCF file to remove low-quality variants and retain only SNPs
    # bcftools filter: Filters the VCF file to remove variants marked as "LowQual"
        # -s LowQual: Specifies the filter to apply. Variants marked as "LowQual" will be removed.
        # bcftools view -v snps: Filters the VCF file to retain only SNPs (single nucleotide polymorphisms).
    # bcftools view: Filters the VCF file to retain only SNPs (single nucleotide polymorphisms).
        # -v snps: Specifies that only SNPs should be retained in the output.
    bcftools filter -s LowQual $vcf | bcftools view -v snps > $f_vcf
done < "$input_fi"