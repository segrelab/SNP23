#!/bin/bash -l
#$ -l h_rt=75:00:00 # Specify the hard time limit for the job
#$ -j y             # Merge the error and output streams into a single file
#$ -o reads_to_vcf.out # Specify the output file
#$ -P hfsp          # Specify the SCC project name you want to use
#$ -N reads_to_vcf  # Give your job a name
#$ -pe omp 5       # Request multiple slots for the Shared Memory application (OpenMP)

module load bwa/0.7.17 
module load samtools
module load trimmomatic/0.36
module load htslib/1.16
module load bcftools/1.16

# Set the path to the metafile
input_fi="metafile.csv"
skip_first_row=true # A flag to skip the first row, set to false if the column names are not included in the metafile

# Define a file name for the read count spreadsheet
output_csv="results/summary.csv"
# Initialize the read count file with headers
echo "sample_name,reference_genome,input_pairs,surviving_pairs,forward_surviving,reverse_surviving,dropped_pairs,total_reads,mapped_reads,unmapped_reads,avg_coverage,min_coverage,max_coverage,SNPs" > $output_csv

# Set a variable of if you want to rerun the analysis
force_rerun=false

# Function to run the analysis pipeline (from trimming to variant calling)
run_analysis_pipeline() {
    local sample_name=$1
    local label=$2
    local ref=$3
    local in_directory=$4
    local out_directory=$5

    # Setting up variables for inputs, file names and etc
    read_1="$in_directory"/200219Seg_"$label"_1_sequence.fastq.gz
    read_2="$in_directory"/200219Seg_"$label"_2_sequence.fastq.gz

    trimmomatic_log="$out_directory/${label}_trimmomatic.log"

    trim_1="$out_directory"/trimmed_"$label"_1_sequence.fastq.gz
    trimu_1="$out_directory"/trimmedu_"$label"_1_sequence.fastq.gz

    trim_2="$out_directory"/trimmed_"$label"_2_sequence.fastq.gz
    trimu_2="$out_directory"/trimmedu_"$label"_2_sequence.fastq.gz

    bam_file="$out_directory"/"$label".bam
    sorted_bam_file="$out_directory"/"$label".sorted.bam
    depth_file="$out_directory/${label}_depth.txt"

    vcf="$out_directory"/"$label".vcf
    f_vcf="$out_directory"/filtered_"$label".vcf

    # If the output files already exist and force_rerun is false, skip this sample, but add the existing results to the summary table
    if [ -f "$f_vcf" ] && [ $force_rerun = false ]; then
        echo "Output files already exist for $sample_name, skipping..."
        # Add the results to the CSV file
        gather_results_and_write_csv $sample_name $ref $trimmomatic_log $bam_file $depth_file $f_vcf
        # Exit the function without running the rest of the processing for this sample
        return
    fi
    
    ## Actually running tools
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
        # /share/pkg/trimmomatic/0.36/src/Trimmomatic-0.36/adapters/NexteraPE-PE.fa: Path to the adapter sequence file, which contains the adapter sequences that need to be clipped. In this case, the file is for NexteraPE adapters.
        # 2:30:10:1:TRUE: These are parameters for the adapter clipping step:
            # 2: Minimum number of seed mismatches to allow in adapter matching.
            # 30: Palindrome clip threshold for identifying "adapter dimer" artifacts, where the forward and reverse adapters are ligated together.
            # 10: Simple clip threshold; this is used to remove simple adapter sequences.
            # 1: Minimum length of a match that will be clipped.
            # TRUE: Specifies that the reads should be clipped only if both reads (forward and reverse) contain adapters.
    trimmomatic PE -threads 10 $read_1 $read_2 \
        $trim_1 $trimu_1 $trim_2 $trimu_2 \
        ILLUMINACLIP:/share/pkg/trimmomatic/0.36/src/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10:1:TRUE 2> $trimmomatic_log

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

    #  Append the sample name and results to the CSV file
    gather_results_and_write_csv $sample_name $ref $trimmomatic_log $sorted_bam_file $depth_file $f_vcf
}

# Function to gather results and write to the CSV file
gather_results_and_write_csv() {
    local sample_name=$1
    local ref=$2
    local trimmomatic_log=$3
    local bam_file=$4
    local depth_file=$5
    local f_vcf=$6

    # Extract key statistics from the Trimmomatic log file
    input_pairs=$(grep -oP '(?<=Input Read Pairs: )\d+' $trimmomatic_log)
    surviving_pairs=$(grep -oP '(?<=Both Surviving: )\d+' $trimmomatic_log)
    forward_surviving=$(grep -oP '(?<=Forward Only Surviving: )\d+' $trimmomatic_log)
    reverse_surviving=$(grep -oP '(?<=Reverse Only Surviving: )\d+' $trimmomatic_log)
    dropped_pairs=$(grep -oP '(?<=Dropped: )\d+' $trimmomatic_log)

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

# Process each sample in the input file
# These variables must match the order/contents of columns in the input file
while IFS=, read -r sample_number sample_name species_name strain_id in_IAMM accession_number source dna_prep ref_genome ref_path pos_cntrl_genome pos_cntrl_path neg_cntrl_genome neg_cntrl_path;do
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

    # Check that the reference genome file exists
    if [ ! -f "$ref" ]; then
        echo "Reference genome file not found for $sample_name: $ref"
        continue
    fi

    # Get the sample label, the sample_name minus the -4500T, eg D20-160027
    label="${sample_name%-4500T}"

    # Get the directory to the raw results (in this repo)
    # TODO: Make this a variable that I set at the top of the script
    directory=results/raw_files/$sample_name

    # Run the analysis pipeline for the current sample
    run_analysis_pipeline $sample_name $label $ref $directory $directory

    # Process positive control if present
     if [ -n "$pos_cntrl_genome" ] && [ -n "$pos_cntrl_path" ]; then
        pos_cntrl_path_trimmed=$(echo "$pos_cntrl_path" | xargs | tr -d '\r')
        pos_cntrl_genome_trimmed=$(echo "$pos_cntrl_genome" | xargs | tr -d '\r')
        pos_cntrl_ref="${pos_cntrl_path_trimmed}/${pos_cntrl_genome_trimmed}"
        pos_cntrl_directory="results/raw_files/$sample_name/controls/positive_control"  # TODO: Make this a variable that I set at the top of the script

        mkdir -p "$pos_cntrl_directory"
        run_analysis_pipeline "${sample_name}_pos_cntrl" "$label" "$pos_cntrl_ref" "$directory" "$pos_cntrl_directory"
    fi

    # Process negative control if present
    if [ -n "$neg_cntrl_genome" ] && [ -n "$neg_cntrl_path" ]; then
        neg_cntrl_path_trimmed=$(echo "$neg_cntrl_path" | xargs | tr -d '\r')
        neg_cntrl_genome_trimmed=$(echo "$neg_cntrl_genome" | xargs | tr -d '\r')
        neg_cntrl_ref="${neg_cntrl_path_trimmed}/${neg_cntrl_genome_trimmed}"
        neg_cntrl_directory="results/raw_files/$sample_name/controls/negative_control"

        mkdir -p "$neg_cntrl_directory"
        echo "Running analysis pipeline for negative control:"
        echo $neg_cntrl_directory
        echo $neg_cntrl_ref
        run_analysis_pipeline "${sample_name}_neg_cntrl" "$label" "$neg_cntrl_ref" "$directory" "$neg_cntrl_directory"
    fi

done < "$input_fi"