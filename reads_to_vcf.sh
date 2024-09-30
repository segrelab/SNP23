#!/bin/bash -l
#$ -l h_rt=75:00:00
#$ -j y
#$ -P hfsp
#$ -pe omp 10

module load bwa/0.7.17 
module load samtools
module load trimmomatic/0.36
module load bcftools
# Set the path to the metafile
input_fi="metafile.csv"

# These must match the order of columns in the input file
while IFS=, read -r sample_number sample_name species_name strain_id in_IAMM accession_number source dna_prep ref_genome;do
    # TODO: Skip the first row (not a real sample, just a header)

    #setting up variables for inputs, file names and etc
    # TODO: Replace with the original results folder /projectnb/hfsp/Strain_Library/Raw_Illumina_200219Seg/
    directory=/projectnb2/hfsp/Challenge21/joseline/data_for_analysis/$label-4500T/
    read_1="$directory"200219Seg_"$label"_1_sequence.fastq.gz
    read_2="$directory"200219Seg_"$label"_2_sequence.fastq.gz

    # TODO: Should these files be version controlled? Or saved in the original folder?
    trim_1="$directory"trimmed_"$label"_1_sequence.fastq.gz
    trimu_1="$directory"trimmedu_"$label"_1_sequence.fastq.gz

    trim_2="$directory"trimmed_"$label"_2_sequence.fastq.gz
    trimu_2="$directory"trimmedu_"$label"_2_sequence.fastq.gz

    bam_file="$directory""$label".bam
    sorted_bam_file="$directory""$label".sorted.bam
    vcf="$directory""$label".vcf
    f_vcf="$directory"filtered_"$label".vcf

  #actually runnning tools 
    trimmomatic PE -threads 10 $read_1 $read_2 \ $trim_1 $trimu_1 \ $trim_2 $trimu_2 \ ILLUMINACLIP:/usr4/bf527/smit2/.conda/pkgs/trimmomatic-0.39-1/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10:1:TRUE

    bwa index $ref
    samtools faidx $ref
    bwa mem $ref $trim_1 $trim_2 | samtools view -S -b > $bam_file
    samtools sort $bam_file -o $sorted_bam_file
    
    samtools faidx $ref
    bcftools mpileup -B -q 30 -g 50 -Ou -f $ref $sorted_bam_file | bcftools call --ploidy 1 -m -A > $vcf

    bcftools filter -s LowQual $vcf | bcftools view -v snps > $f_vcf
done < $input_fi