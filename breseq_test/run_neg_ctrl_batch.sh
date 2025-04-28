#!/bin/bash -l

#$ -P hfsp              # Specify the SCC project name you want to use
#$ -l h_rt=12:00:00     # Specify the hard time limit for the job
#$ -N hot1a3-neg-ctrl   # Give job a name
#$ -j y                 # Merge the error and output streams into a single file

module load miniconda
conda activate /projectnb/hfsp/SNP23/envs/breseq

# Commenting out because I already ran it
# breseq -r /projectnb/hfsp/IAMM_genomes_from_Luca/D20-160028_assembly.fasta -o /projectnb/hfsp/SNP23/breseq_test/05-HOT1A3_neg_cntrl results/raw_files/D20-160028-4500T/trimmed_D20-160028_1_sequence.fastq.gz results/raw_files/D20-160028-4500T/trimmed_D20-160028_2_sequence.fastq.gz

breseq -r /projectnb/hfsp/IAMM_genomes_from_Luca/D20-160033_assembly.fasta -o /projectnb/hfsp/SNP23/breseq_test/06-Luca5_neg_ctrl /projectnb/hfsp/SNP23/results/raw_files/D20-160033-4500T/trimmed_D20-160033_1_sequence.fastq.gz /projectnb/hfsp/SNP23/results/raw_files/D20-160033-4500T/trimmed_D20-160033_2_sequence.fastq.gz

breseq -r /projectnb/hfsp/IAMM_genomes_from_Luca/D20-160039_assembly.fasta -o /projectnb/hfsp/SNP23/breseq_test/07-citrea_neg_ctrl /projectnb/hfsp/SNP23/results/raw_files/D20-160039-4500T/trimmed_D20-160039_1_sequence.fastq.gz /projectnb/hfsp/SNP23/results/raw_files/D20-160039-4500T/trimmed_D20-160039_2_sequence.fastq.gz

breseq -r /projectnb/hfsp/IAMM_genomes_from_Luca/D20-160042_assembly.fasta -o /projectnb/hfsp/SNP23/breseq_test/08-dss3_neg_ctrl /projectnb/hfsp/SNP23/results/raw_files/D20-160042-4500T/trimmed_D20-160042_1_sequence.fastq.gz /projectnb/hfsp/SNP23/results/raw_files/D20-160042-4500T/trimmed_D20-160042_2_sequence.fastq.gz
