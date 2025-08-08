#!/bin/bash -l
#$ -l h_rt=12:00:00                     # Specify the hard time limit for the job
#$ -j y                                 # Merge the error and output streams into a single file
#$ -o run_missing_samples.$JOB_ID.out   # Specify the output file
#$ -P hfsp                              # Specify the SCC project name you want to use
#$ -N breseq-man                        # Give your job a name

module load miniconda
conda activate /projectnb/hfsp/SNP23/envs/breseq

# hot5_f3
echo "================================================================================"
echo "Processing hot5_f3..."
echo "================================================================================"
breseq -r /projectnb/hfsp/IAMM_reference_files/prokka_results/HOT5_F03/HOT5_F03.gbk -o breseq_results/raw_files/hot5_f3/plasmidsaurus_vs_ref /projectnb/hfsp/Plasmidsaurus25/archive/NDPXLJ_polished_results/NDPXLJ_6_polish/NDPXLJ_6_R1_001.fastq.gz /projectnb/hfsp/Plasmidsaurus25/archive/NDPXLJ_polished_results/NDPXLJ_6_polish/NDPXLJ_6_R2_001.fastq.gz

# ksed
echo "================================================================================"
echo "Processing ksed..."
echo "================================================================================"
breseq -r /projectnb/hfsp/IAMM_reference_files/prokka_results/644736380/644736380.gbk -o breseq_results/raw_files/ksed/plasmidsaurus_vs_ref /projectnb/hfsp/Plasmidsaurus25/archive/TKDHS4_polished_results/TKDHS4_27_polish/TKDHS4_27_R1_001.fastq.gz /projectnb/hfsp/Plasmidsaurus25/archive/TKDHS4_polished_results/TKDHS4_27_polish/TKDHS4_27_R2_001.fastq.gz

# roseo
echo "================================================================================"
echo "Processing roseo..."
echo "================================================================================"
breseq -r /projectnb/hfsp/IAMM_reference_files/prokka_results/Roseovarius_mucosus_SMR3/Roseovarius_mucosus_SMR3.gbk -o breseq_results/raw_files/roseo/plasmidsaurus_vs_ref /projectnb/hfsp/Plasmidsaurus25/archive/NDPXLJ_polished_results/NDPXLJ_7_polish/NDPXLJ_7_R1_001.fastq.gz /projectnb/hfsp/Plasmidsaurus25/archive/NDPXLJ_polished_results/NDPXLJ_7_polish/NDPXLJ_7_R2_001.fastq.gz

# sedthio
echo "================================================================================"
echo "Processing sedthio..."
echo "================================================================================"
breseq -r /projectnb/hfsp/IAMM_reference_files/prokka_results/2609459601/2609459601.gbk -o breseq_results/raw_files/sedthio/plasmidsaurus_vs_ref /projectnb/hfsp/Plasmidsaurus25/archive/6LTRG8_polished_results/6LTRG8_2_polish/6LTRG8_2_R1_001.fastq.gz /projectnb/hfsp/Plasmidsaurus25/archive/6LTRG8_polished_results/6LTRG8_2_polish/6LTRG8_2_R2_001.fastq.gz

# t6c
echo "================================================================================"
echo "Processing t6c..."
echo "================================================================================"
breseq -r /projectnb/hfsp/IAMM_reference_files/prokka_results/2606217322/2606217322.gbk -o breseq_results/raw_files/t6c/plasmidsaurus_vs_ref /projectnb/hfsp/Plasmidsaurus25/archive/TKDHS4_polished_results/TKDHS4_14_polish/TKDHS4_14_S46_R1_001.fastq.gz /projectnb/hfsp/Plasmidsaurus25/archive/TKDHS4_polished_results/TKDHS4_14_polish/TKDHS4_14_S46_R2_001.fastq.gz
