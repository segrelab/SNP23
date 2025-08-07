#!/bin/bash -l
#$ -l h_rt=12:00:00                     # Specify the hard time limit for the job
#$ -j y                                 # Merge the error and output streams into a single file
#$ -o run_missing_samples.$JOB_ID.out   # Specify the output file
#$ -P hfsp                              # Specify the SCC project name you want to use
#$ -N breseq-man                        # Give your job a name

module load miniconda
conda activate /projectnb/hfsp/SNP23/envs/breseq

# t6c
breseq -r /projectnb/hfsp/IAMM_reference_files/prokka_results/2606217322/2606217322.gbk -o breseq_results/raw_files/t6c/plasmidsaurus_vs_ref /projectnb/hfsp/Plasmidsaurus25/archive/TKDHS4_polished_results/TKDHS4_14_polish/TKDHS4_14_S46_R1_001.fastq.gz /projectnb/hfsp/Plasmidsaurus25/archive/TKDHS4_polished_results/TKDHS4_14_polish/TKDHS4_14_S46_R2_001.fastq.gz
