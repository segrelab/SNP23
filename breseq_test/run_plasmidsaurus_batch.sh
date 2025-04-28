#!/bin/bash -l

#$ -P hfsp              # Specify the SCC project name you want to use
#$ -l h_rt=12:00:00     # Specify the hard time limit for the job
#$ -N plasmidsaurus-test   # Give job a name
#$ -j y                 # Merge the error and output streams into a single file

module load miniconda
conda activate /projectnb/hfsp/SNP23/envs/breseq

# HOT1A3
# Commenting this out because I already ran it
# breseq -r /projectnb/hfsp/Plasmidsaurus25/genomes/28108.53.fna -o breseq_test/09-hot1a3-plasmidsaurus /projectnb/hfsp/Plasmidsaurus25/fastqs/TKDHS4_7_S39_R1_001.fastq.gz /projectnb/hfsp/Plasmidsaurus25/fastqs/TKDHS4_7_S39_R2_001.fastq.gz

# HOT1A3 negative control
breseq -r /projectnb/hfsp/Plasmidsaurus25/archive/TKDHS4_polished_results/TKDHS4_7_polish/TKDHS4_7_polished.fasta -o breseq_test/11-hot1a3-plasmidsaurus-neg-cntrl /projectnb/hfsp/Plasmidsaurus25/fastqs/TKDHS4_7_S39_R1_001.fastq.gz /projectnb/hfsp/Plasmidsaurus25/fastqs/TKDHS4_7_S39_R2_001.fastq.gz

# HOT5_e6
# Commenting this out because I already ran it
# breseq -r /projectnb/hfsp/Plasmidsaurus25/genomes/E06.fna -o breseq_test/10-hot5_e6-plasmidsaurus /projectnb/hfsp/Plasmidsaurus25/fastqs/TKDHS4_25_S57_R1_001.fastq.gz /projectnb/hfsp/Plasmidsaurus25/fastqs/TKDHS4_25_S57_R2_001.fastq.gz

# HOT5_e6 negative control
breseq -r /projectnb/hfsp/Plasmidsaurus25/archive/TKDHS4_polished_results/TKDHS4_25_polish/TKDHS4_25_polished.fasta -o breseq_test/12-hot5_e6-plasmidsaurus-neg-cntrl /projectnb/hfsp/Plasmidsaurus25/fastqs/TKDHS4_25_S57_R1_001.fastq.gz /projectnb/hfsp/Plasmidsaurus25/fastqs/TKDHS4_25_S57_R2_001.fastq.gz
