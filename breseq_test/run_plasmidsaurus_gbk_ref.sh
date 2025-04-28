#!/bin/bash -l

#$ -P hfsp              # Specify the SCC project name you want to use
#$ -l h_rt=12:00:00     # Specify the hard time limit for the job
#$ -N gbk-ref-test   # Give job a name
#$ -j y                 # Merge the error and output streams into a single file
#$ -m e                # Send email when the job ends
#$ -pe omp 8          # Request 8 cores for the job (breseq is a single node, multi-threaded program)

module load miniconda
conda activate /projectnb/hfsp/SNP23/envs/breseq

# HOT1A3
# Commenting this out becuase I don't know what GBK file to use
# breseq -r /projectnb/hfsp/Plasmidsaurus25/genomes/28108.53.fna -o breseq_test/09-hot1a3-plasmidsaurus /projectnb/hfsp/Plasmidsaurus25/fastqs/TKDHS4_7_S39_R1_001.fastq.gz /projectnb/hfsp/Plasmidsaurus25/fastqs/TKDHS4_7_S39_R2_001.fastq.gz

# HOT5_e6
breseq -r /projectnb/hfsp/Library_genomes/gbk/HOT5_E06.gbk -o breseq_test/13-hot5_e6-plasmidsaurus-gbk /projectnb/hfsp/Plasmidsaurus25/fastqs/TKDHS4_25_S57_R1_001.fastq.gz /projectnb/hfsp/Plasmidsaurus25/fastqs/TKDHS4_25_S57_R2_001.fastq.gz
