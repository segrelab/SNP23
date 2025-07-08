#!/bin/bash -l
#$ -l h_rt=12:00:00                 # Specify the hard time limit for the job
#$ -j y                             # Merge the error and output streams into a single file
#$ -o run_de1_cntrl.$JOB_ID.out     # Specify the output file
#$ -P hfsp                          # Specify the SCC project name you want to use
#$ -N de1                           # Give your job a name

module load miniconda
conda activate /projectnb/hfsp/SNP23/envs/breseq

pos_cntrl_ref="/projectnb/hfsp/IAMM_reference_files/gbk/GCA_000020585.gbk"
out_dir="breseq_results/raw_files/de1"
in_dir="/projectnb/hfsp/Plasmidsaurus25/archive/B9NRXP_polished_results/B9NRXP_8_polish"
plasmidsaurus_id="B9NRXP_8_S112"

breseq -r $pos_cntrl_ref -o ${out_dir}/pos_cntrl $in_dir/${plasmidsaurus_id}_R1_001.fastq.gz $in_dir/${plasmidsaurus_id}_R2_001.fastq.gz && gdtools ANNOTATE -o ${out_dir}/pos_cntrl/output/output.tsv -r $pos_cntrl_ref -f TSV ${out_dir}/pos_cntrl/output/output.gd