#!/bin/bash -l
#$ -l h_rt=75:00:00 # Specify the hard time limit for the job
#$ -j y             # Merge the error and output streams into a single file
#$ -o run_all_breseq.out # Specify the output file
#$ -P hfsp          # Specify the SCC project name you want to use
#$ -N breseq  # Give your job a name
#$ -pe omp 8       # Request multiple slots for the Shared Memory application (OpenMP)

module load miniconda
conda activate /projectnb/hfsp/SNP23/envs/breseq

# Set the path to the metafile
input_fi="plasmidsaurus_metafile.csv"
skip_first_row=true # A flag to skip the first row, set to false if the column names are not included in the metafile

# Set a variable of if you want to rerun the analysis
force_rerun=false

# Process each sample in the input file
# These variables must match the order/contents of columns in the input file
while IFS=, read -r strain_id species_name plasmidsaurus_id gbk_ref_genome path_to_gbk pos_cntrl_genome pos_cntrl_genome_path neg_cntrl_genome neg_cntrl_genome_path Notes;do
    # Skip the first row (not a real sample, just a header)
    if $skip_first_row; then
        skip_first_row=false
        continue
    fi

    # Trim any leading/trailing whitespace including special characters
    ref_path_trimmed=$(echo "$path_to_gbk" | xargs | tr -d '\r')
    ref_genome_trimmed=$(echo "$gbk_ref_genome" | xargs | tr -d '\r')

    # Remove any trailing slashes from ref_path_trimmed
    ref_path_trimmed="${ref_path_trimmed%/}"

    # If the reference genome file name is empty, skip it
    if [ -z "$ref_genome_trimmed" ]; then
        continue
    fi

    # Create a single variable for the reference genome path and file name
    ref="${ref_path_trimmed}/${ref_genome_trimmed}"

    # Check that the reference genome file exists
    if [ ! -f "$ref" ]; then
        echo "Reference genome file not found for $plasmidsaurus_id: $ref"
        continue
    fi

    # Get the batch label, the plasmidsaurus_id before the underscore eg TKDHS4
    batch=$(echo "$plasmidsaurus_id" | cut -d'_' -f1)

    # Get the polished sample ID (not actually sure what this means)
    # By taking the first two parts of the plasmidsaurus_id
    sample_name=$(echo "$plasmidsaurus_id" | cut -d'_' -f1-2)

    # Get the directory to the raw results
    # TODO: Make this a variable that I set at the top of the script
    in_dir=/projectnb/hfsp/Plasmidsaurus25/archive/${batch}_polished_results/${sample_name}_polish

    # Set the directiory for the results
    out_dir="breseq_results/raw_files/${plasmidsaurus_id}"

    # Run the breseq for the current sample
    breseq -r $ref -o $out_dir $in_dir/${plasmidsaurus_id}_R1_001.fastq.gz $in_dir/${plasmidsaurus_id}_R2_001.fastq.gz

    # # Process positive control if present
    #  if [ -n "$pos_cntrl_genome" ] && [ -n "$pos_cntrl_path" ]; then
    #     pos_cntrl_path_trimmed=$(echo "$pos_cntrl_path" | xargs | tr -d '\r')
    #     pos_cntrl_genome_trimmed=$(echo "$pos_cntrl_genome" | xargs | tr -d '\r')
    #     pos_cntrl_ref="${pos_cntrl_path_trimmed}/${pos_cntrl_genome_trimmed}"
    #     pos_cntrl_directory="results/raw_files/$sample_name/controls/positive_control"  # TODO: Make this a variable that I set at the top of the script

    #     mkdir -p "$pos_cntrl_directory"
    #     run_analysis_pipeline "${sample_name}_pos_cntrl" "$label" "$pos_cntrl_ref" "$directory" "$pos_cntrl_directory"
    # fi

    # # Process negative control if present
    # if [ -n "$neg_cntrl_genome" ] && [ -n "$neg_cntrl_path" ]; then
    #     neg_cntrl_path_trimmed=$(echo "$neg_cntrl_path" | xargs | tr -d '\r')
    #     neg_cntrl_genome_trimmed=$(echo "$neg_cntrl_genome" | xargs | tr -d '\r')
    #     neg_cntrl_ref="${neg_cntrl_path_trimmed}/${neg_cntrl_genome_trimmed}"
    #     neg_cntrl_directory="results/raw_files/$sample_name/controls/negative_control"

    #     mkdir -p "$neg_cntrl_directory"
    #     echo "Running analysis pipeline for negative control:"
    #     echo $neg_cntrl_directory
    #     echo $neg_cntrl_ref
    #     run_analysis_pipeline "${sample_name}_neg_cntrl" "$label" "$neg_cntrl_ref" "$directory" "$neg_cntrl_directory"
    # fi

done < "$input_fi"