#!/bin/bash

echo "strain_id,ref_genome,num_contigs" > contig_counts.csv

while IFS=, read -r sample_number sample_name species_name strain_id in_IAMM accession_number source dna_prep ref_genome _
do
    # Skip header
    [[ $strain_id == "strain_id" ]] && continue
    [[ -z "$ref_genome" ]] && continue

    # Change the ref genome file name from ending in .fna to .gb
    ref_genome=${ref_genome%.fna}.gbk

    # Add the path to the ref genome file
    ref_genome="/projectnb/hfsp/IAMM_reference_files/gbk/$ref_genome"

    if [[ -f "$ref_genome" ]]; then
        count=$(grep -c "^LOCUS" "$ref_genome")
        echo "$strain_id,$ref_genome,$count" >> contig_counts.csv
    fi
done < metafile.csv