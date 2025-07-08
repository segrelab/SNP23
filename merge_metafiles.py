import pandas as pd
import sys
import os

def merge_genome_files(plasmidsaurus_path, iamm_path, output_path):
    """
    Merges two CSV files based on specified columns.

    Args:
        plasmidsaurus_path (str): Path to plasmidsaurus_metafile.csv.
        iamm_path (str): Path to iamm_reference_genomes_with_paths.csv.
        output_path (str): Path for the output merged CSV file.
    """
    # --- 1. Input Validation ---
    if not os.path.exists(plasmidsaurus_path):
        print(f"Error: File not found at '{plasmidsaurus_path}'", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(iamm_path):
        print(f"Error: File not found at '{iamm_path}'", file=sys.stderr)
        sys.exit(1)

    # --- 2. Read CSV files into pandas DataFrames ---
    print("Reading input files...")
    try:
        plasmidsaurus_df = pd.read_csv(plasmidsaurus_path)
        iamm_df = pd.read_csv(iamm_path)
    except Exception as e:
        print(f"Error reading CSV files: {e}", file=sys.stderr)
        sys.exit(1)

    # --- 3. Check for required columns ---
    if 'strain_id' not in plasmidsaurus_df.columns:
        print(f"Error: 'strain_id' column not found in {plasmidsaurus_path}", file=sys.stderr)
        sys.exit(1)
    if 'strain_id' not in iamm_df.columns:
        print(f"Error: 'strain_id' column not found in {iamm_path}", file=sys.stderr)
        sys.exit(1)

    # --- 4. Perform the merge ---
    print("Merging dataframes on 'strain_id'...")
    # Using a 'left' merge to keep all rows from the plasmidsaurus file
    merged_df = pd.merge(
        plasmidsaurus_df,
        iamm_df,
        left_on='strain_id',
        right_on='strain_id',
        how='left'
    )

    # Remove rows where 'plasmidsaurus_id' is NaN
    merged_df = merged_df[merged_df['plasmidsaurus_id'].notna()]

    # --- 5. Save the result ---
    merged_df.to_csv(output_path, index=False)
    print(f"\nSuccessfully merged files.")
    print(f"Output saved to: {output_path}")


if __name__ == "__main__":
    # Define file paths
    plasmidsaurus_file = "/projectnb/hfsp/SNP23/plasmidsaurus_metafile.csv"
    iamm_file = "/projectnb/hfsp/iamm-collection/iamm_references.csv"
    output_file = "/projectnb/hfsp/SNP23/breseq_results/merged_metafiles.csv"
    
    merge_genome_files(plasmidsaurus_file, iamm_file, output_file)