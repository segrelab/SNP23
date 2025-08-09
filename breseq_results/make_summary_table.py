import pandas as pd
import glob
import os
from Bio import SeqIO
import json

def parse_genbank_for_genes_and_lengths(gbk_path):
    """Parses a GenBank file to extract lengths of all genes by locus_tag."""
    genbank_info = []
    try:
        for record in SeqIO.parse(gbk_path, "genbank"):
            for feature in record.features:
                if feature.type == 'CDS':
                    gene_name = feature.qualifiers.get('gene', [None])[0]
                    if gene_name:
                        feature_info = {}
                        feature_info['gene_name'] = gene_name
                        feature_info['length'] = len(feature)
                        feature_info['start_pos'] = feature.location.start
                        feature_info['end_pos'] = feature.location.end
                        # Add the feature info to the list
                        genbank_info.append(feature_info)
        # Convert the list of dictionaries to a DataFrame
        genbank_df = pd.DataFrame(genbank_info)
    except FileNotFoundError:
        print(f"  - WARNING: GenBank file not found at {gbk_path}")
    except Exception as e:
        print(f"  - WARNING: Could not parse GenBank file {gbk_path}. Error: {e}")
    return genbank_df


def get_gene_length(row, length_map):
    """
    Safely gets a gene's length from the pre-computed map.
    Returns None if the mutation is intergenic or the gene is not found.
    """
    # 1. If the mutation is intergenic, return None immediately.
    if 'intergenic' in str(row.get('gene_position', '')):
        return None
    
    gene_name = row.get('gene_name')

    # 2. Check if the gene_name exists and is a key in the map.
    if pd.notna(gene_name) and gene_name in length_map:
        # 3. If it exists, return the length.
        return length_map[gene_name]
    
    # 4. If it's not found, print a warning and return None.
    # Avoid printing warnings for common non-gene entries like '–/–'.
    if pd.notna(gene_name) and gene_name != '–/–':
        print(f"  - WARNING: Gene '{gene_name}' not found in the reference GenBank length map.")
    return None


def analyze_breseq_outputs(base_dir, output_csv_path, meta_file_path, gbk_base_dir, detailed_output_dir):
    """
    Analyzes breseq output files, saves a detailed table with gene lengths,
    and creates a final summary of SNP counts and most frequent genes.
    """
    try:
        meta_df = pd.read_csv(meta_file_path)
        meta_df.set_index('strain_id', inplace=True)
        print("Successfully loaded metadata file.")
    except (FileNotFoundError, KeyError) as e:
        print(f"FATAL: Could not load or process metadata file {meta_file_path}. Error: {e}")
        return

    os.makedirs(detailed_output_dir, exist_ok=True)
    
    search_pattern = os.path.join(base_dir, '*/plasmidsaurus_vs_ref/output/output.tsv')
    file_paths = glob.glob(search_pattern)

    if not file_paths:
        print(f"Warning: No 'output.tsv' files found matching the pattern: {search_pattern}")
        return

    print(f"Found {len(file_paths)} samples to analyze...")
    summary_data = []
    gbk_cache = {}

    for path in file_paths:
        try:
            parts = path.split(os.sep)
            sample_id = parts[-4]
            print(f"\nProcessing sample: {sample_id}")

            df = pd.read_csv(path, sep='\t', on_bad_lines='skip')

            # Pull out the total number of mutations found (the number of rows in the DataFrame)
            snp_count = len(df)

            # Get the counts of each unique value in the 'mutation_category' column
            mutation_type_counts = df['mutation_category'].value_counts().to_dict()

            # # Get the path to the GenBank file and make sure it exists
            # try:
            #     gbk_filename = meta_df.loc[sample_id, 'ref_genome_name']
            #     if pd.notna(gbk_filename):
            #         gbk_path = os.path.join(gbk_base_dir, gbk_filename, f"{gbk_filename}.gbk")
            #     else:
            #         print(f"  - WARNING: No GenBank file listed for sample {sample_id} in metadata.")
            # except KeyError:
            #     print(f"  - WARNING: Sample ID '{sample_id}' not found in metadata file.")

            # # Get the locus gene names and lengths from the GenBank file
            # if gbk_path not in gbk_cache:
            #     print(f"  - Parsing GenBank file: {gbk_path}")
            #     gbk_cache[gbk_path] = parse_genbank_for_genes_and_lengths(gbk_path)
            # genbank_df = gbk_cache[gbk_path]
                
            # if not genbank_df.empty:
            #     # Create the gene length mapping dictionary once for efficiency
            #     gene_length_map = genbank_df.set_index('gene_name')['length'].to_dict()

            #     # Add gene lengths, but only if the mutation is not 'intergenic'
            #     df['gene_length'] = df.apply(
            #         lambda row: get_gene_length(row, gene_length_map),
            #         axis=1
            #     )
            #     detailed_filename = os.path.join(detailed_output_dir, f"{sample_id}_details_with_lengths.csv")
            #     df.to_csv(detailed_filename, index=False)
            #     print(f"  - Saved detailed analysis with gene lengths to: {detailed_filename}")
            # else:
            #     print("  - Skipping gene length analysis due to missing GenBank info.")

            # # --- Calculate SNP density to find the gene with the most mutations per base pair ---
            # gene_with_highest_density = 'N/A'
            # highest_density = 0.0

            # # Proceed only if the gene_length column was successfully added
            # if 'gene_length' in df.columns:
            #     # Create a working copy, filtering for rows with valid gene names and lengths
            #     coding_mutations = df.dropna(subset=['gene_name', 'gene_length']).copy()
                
            #     # Ensure gene_length is a numeric type, converting non-numeric values (like 'N/A') to NaN
            #     coding_mutations['gene_length'] = pd.to_numeric(coding_mutations['gene_length'], errors='coerce')
            #     coding_mutations.dropna(subset=['gene_length'], inplace=True)
                
            #     # Ensure gene lengths are > 0 to avoid division-by-zero errors
            #     coding_mutations = coding_mutations[coding_mutations['gene_length'] > 0]

            #     if not coding_mutations.empty:
            #         # 1. Get SNP counts for each gene
            #         gene_counts = coding_mutations['gene_name'].value_counts()

            #         # 2. Get the unique length for each gene
            #         gene_lengths = coding_mutations.drop_duplicates('gene_name').set_index('gene_name')['gene_length']
                    
            #         # 3. Align counts and lengths to ensure we only divide where both values exist
            #         aligned_counts, aligned_lengths = gene_counts.align(gene_lengths, join='inner')

            #         if not aligned_counts.empty:
            #             # 4. Calculate density (SNPs per base pair)
            #             snp_density = aligned_counts / aligned_lengths
                        
            #             if not snp_density.empty:
            #                 gene_with_highest_density = snp_density.idxmax()
            #                 highest_density = snp_density.max()

            # Get the total number of reads and the number of contigs and the fit_mean for the longest contig from the sample's summary.json file
            summary_json_path = os.path.join(os.path.dirname(path), 'summary.json')
            num_contigs = None
            coverage_of_longest = None

            if os.path.exists(summary_json_path):
                try:
                    with open(summary_json_path, 'r') as f:
                        summary_json = json.load(f)

                    # Get the total number of reads from "reads"
                    total_reads = summary_json.get('reads', {}).get('total_reads', 0)
                    
                    # Safely get the dictionary of references
                    references_dict = summary_json.get('references', {}).get('reference', {})

                    if references_dict:
                        # Get the total number of contigs (references)
                        num_contigs = len(references_dict)

                        # Find the dictionary of the reference with the maximum length
                        # The .get('length', 0) makes it safe if a reference is missing a 'length' key
                        longest_ref_data = max(references_dict.values(), key=lambda ref: ref.get('length', 0))
                        
                        # Get the coverage average for that longest reference
                        coverage_of_longest = longest_ref_data.get('coverage_average')
                    else:
                        total_reads = 0
                        num_contigs = 0
                        coverage_of_longest = None

                except Exception as e:
                    print(f"  - WARNING: Could not read or process summary.json for {sample_id}. Error: {e}")
            else:
                print(f"  - WARNING: No summary.json found for {sample_id}.")

            # Prepare the base summary data for the current sample
            sample_summary = {
                'sample_id': sample_id,
                'total_predicted_mutations': snp_count,
                # 'gene_with_highest_density': gene_with_highest_density,
                # 'mutation_density_mutations_per_kb': highest_density * 1000,
                'total_reads': total_reads,
                'num_contigs': num_contigs,
                'avg_coverage_longest_contig': coverage_of_longest
            }

            # Add the dynamic mutation category counts to the summary
            sample_summary.update(mutation_type_counts)
            
            # Add the sample summary to the list
            summary_data.append(sample_summary)

        except Exception as e:
            print(f"Error processing file {path}: {e}")

    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(output_csv_path, index=False)
        print(f"\nAnalysis complete. Summary saved to: {output_csv_path}")
    else:
        print("\nNo data was processed. Output file will not be created.")

if __name__ == "__main__":
    breseq_base_dir = "/projectnb/hfsp/SNP23/breseq_results/raw_files"
    meta_file = "/projectnb/hfsp/SNP23/breseq_results/merged_metafiles.csv"
    gbk_dir = "/projectnb/hfsp/IAMM_reference_files/prokka_results"
    summary_output_file = "/projectnb/hfsp/SNP23/breseq_results/snp_summary.csv"
    detailed_dir = "/projectnb/hfsp/SNP23/breseq_results/detailed_length_reports"
    
    print("Starting analysis. Make sure you have installed biopython: pip install biopython")
    analyze_breseq_outputs(breseq_base_dir, summary_output_file, meta_file, gbk_dir, detailed_dir)