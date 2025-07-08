import csv
import os
import sys

def find_genome_files(csv_path, search_dir):
    """
    Finds files in search_dir that start with names from a CSV file
    and writes the results to a new CSV file.
    """
    if not os.path.isdir(search_dir):
        print(f"Error: Directory not found at '{search_dir}'", file=sys.stderr)
        sys.exit(1)

    try:
        with open(csv_path, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            data = list(reader)
            fieldnames = reader.fieldnames
    except FileNotFoundError:
        print(f"Error: CSV file not found at '{csv_path}'", file=sys.stderr)
        sys.exit(1)
    except KeyError:
        print(f"Error: 'ref_genome_name' column not found in '{csv_path}'", file=sys.stderr)
        sys.exit(1)

    print(f"Searching for files in: {search_dir}\n")
    
    # To make search faster, list directory contents once
    files_in_dir = os.listdir(search_dir)
    found_count = 0

    for row in data:
        name_to_find = row.get('ref_genome_name')
        found_path = 'NOT_FOUND'
        if name_to_find and name_to_find != 'NONE':
            for filename in files_in_dir:
                if filename.startswith(name_to_find):
                    found_path = os.path.join(search_dir, filename)
                    found_count += 1
                    break # Stop after finding the first match
        row['found_path'] = found_path

    # Define the output file path
    output_fieldnames = fieldnames + ['found_path']
    base, ext = os.path.splitext(csv_path)
    output_csv_path = f"{base}_with_paths{ext}"

    # Write the updated data to a new CSV file
    with open(output_csv_path, 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=output_fieldnames)
        writer.writeheader()
        writer.writerows(data)

    print(f"Found {found_count} matching files.")
    print(f"Results with file paths written to: {output_csv_path}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} <path_to_csv> <directory_to_search>")
        sys.exit(1)

    find_genome_files(sys.argv[1], sys.argv[2])