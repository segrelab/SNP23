# Run the genbank with one trimmed file (this is what I think Daniel did)
# Ran 3:10 PM - no later than 3:19 PM
# The results are not the exact same as what Daniel got, maybe he used a different read file
breseq -r /projectnb/hfsp/Library_genomes/gbk/2687453488.gbk -o breseq_test/01-HOT1A3_genbank_read1 results/raw_files/D20-160028-4500T/trimmed_D20-160028_1_sequence.fastq.gz

# Then the other trimmed read file
# The results are very similar, but not the exact same as the one generated above
breseq -r /projectnb/hfsp/Library_genomes/gbk/2687453488.gbk -o breseq_test/02-HOT1A3_genbank_read2 results/raw_files/D20-160028-4500T/trimmed_D20-160028_2_sequence.fastq.gz

# Try both read files
# Ran fine with one reference genome and two read files
breseq -r /projectnb/hfsp/Library_genomes/gbk/2687453488.gbk -o breseq_test/03-HOT1A3_genbank_both_reads results/raw_files/D20-160028-4500T/trimmed_D20-160028_1_sequence.fastq.gz results/raw_files/D20-160028-4500T/trimmed_D20-160028_2_sequence.fastq.gz

# Both read files with the fast reference genome
breseq -r /projectnb/hfsp/IAMM_genomes_from_Luca/2687453488.fna -o breseq_test/04-HOT1A3_fast_both_reads results/raw_files/D20-160028-4500T/trimmed_D20-160028_1_sequence.fastq.gz results/raw_files/D20-160028-4500T/trimmed_D20-160028_2_sequence.fastq.gz

# HOT1A3 negative control
# Using the fasta file means I can't get the gene information, but I think that's okay
breseq -r /projectnb/hfsp/IAMM_genomes_from_Luca/D20-160028_assembly.fasta -o breseq_test/05-HOT1A3_neg_cntrl results/raw_files/D20-160028-4500T/trimmed_D20-160028_1_sequence.fastq.gz results/raw_files/D20-160028-4500T/trimmed_D20-160028_2_sequence.fastq.gz
