import argparse
from Bio import SeqIO
import re
import subprocess
import pysam
import pandas as pd
import concurrent.futures
import csv
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import numpy as np
from matplotlib.lines import Line2D
import os



def get_complement(sequence):
    """
    Calculate the complementary repeat strand for a given repeat sequence.
    
    Args:
    sequence (str): A DNA sequence (string of A, T, C, G, and possibly N for any base).

    Returns:
    str: The complementary DNA sequence.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}  # Mapping from nucleotide to its complement
    return ''.join([complement[base] for base in sequence.upper()])  # Build and return the complementary sequence


def reverse_sequence(sequence):
    """
    Reverse the order of bases in a repeat sequence.
    
    Args:
    sequence (str): A repeat sequence to be reversed.

    Returns:
    str: The reversed DNA sequence.
    """
    return sequence[::-1]  # Return the sequence in reverse order


def process_repeats(repeat):
    """
    Generate all biologically relevant variations of a DNA sequence: original, complement,
    reversed, and reversed complement.
    
    Args:
    repeat (str): A DNA sequence.

    Returns:
    list: A list containing the original sequence, its complement, the reversed sequence,
          and the reversed complement sequence.
    """
    repeat1 = repeat  # Original sequence
    repeat2 = get_complement(repeat1)  # Complementary sequence
    repeat3 = reverse_sequence(repeat1)  # Reversed sequence
    repeat4 = get_complement(reverse_sequence(repeat1))  # Reversed complement sequence
    return [repeat1, repeat2, repeat3, repeat4]  # Return all variants


def extract_sequences(input_fastq_file, output_fasta1_file, output_fasta2_file, suffix1, suffix2, repeat_patterns):
    """
    Extract sequences located between specific repeated patterns from a FASTQ file and save them to a FASTA file.
    
    Parameters:
    input_fastq_file: Path to the input FASTQ file.
    output_fasta1_file: Path to the output the first gRNA file.
    output_fasta2_file: Path to the output the second gRNA  file.
    suffix1 and suffix2: Suffix string to append to the sequence description.
    repeat_patterns: A list containing regular expressions for the repeated sequences to search for.
    """
    repeat_patterns = [
        r'ATCTACAACAGTAGAAATTC',
        r'GAATTTCTACTGTTGTAGAT',
        r'TAGATGTTGTCATCTTTAAG',
        r'CTTAAAGATGACAACATCTA'
    ]
    fasta_line = 0
    no_three_line = 0
    with open(output_fasta1_file, 'w', encoding='UTF-8') as file1_out, open(output_fasta2_file, 'w', encoding='UTF-8') as file2_out:
        # Parse the FASTQ file and read each sequence record
        for seq_record in SeqIO.parse(input_fastq_file, 'fastq'):
            # Find positions of all matching repeat sequences
            repeat_matches = [
                (match.start(), pattern)
                for pattern in repeat_patterns
                for match in re.finditer(pattern, str(seq_record.seq))
            ]
            # Ensure exactly two repeat sequences are found
            if 2<= len(repeat_matches) <= 3:
                fasta_line += 1
                # Sort the matches by their positions (to ensure they are in order of appearance)
                #repeat_matches.sort(key=lambda x: x[0])改
                # Calculate the start and end positions for the sequence extraction
                start_position = repeat_matches[0][0]
                end_position = repeat_matches[1][0]
                start_repeat_length = len(repeat_matches[0][1])
                sequence_extract = seq_record.seq[start_position + start_repeat_length:end_position] + '\n'
                # Construct the header for the FASTA format
                fasta_header = ">" + seq_record.description + suffix1 + "\n"
                # Write the header and the extracted sequence to the file
                file1_out.writelines([fasta_header, str(sequence_extract)])

                try:
                    start_position2 = repeat_matches[1][0]
                    end_position2 = repeat_matches[2][0]
                    start_repeat_length2 = len(repeat_matches[1][1])
                    sequence_extract2 = seq_record.seq[start_position2 + start_repeat_length2:end_position2] + '\n'
                    # Construct the header for the FASTA format
                    fasta_header = ">" + seq_record.description + suffix2 + "\n"
                    # Write the header and the extracted sequence to the file
                    file2_out.writelines([fasta_header, str(sequence_extract2)])
                except Exception as e:
                    no_three_line += 1 
                    continue
        print(f"-----------------fasta_line---------------------: {fasta_line}")
        print(f"----------------------no_three_line- --------------------{no_three_line}")
    return no_three_line * 3 - fasta_line

def build_bowtie2_index(genome_fa, prefix):
    """
    Build a Bowtie2 index for the specified genome file.
    
    Parameters:
    genome_fa (str): The path to the genome fasta file.
    prefix (str): The prefix for the index files to be created.

    The function attempts to create a Bowtie2 index using the provided genome fasta file.
    It will report success or catch and report any errors encountered during the process.
    """
    try:
        # Run the Bowtie2 build command with subprocess.run
        subprocess.run(['bowtie2-build', genome_fa, prefix], check=True)
        print(f"Bowtie2 index for {genome_fa} created successfully with prefix {prefix}.")
    except subprocess.CalledProcessError as e:
        print(f"Error during bowtie2 index building: {e}")


def run_bowtie2(input_fasta_files, index_prefix, output_sam_files, single_or_dual):
    """
    Execute the Bowtie2 alignment tool to align sequences from a FASTA file against a reference genome.

    Parameters:
    input_fasta_files (list): The path to the FASTA file containing sequences to align.
    index_prefix (str): The prefix for the Bowtie2 index files corresponding to the reference genome.
    output_sam_files (list): The path where the SAM file containing the alignment results will be saved.

    This function runs Bowtie2 with specified options, capturing any errors that occur during execution.
    """
    try:
        if single_or_dual >= 0:
            # Execute Bowtie2 alignment with specified parameters
            result = subprocess.run(['bowtie2', '-f', '-a', '-p', '10', '--local', '-x', index_prefix, '-U', input_fasta_files[0], '-S', output_sam_files[0], '-N', '1', '-L22'], check=True)
            # Check the return code to confirm success
            if result.returncode == 0:
                print(f"Bowtie2 alignment completed successfully. Output saved to '{output_sam_files}'.")
        elif single_or_dual < 0:
            # Execute Bowtie2 alignment with specified parameters
            result1 = subprocess.run(['bowtie2', '-f', '-a', '-p', '10', '--local', '-x', index_prefix, '-U', input_fasta_files[0], '-S', output_sam_files[0], '-N', '1', '-L22'], check=True)
            result2 = subprocess.run(['bowtie2', '-f', '-a', '-p', '10', '--local', '-x', index_prefix, '-U', input_fasta_files[1], '-S', output_sam_files[1], '-N', '1', '-L22'], check=True)
            # Check the return code to confirm success
            if result1.returncode == 0:
                print(f"Bowtie2 alignment completed successfully. Output saved to '{output_sam_files[0]}'.")
            if result2.returncode == 0:
                print(f"Bowtie2 alignment completed successfully. Output saved to '{output_sam_files[1]}'.")
    except subprocess.CalledProcessError as e:
        # Print error message if Bowtie2 fails
        print(f"Error during Bowtie2 alignment for '{input_fasta_files}': {e}")


def separate_and_sum(cigar_string):
    """
    Parse CIGAR string to calculate the total length of matches (M) and deletions (D).

    Args:
    cigar_string (str): The CIGAR string from an alignment.

    Returns:
    int: The total length of matches and deletions.
    """
    # Find all occurrences of numbers followed by 'M' or 'D'
    operations = re.findall(r'\d+[MD]', cigar_string)
    # Sum the lengths for 'M' and 'D' operations
    sum_M = sum(int(item[:-1]) for item in operations if item.endswith('M'))
    sum_D = sum(int(item[:-1]) for item in operations if item.endswith('D'))
    return sum_M + sum_D


def process_sam(input_sam_files, output_mapping_files, single_or_dual):
    """
    Process a SAM file to extract reads that meet specific CIGAR string criteria and
    directly write the filtered reads to a file, optimizing memory and performance.

    Args:
    input_sam_files (list): Path to the SAM file.
    output_mapping_files (list): Path to the output file.
    """
    headers = ['Name', 'Flag', 'gRNA_id', 'Location', 'Mapping quality', 'CIGAR string', 'Sequence']
    if single_or_dual >= 0:
        sam_names = []
        sam_flags = []
        sam_reference_names = []
        sam_poses = []
        sam_mapping_qualities = []
        sam_cigarstrings = []
        sam_sequences = []
        with pysam.AlignmentFile(input_sam_files[0], 'r') as samfile:
            for read in samfile:
                if read.mapping_quality == 255 and 22 <= separate_and_sum(read.cigarstring) <= 24:
                    sam_names.append(read.query_name if read.query_name else 'N')
                    sam_flags.append(read.flag if read.flag is not None else 'N')
                    sam_reference_names.append(read.reference_name if read.reference_name is not None else 'N')
                    sam_poses.append(read.pos if read.pos is not None else 'N')
                    sam_sequences.append(read.query_sequence if read.query_sequence else 'N')
                    sam_cigarstrings.append(read.cigarstring if read.cigarstring else 'N')
                    sam_mapping_qualities.append(read.mapping_quality if read.mapping_quality is not None else 'N')
        output = list(zip(sam_names, sam_flags, sam_reference_names, sam_poses, sam_mapping_qualities, sam_cigarstrings, sam_sequences))
        with open(output_mapping_files[0], 'w') as f:
            f.write(','.join(headers) + '\n')  # Write the headers
            for item in output:
                line = ','.join(map(str, item))
                f.write(line + '\n')

    elif single_or_dual < 0:
        sam_names1 = []
        sam_flags1 = []
        sam_reference_names1 = []
        sam_poses1 = []
        sam_mapping_qualities1 = []
        sam_cigarstrings1 = []
        sam_sequences1 = []
        with pysam.AlignmentFile(input_sam_files[0], 'r') as samfile:
            for read in samfile:
                if read.mapping_quality == 255 and 22 <= separate_and_sum(read.cigarstring) <= 24:
                    sam_names1.append(read.query_name if read.query_name else 'N')
                    sam_flags1.append(read.flag if read.flag is not None else 'N')
                    sam_reference_names1.append(read.reference_name if read.reference_name is not None else 'N')
                    sam_poses1.append(read.pos if read.pos is not None else 'N')
                    sam_sequences1.append(read.query_sequence if read.query_sequence else 'N')
                    sam_cigarstrings1.append(read.cigarstring if read.cigarstring else 'N')
                    sam_mapping_qualities1.append(read.mapping_quality if read.mapping_quality is not None else 'N')
        output1 = list(zip(sam_names1, sam_flags1, sam_reference_names1, sam_poses1, sam_mapping_qualities1, sam_cigarstrings1, sam_sequences1))
        with open(output_mapping_files[0], 'w') as f:
            f.write(','.join(headers) + '\n')  # Write the headers
            for item in output1:
                line = ','.join(map(str, item))
                f.write(line + '\n')
        sam_names2 = []
        sam_flags2 = []
        sam_reference_names2 = []
        sam_poses2 = []
        sam_mapping_qualities2 = []
        sam_cigarstrings2 = []
        sam_sequences2 = []
        with pysam.AlignmentFile(input_sam_files[1], 'r') as samfile:
            for read in samfile:
                if read.mapping_quality == 255 and 22 <= separate_and_sum(read.cigarstring) <= 24:
                    sam_names2.append(read.query_name if read.query_name else 'N')
                    sam_flags2.append(read.flag if read.flag is not None else 'N')
                    sam_reference_names2.append(read.reference_name if read.reference_name is not None else 'N')
                    sam_poses2.append(read.pos if read.pos is not None else 'N')
                    sam_sequences2.append(read.query_sequence if read.query_sequence else 'N')
                    sam_cigarstrings2.append(read.cigarstring if read.cigarstring else 'N')
                    sam_mapping_qualities2.append(read.mapping_quality if read.mapping_quality is not None else 'N')
        output2 = list(zip(sam_names2, sam_flags2, sam_reference_names2, sam_poses2, sam_mapping_qualities2, sam_cigarstrings2, sam_sequences2))
        with open(output_mapping_files[1], 'w') as f:
            f.write(','.join(headers) + '\n')  # Write the headers
            for item in output2:
                line = ','.join(map(str, item))
                f.write(line + '\n')


def merge_same_gRNA(input_gRNA_file1, input_gRNA_file2, output_process_file):
    headers = ['Name', 'gRNA_id']
    # Function to read data using pandas for simplicity and efficiency
    def read_data(filename):
        data_dict = {}
        with open(filename, 'r') as file:
            next(file)  # Skip the header in the first file
            for line in file:
                values = line.strip().split(',')
                name = values[0]
                val = values[2].split('***')[0]
                if name in data_dict:
                    data_dict[name].append(val)
                else:
                    data_dict[name] = [val]
        return data_dict

    # Read data from both files
    file1_data = read_data(input_gRNA_file1)
    file2_data = read_data(input_gRNA_file2)
    # Merge the data
    merged_data = {}
    for name, values in file1_data.items():
        if name in file2_data:
            merged_data[name] = '&'.join(values + file2_data[name])
    # Write merged data to a new file
    with open(output_process_file, 'w') as merged_file:
        merged_file.write(','.join(headers) + '\n')
        for name, value in merged_data.items():
            merged_file.write(f"{name},{value}\n")

    print(f"Merged data has been written to {output_process_file}.")


def count_enrichment(mapping_file, count_file):
    """
    Count the occurrences of each gRNA_id in the mapping file and save the counts to another file.

    Args:
    mapping_file (str): The file path to the input file containing gRNA mappings.
    count_file (str): The file path where the count results will be saved.

    The function reads the mapping file, counts occurrences of each gRNA_id, and writes
    these counts to the count_file.
    """
    # Load the mapping data
    gRNA_mapping = pd.read_csv(mapping_file, sep=',', engine='python')
    # Count the occurrences of each gRNA_id
    gRNA_id_counts = gRNA_mapping['gRNA_id'].value_counts(normalize=False)
    # Save the counts
    gRNA_id_counts.to_csv(count_file, sep=',', header=True)


def process_fastq(input_fastq_file, index_prefix, suffix1, suffix2, repeat_patterns):
    """
    Process a FASTQ file through various steps to analyze genetic sequences.
    
    Args:
    input_fastq_file (str): Path to the input FASTQ file.
    index_prefix (str): Prefix for the Bowtie2 index files.
    suffix (str): Suffix to append to sequence descriptions in the FASTA file.
    repeat_patterns (list): List of repeat patterns used in sequence extraction.

    The function automates the following steps:
    - Converts FASTQ to FASTA and extracts specific sequences (currently commented out).
    - Aligns the sequences using Bowtie2.
    - Processes the SAM file to filter and capture relevant data.
    - Counts the enrichment of specific features from the processed data.
    """
    # Define file paths for outputs based on input FASTQ file
    output_fasta1_file = input_fastq_file.replace('.fq', '_g1.fasta')
    output_fasta2_file = input_fastq_file.replace('.fq', '_g2.fasta')
    output_fasta_files = [output_fasta1_file, output_fasta2_file]
    output_sam1_file = input_fastq_file.replace('.fq', '_g1.sam')
    output_sam2_file = input_fastq_file.replace('.fq', '_g2.sam')
    output_sam_files = [output_sam1_file, output_sam2_file]
    processed_txt = input_fastq_file.replace('.fq', '_processed.txt')
    processed_txt1 = input_fastq_file.replace('.fq', '_processed_g1.txt')
    processed_txt2 = input_fastq_file.replace('.fq', '_processed_g2.txt')
    processed_txts = [processed_txt1, processed_txt2]
    count_txt = input_fastq_file.replace('.fq', '_count_enrichment.txt')

    # Uncomment the following line if sequence extraction is required
    single_or_dual = extract_sequences(input_fastq_file, output_fasta1_file, output_fasta2_file, suffix1, suffix2, repeat_patterns)
    # Run Bowtie2 alignment
    run_bowtie2(output_fasta_files, index_prefix, output_sam_files, single_or_dual)
    # Process the generated SAM file
    process_sam(output_sam_files, processed_txts, single_or_dual)
    if single_or_dual < 0:
        merge_same_gRNA(processed_txts[0], processed_txts[1], processed_txt)
        # Count and record the enrichment
        count_enrichment(processed_txt, count_txt)
    elif single_or_dual >= 0:
        count_enrichment(processed_txt1, count_txt)
    return  single_or_dual

def merge_count_enrichments(input_count_enrichment_1, input_count_enrichment_2, output_merge_file, parallel_name):
    """
    Merge two enrichment count files and write the combined results to a new file,
    excluding specific gRNA_ids based on the parallel_name.

    Args:
    input_count_enrichment_1 (str): Path to the first count enrichment file.
    input_count_enrichment_2 (str): Path to the second count enrichment file.
    output_merge_file (str): Path to the output file where merged results will be saved.
    parallel_name (str): The name of the parallel to determine which gRNA_ids to exclude.
    """
    headers = ['gRNA_id', 'reads']
    with open(input_count_enrichment_1, 'r') as f1, open(input_count_enrichment_2, 'r') as f2, open(output_merge_file, 'w') as output:
        next(f1)  # Skip the header in the first file
        next(f2)  # Skip the header in the second file

        # Create dictionaries from the files
        count_dict_1 = {line.strip().split(',')[0]: int(line.strip().split(',')[1]) for line in f1}
        count_dict_2 = {line.strip().split(',')[0]: int(line.strip().split(',')[1]) for line in f2}

        # Combine the keys from both dictionaries
        all_gRNA_ids = set(count_dict_1.keys()) | set(count_dict_2.keys())
        
        # Write the headers to the output file
        output.write(','.join(headers) + '\n')
        
        # Iterate over all gRNA_ids and merge counts
        for gRNA_id in all_gRNA_ids:
            # Merge counts from both files, defaulting to 0 if not found
            total_count = count_dict_1.get(gRNA_id, 0) + count_dict_2.get(gRNA_id, 0)
            output.write(f"{gRNA_id},{total_count}\n")


def gene_count_enrichment(input_gRNA_enrichment, output_gene_enrichment, single_or_dual):
    """
    Aggregate gRNA enrichment counts to gene-level enrichment counts and save to a file.

    Args:
    input_gRNA_enrichment (str): Path to the input file containing gRNA-based enrichment counts.
    output_gene_enrichment (str): Path to the output file where gene-level enrichment counts will be saved.

    The function processes the gRNA enrichment data by extracting the gene ID from each gRNA ID,
    aggregates counts by gene, sorts the results by total counts in descending order, and
    saves the results to the specified output file.
    """
    # Load the gRNA enrichment data
    enrichment = pd.read_csv(input_gRNA_enrichment, sep = ',')
    def generate_gene_id(x):
        parts = x.split('-')
        gene_1 = f"{parts[0]}&{x.split('&')[1].split('-')[0]}"
        gene_2 = f"{x.split('&')[1].split('-')[0]}&{parts[0]}"
        return gene_1 if gene_1 < gene_2 else gene_2
    if single_or_dual >= 0:
        # Extract gene ID from gRNA ID (Gene ID is before the first '-')
        enrichment['gene_id'] = enrichment['gRNA_id'].str.split('-').str[0]
        print(f"single_or_dual1: {single_or_dual}")
    elif single_or_dual < 0:
        enrichment['gene_id'] = enrichment['gRNA_id'].apply(generate_gene_id)
        print(f"single_or_dual2: {single_or_dual}")
    # Remove the gRNA_id column as it is no longer needed
    enrichment = enrichment.drop(columns=['gRNA_id'])
    # Reorder columns so that 'gene_id' is first
    cols = ['gene_id'] + [col for col in enrichment if col != 'gene_id']
    enrichment = enrichment[cols]
    # Group by gene_id, summing all other numeric columns, and reset index for later processing
    enrichment_grouped = enrichment.groupby('gene_id').sum().reset_index()
    # Sort the grouped data by 'reads' in descending order to see the most enriched genes on top
    enrichment_grouped = enrichment_grouped.sort_values(by='reads', ascending=False)
    # Save the gene-level enrichment data to a CSV file
    enrichment_grouped.to_csv(output_gene_enrichment, index=False)


def normalize_counts(input_enrichment_file, output_normalization_file):
    """
    Normalize the read counts of genes in a file and save the normalized counts to another CSV file.

    Args:
    input_enrichment_file (str): Path to the input file containing gene ids and their read counts.
    output_normalization_file (str): Path to the output file where normalized read counts will be saved.

    The function reads the input file, calculates the total read count, and then computes the
    normalized count for each gene by dividing its read count by the total count. The results
    are then saved to the output file.
    """
    # Load data from the input file
    df = pd.read_csv(input_enrichment_file, sep=',')
    # Calculate the total count of reads across all genes
    total_count = df.iloc[:, 1].sum()
    # Normalize the read counts for each gene
    df.iloc[:, 1] = df.iloc[:, 1] / total_count  # Divide each gene's counts by the total to get the proportion
    # Update the column names to reflect the new data
    df.columns = ['gene_id', 'normalized_reads']
    # Save the normalized data to a CSV file
    df.to_csv(output_normalization_file, sep=',', index=False)


def find_common_prefix(str_list):
    """
    Find the longest common prefix among a list of strings.

    Args:
    str_list (list of str): The list of strings to evaluate.

    Returns:
    str: The longest common prefix shared by all the strings in the list. If the list is empty
         or there is no common prefix, returns an empty string.
    """
    if not str_list:
        return "Find_common_prefix (str_list): The list is empty"

    # Start with the first string in the list as the initial prefix
    prefix = str_list[0]
    
    # Iterate through the list of strings
    for s in str_list:
        # Reduce the prefix until it matches the beginning of s
        while s[:len(prefix)] != prefix and prefix:
            prefix = prefix[:-1]  # Remove the last character of the prefix
            if not prefix:  # If the prefix is empty, return it immediately
                return "Find_common_prefix (str_list): The prefix is empty"
    
    return prefix  # Return the found common prefix


def process_gene(input_gRNA_enrichment, parallel_name, single_or_dual):
    """
    Process gene enrichment by calculating merged gene counts and normalizing them.

    Args:
    input_gRNA_enrichment (str): Path to the input file containing merged raw counts.
    parallel_name (str): Identifier for the parallel experiment, used to generate output file names.

    This function performs the following steps:
    1. Calculate gene-level enrichment counts from merged output data.
    2. Normalize these gene counts to create a normalized gene count file.
    """
    # Generate file names based on the parallel experiment name
    output_gene_enrichment = f"{parallel_name}_merged_gene_enrichment.txt"
    normalized_file = f"{parallel_name}_merged_gene_normalization.txt"

    gene_count_enrichment(input_gRNA_enrichment, output_gene_enrichment, single_or_dual)
    normalize_counts(output_gene_enrichment, normalized_file)


def merge_and_average_counts(input_normalized_files, output_final_file):
    """
    Merge and average normalized counts from multiple files and write the results, including a specified
    gene entry, to an output file.

    Args:
    input_normalized_files (list): List of file paths containing normalized gene read counts.
    output_final_file (str): Path to the output file where the merged and averaged counts will be saved.
    """
    headers = ['gene_id', 'reads_normalization']
    all_counts = []  # List to store counts data from each file
    gene_sets = []   # List to store sets of gene_ids from each file
    
    # Read each file and store data
    for file_path in input_normalized_files:
        current_file_counts = {}
        with open(file_path, 'r') as file:
            next(file)   # Skip header
            for line in file:
                gene_id, norm_count = line.strip().split(',')
                current_file_counts[gene_id] = float(norm_count)
        all_counts.append(current_file_counts)
        gene_sets.append(set(current_file_counts.keys()))

    # Find common genes across all files
    common_genes = set.intersection(*gene_sets)
    
    # Prepare data structure for combined and average counts
    combined_data = {gene: [] for gene in common_genes}
    for gene in common_genes:
        total_count = 0
        for file_counts in all_counts:
            total_count += file_counts[gene]
            combined_data[gene].append(file_counts[gene])
        average_count = total_count / len(input_normalized_files)
        combined_data[gene].append(average_count)
    
    # Write the combined and averaged data to the output file
    with open(output_final_file, 'w') as output:
        header_line = ['gene_id'] + [f'file_{i+1}_norm' for i in range(len(input_normalized_files))] + ['average_norm']
        output.write(','.join(header_line) + '\n')
        for gene, counts in combined_data.items():
            output.write(f"{gene},{','.join(map(str, counts))}\n")


def analyze_data(input_parallel1_file, input_parallel2_file, output_log2fc_file):
    """
    Analyze gene expression data from two parallel experiments to compute log2 fold changes and
    statistical significance.

    Args:
    input_parallel1_file (str): Path to the file containing normalized counts for parallel 1.
    input_parallel2_file (str): Path to the file containing normalized counts for parallel 2.
    output_log2fc_file (str): Path to the output file where the results will be saved.

    The function performs a t-test for each common gene between the two datasets, adjusts p-values
    for multiple testing, calculates log2 fold changes, and saves the results in the specified file.
    """
    def load_data(file_path):
        """Helper function to load data from a file into a dictionary."""
        data = {}
        with open(file_path, 'r') as file:
            csv_reader = csv.DictReader(file)
            for row in csv_reader:
                gene_id = row['gene_id']
                norms = [float(value) for key, value in row.items() if 'norm' in key]
                data[gene_id] = norms
        return data

    # Load data from both files
    data1 = load_data(input_parallel1_file)
    data2 = load_data(input_parallel2_file)

    # Determine the common genes in both datasets
    common_genes = set(data1.keys()) & set(data2.keys())

    # Analyze data: t-test for each common gene
    results = []
    for gene in common_genes:
        norms1 = data1[gene]
        norms2 = data2[gene]
        t_stat, p_value = stats.ttest_ind(norms1, norms2)
        avg_norm1 = np.mean(norms1)
        avg_norm2 = np.mean(norms2)
        results.append((gene, t_stat, p_value, avg_norm1, avg_norm2))

    # Apply multiple test correction to adjust p-values
    p_values = [result[2] for result in results]
    reject, p_values_corrected, _, _ = multipletests(p_values, method='fdr_bh')

    # Compute log2 fold changes and prepare final results
    final_results = []
    for i, (gene, t_stat, p_value, avg_norm1, avg_norm2) in enumerate(results):
        log2fc = np.log2(avg_norm1 / avg_norm2) if avg_norm2 != 0 else float('inf')
        log10p = -np.log10(p_value) if p_value > 0 else float('inf')
        final_results.append((gene, log2fc, t_stat, p_value, p_values_corrected[i], log10p))

    # Write results to the output file
    with open(output_log2fc_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Gene ID', 'Log2fc', 't-statistic', 'p-value', 'FDR', 'Log10p'])
        writer.writerows(final_results)


def plot_data(input_log2fc_file):
    """
    Read data from a file, plot a volcano plot, and save the plot to a file.

    Args:
    input_log2fc_file (str): Path to the input file containing gene data with Log2fc and Log10p values.

    The function generates a volcano plot to visualize the statistical significance and
    log2 fold change of gene expression, and saves the plot as a PNG file.
    """
    # Read file to get data
    result = pd.read_csv(input_log2fc_file)
    max_value = result[result["Log10p"] != 10]["Log10p"].replace(np.inf, np.nan).max()
    result["Log10p"] = result["Log10p"].replace(np.inf, max_value)

    # Plot settings
    plt.figure(figsize=(10, 7))  # Set figure size for the plot
    FDR_threshold = 0.05         # Significance threshold for False Discovery Rate (FDR)
    x_threshold = 0              # Threshold for Log2fc to consider a gene significant

    # Define groups based on thresholds for coloring the scatter plot
    result['group'] = 'Not significant'
    result.loc[(result['FDR'] < FDR_threshold) & (result['Log2fc'] < x_threshold), 'group'] = 'Down-regulated genes'
    result.loc[(result['FDR'] < FDR_threshold) & (result['Log2fc'] > x_threshold), 'group'] = 'Up-regulated genes'

    # Define color palette for the groups
    palette = {'Down-regulated genes': '#ED723F', 'Not significant': '#564232', 'Up-regulated genes': '#5F479A'}

    # Plot all genes with scatter plot
    ax = sns.scatterplot(x='Log2fc', y='Log10p', hue='group', data=result, palette=palette, edgecolor=None, s=5)
    # Remove default legend and create a custom one
    ax.legend_.remove()
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='Down-regulated genes', markersize=10, markerfacecolor='#ED723F'),
        Line2D([0], [0], marker='o', color='w', label='Not significant', markersize=10, markerfacecolor='#564232'),
        Line2D([0], [0], marker='o', color='w', label='Up-regulated genes', markersize=10, markerfacecolor='#5F479A')
    ]
    plt.legend(handles=legend_elements, loc='best')

    # Style the plot with customized appearance
    plt.xlim(float(result['Log2fc'].min()) - 0.5, float(result['Log2fc'].max()) + 0.5)  # Set x-axis limits
    plt.ylim(float(result['Log10p'].min()) - 0.5, float(result['Log10p'].max()) + 0.5)  # Set y-axis limits
    mpl.rcParams['font.weight'] = 'bold'  # Set global font weight to bold
    plt.xlabel('Gene Fitness', fontproperties='Arial', fontsize=18, fontweight='bold')  # X-axis label
    plt.ylabel('-Log$_{10}P$value', fontproperties='Arial', fontsize=18)  # Y-axis label


    # Style the plot borders and tick labels
    for spine in ['bottom', 'left', 'top', 'right']:
        plt.gca().spines[spine].set_linewidth(2)  # Set border line width
    plt.tick_params(width=1.5, labelsize=15)  # Adjust tick width and label size
    for tick_label in plt.gca().get_xticklabels() + plt.gca().get_yticklabels():
        tick_label.set_fontname('Arial')  # Set font style
        tick_label.set_weight('bold')    # Set font weight

    # Save and close the plot
    output_plot_file = input_log2fc_file.split('.')[0] + '.png'
    plt.savefig(output_plot_file, dpi=750, format='png', bbox_inches='tight')  # Save the plot
    plt.close()


def main():
    # Create an argument parser object
    parser = argparse.ArgumentParser(description="Complete workflow including sequence extraction, "
                                                 "bowtie2 index building, alignment, SAM file processing, "
                                                 "and handling repeat sequences.")
    
    # Required and optional arguments
    parser.add_argument("--input_fastq1", required=True, 
                        help="First input FASTQ file for single-end or first paired-end sequencing data.")
    parser.add_argument("--input_fastq2", 
                        help="Second input FASTQ file for paired-end sequencing data. Optional.")
    parser.add_argument("--genome_fa", required=True, 
                        help="Genome FASTA file for Bowtie2 index building.")
    parser.add_argument("--repeat", required=True, 
                        help="The primary DNA repeat sequence to be analyzed.")
    parser.add_argument("--suffix1", default=" I", 
                        help="Suffix to append to sequence IDs in the FASTA file, default is ' I'.")
    parser.add_argument("--suffix2", default=" II", 
                        help="Suffix to append to sequence IDs in the FASTA file, default is ' I'.")
    parser.add_argument("--samples", nargs='+', 
                        help="List of sample identifiers for processing multiple samples.")
    parser.add_argument("--samples_control", nargs='+', 
                        help="List of sample identifiers for the control group in comparative studies.")
    parser.add_argument("--samples_experiment", nargs='+', 
                        help="List of sample identifiers for the experimental group in comparative studies.")

    args = parser.parse_args()

    # Initialize variables for processing
    fastq1 = args.input_fastq1
    if args.input_fastq2:
        fastq2 = args.input_fastq2
    input_fastq_file = process_repeats(args.repeat)
    index_prefix = args.genome_fa.replace('.fasta', '')
    suffix1 = args.suffix1
    suffix2 = args.suffix2
    repeat_patterns = process_repeats(args.repeat)
    samples = args.samples
    samples_control = args.samples_control
    samples_experiment = args.samples_experiment


    # Build bowtie2 index
    build_bowtie2_index(args.genome_fa, index_prefix)


    # Concurrent processing of FASTQ files using ThreadPoolExecutor
    with concurrent.futures.ThreadPoolExecutor() as executor:
        # Start a task for each FASTQ file
        future_to_fastq = {
            executor.submit(process_fastq, fastq1, index_prefix, suffix1, suffix2, repeat_patterns): fastq1
        }
        if args.input_fastq2:
            future_to_fastq[executor.submit(process_fastq, fastq2, index_prefix, suffix1, suffix2, repeat_patterns)] = fastq2

        # Wait for all tasks to be completed
        for future in concurrent.futures.as_completed(future_to_fastq):
            fastq = future_to_fastq[future]
            try:
                single_or_dual = future.result()  # Obtain the function result and check whether there are any exceptions
            except Exception as exc:
                print(f'{fastq} generated an exception: {exc}')

        print(f"single_or_dual: {single_or_dual}")
        # Execute the merge count file after all tasks are completed
        count_txt1 = fastq1.replace('.fq', '_count_enrichment.txt')
        parallel_name = fastq1.rsplit("_")[0]
        merged_output_file = parallel_name + "_merged_count_enrichment.txt"
        
        # If a second FASTQ file is provided, try merging the count files
        if args.input_fastq2:
            count_txt2 = fastq2.replace('.fq', '_count_enrichment.txt')
            merge_count_enrichments(count_txt1, count_txt2, merged_output_file, parallel_name)
        else:
            # If there is only one file, simply rename it to merge output file
            os.rename(count_txt1, merged_output_file, parallel_name)


    # Concurrent processing of gRNA enrichment files using ThreadPoolExecutor
    with concurrent.futures.ThreadPoolExecutor() as executor:
        # Start a task for each gRNA enrichment file
        future_to_gene = {
            executor.submit(process_gene, merged_output_file, parallel_name, single_or_dual): merged_output_file
        }

        # Wait for all tasks to be completed
        for future in concurrent.futures.as_completed(future_to_gene):
            gene = future_to_gene[future]
            try:
                future.result()  # Obtain the function result and check whether there are any exceptions
            except Exception as exc:
                print(f'{gene} generated an gene exception: {exc}')

        normalized_files = []
        if parallel_name in samples_control:
            parallels = samples_control
        elif parallel_name in samples_experiment:
            parallels = samples_experiment 
        for parallel in parallels:
                normalized_file = parallel + "_merged_gene_normalization.txt"
                normalized_files.append(normalized_file)
        sample_name = find_common_prefix(parallels)
        merged_average_file = sample_name + '_final_average_counts.txt'


        with concurrent.futures.ThreadPoolExecutor() as executor:
            future_to_final = {
                executor.submit(merge_and_average_counts, normalized_files, merged_average_file): normalized_files
            }

            # Wait for all tasks to be completed
            for future in concurrent.futures.as_completed(future_to_final):
                final = future_to_final[future]
                try:
                    future.result()  # Obtain the function result and check whether there are any exceptions
                except Exception as exc:
                    print(f'{final} generated an gene exception: {exc}')

            sample_control_name_fore = find_common_prefix(samples_control)
            sample_experiment_name_fore = find_common_prefix(samples_experiment)
            sample_control_name,  sample_experiment_name = sample_control_name_fore + '_final_average_counts.txt', sample_experiment_name_fore + '_final_average_counts.txt'
            log2fc_file = sample_experiment_name_fore + "_vs_" + sample_control_name_fore + "_" + "log2fc.txt"
            analyze_data(sample_experiment_name, sample_control_name, log2fc_file)

            plot_data(log2fc_file)

if __name__ == "__main__":
    main()