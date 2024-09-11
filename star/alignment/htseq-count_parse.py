#!/usr/bin/env python

def extract_gene_counts(slurm_file, output_files):
    htseq_count = 0  # Counter for tracking the number of HTSeq-count outputs
    gene_counts = []  # List to hold the gene count lines

    with open(slurm_file, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if line.startswith("ENSMUSG"):
            # Start collecting gene counts
            gene_counts.append(line)

            # Continue collecting until we hit the "__no_feature" line
            j = i + 1
            while not lines[j].startswith("__no_feature"):
                gene_counts.append(lines[j])
                j += 1

            # Write to the appropriate output file
            with open(output_files[htseq_count], 'w') as output_file:
                output_file.writelines(gene_counts)

            # Reset and move to the next HTSeq-count section
            htseq_count += 1
            gene_counts = []  # Clear the list for the next set of gene counts

            # Stop once all four HTSeq-count outputs are processed
            if htseq_count == 4:
                break


# File paths
slurm_file = 'slurm-16039341.out'
output_files = [
    'stranded_3_2B_control_S3_L008.genecount',
    'reverse_3_2B_control_S3_L008.genecount',
    'stranded_32_4G_both_S23_L008.genecount',
    'reverse_32_4G_both_S23_L008.genecount'
]

# Call the function to extract gene counts
extract_gene_counts(slurm_file, output_files)
