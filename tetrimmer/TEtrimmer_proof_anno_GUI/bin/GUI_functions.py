from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import subprocess
import pandas as pd
import click
import traceback

special_character = {
        '/': '__',
        '\\': '__',
        ' ': '_',
        '|': '_',
        '<': '_',
        ':': '_',
        '"': '_',
        '*': '_',
        '?': '_'
    }


def sanitize_name(name):

    for char, replacement in special_character.items():
        name = name.replace(char, replacement)
    return name


def separate_sequences(input_file, output_dir):
    """
    Separates input file into single separate FASTA files
    """
    os.makedirs(output_dir, exist_ok=True)

    with open(input_file, 'r') as fasta_file:

        id_list = []
        for record in SeqIO.parse(fasta_file, 'fasta'):

            # Modify record.id to remove special characters
            sanitized_id = sanitize_name(record.id)

            # double check if sanitized_id is unique. If not, modify sanitized_
            if sanitized_id not in id_list:
                id_list.append(sanitized_id)
            else:
                id_list.append(sanitized_id)
                count = id_list.count(sanitized_id)
                sanitized_id = f"{sanitized_id}_{count}"

            # Define output file name
            output_filename = os.path.join(output_dir, f"{sanitized_id}.fa")

            # Convert sequence name to sanitized_id
            record.id = sanitized_id
            record.description = ''

            # Write single FASTA file using sanitized name
            # the record.id is now same as sanitized_id
            with open(output_filename, 'w') as output_file:
                SeqIO.write(record, output_file, 'fasta')


def blast(input_file, genome_file, blast_out_dir, e_value=1e-40,  bed_file=False):

    # Define blast out without header
    blast_out_file = os.path.join(blast_out_dir, f"{os.path.basename(input_file)}_no_header.b")
    # Define blast out with header
    blast_out_file_header = os.path.join(blast_out_dir, f"{os.path.basename(input_file)}_blast.txt")

    # -outfmt 6 use "\t" as deliminator. -outfmt 10 use "," as deliminator
    blast_cmd = [
        "blastn",
        "-query", input_file,
        "-db", genome_file,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand evalue bitscore",
        "-evalue", str(e_value),
        "-out", blast_out_file
    ]

    subprocess.run(blast_cmd, check=True, text=True)

    # Check if the blast hit number is 0
    if os.path.getsize(blast_out_file) == 0:
        return False, False

    # Add header to blast out
    # Define the header
    header = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "sstrand", "evalue", "bitscore"
    ]

    # Read the BLAST output into a DataFrame
    blast_out_file_df = pd.read_csv(blast_out_file, delimiter='\t', header=None, names=header)

    # Save the DataFrame back to the file
    blast_out_file_df.to_csv(blast_out_file_header, sep='\t', index=False)

    try:
        # Run the BLAST command
        subprocess.run(blast_cmd, check=True, capture_output=True, text=True)

    except subprocess.CalledProcessError as e:
        click.echo(f"An error occurred during BLAST: \n {traceback.format_exc()}")
        click.echo(f"\nBLAST failed with error code {e.returncode}")
        click.echo(e.stderr)
        return False, False

    if not bed_file:
        # Define the header
        col_n = ["V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"]

        # Read the BLAST output file into a pandas DataFrame
        try:
            blast_out_teaid_file = f"{blast_out_file}_TEAid.b"
            blast_data_teaid_df = pd.read_csv(blast_out_file, delimiter="\t", header=None, names=col_n)

            # Save the DataFrame back to the file with the header, use "," as deliminator for TE-Aid
            blast_data_teaid_df.to_csv(blast_out_teaid_file, sep=',', index=False)

            return blast_out_teaid_file, blast_out_file_header

        except Exception as e:
            click.echo(f"An error occurred while processing the BLAST output for TE-Aid : \n {traceback.format_exc()}")
            return False, False
    else:
        try:
            blast_out_bed_file = f"{blast_out_file}.bed"

            blast_data_bed_df = pd.read_csv(blast_out_file, delimiter='\t', header=None)

            # Define a function to determine the strand and adjust columns accordingly
            def process_row(row):
                if 'plus' in row[10]:
                    return [row[1], row[8], row[9], row[0], 0, '+']
                else:
                    return [row[1], row[9], row[8], row[0], 0, '-']

            # Apply the function to each row and create a new DataFrame
            blast_data_bed_df_processed = blast_data_bed_df.apply(process_row, axis=1, result_type='expand')

            # Save bed file
            blast_data_bed_df_processed.to_csv(blast_out_bed_file, sep='\t', index=False, header=False)

            return blast_out_bed_file, blast_out_file_header

        except Exception as e:
            click.echo(f"An error occurred while processing the BLAST output for BED file : \n {traceback.format_exc()}")
            return False, False


# Generate bed file based on fasta header
# The fasta header could look like
"""
> scaffold_1:343622-349068(+)
> scaffold_1:346171-347467(+)
> scaffold_1:346385-347665(+)
> scaffold_1:346200-347196(+)
"""
def fasta_header_to_bed(input_file, output_file):
    unique_lines = set()

    with open(input_file, 'r') as fasta, open(output_file, 'w') as bed:
        for line in fasta:
            if line.startswith('>'):
                header = line.strip().lstrip('>')  # Remove '>'
                parts = header.rsplit(':', 2)  # Split from the right side
                scaffold = parts[0]
                range_strand = parts[1].split('(')
                range_part = range_strand[0]
                strand = range_strand[1].split(')')[0]
                start, end = range_part.split('-')
                bed_line = f'{scaffold}\t{start}\t{end}\tTEtrimmer\t0\t{strand}\n'

                # Add the line to the set if it's not already present
                if bed_line not in unique_lines:
                    bed.write(bed_line)
                    unique_lines.add(bed_line)

    return output_file


# Avoid to use bedtools slop and getfasta to make it compatible with Windows system
def extend_bed_regions(bed_file, left_extension, right_extension, chrom_sizes, output_bed):

    with open(bed_file, 'r') as infile, open(output_bed, 'w') as outfile:
        for line in infile:
            parts = line.strip().split()
            chrom, start, end, strand = parts[0], int(parts[1]), int(parts[2]), parts[5]

            # Adjust extensions based on strand
            if strand == '+':
                new_start = max(0, start - left_extension)
                new_end = min(chrom_sizes[chrom], end + right_extension)
            else:  # For anti sense strand
                new_start = max(0, start - right_extension)  # Extend "start" less, as it's the "end"
                new_end = min(chrom_sizes[chrom], end + left_extension)  # Extend "end" more, as it's the "start"

            # Write the extended region to the output BED file
            outfile.write(f"{chrom}\t{new_start}\t{new_end}\tTEtrimmer\t0\t{strand}\n")
    return output_bed


def read_bed(bed_file):

    with open(bed_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            yield {
                'chrom': parts[0],
                'start': int(parts[1]),
                'end': int(parts[2]),
                'strand': parts[5]
            }


def extract_fasta_from_bed(genome_fasta, bed_file, output_fasta):

    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, 'fasta'))
    with open(output_fasta, 'w') as out_fasta:
        for region in read_bed(bed_file):
            chrom = region['chrom']
            start = region['start']
            end = region['end']
            strand = region['strand']
            if chrom in genome:
                sequence = genome[chrom].seq[start:end]
                if strand == '-':
                    sequence = sequence.reverse_complement()  # Reverse complement if on negative strand
                seq_record = SeqRecord(sequence, id=f"{chrom}:{start}-{end}({strand})", description="")
                SeqIO.write(seq_record, out_fasta, 'fasta')
    return output_fasta