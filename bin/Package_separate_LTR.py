import click
import os
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


@click.command()
@click.option('--input_file', '-i',  required=True, type=str, help='Path to the input sequence FASTA file.')
@click.option('--output_dir', '-o', default=None, type=str, help='Directory to save output files.')
def process_sequences(input_file, output_dir):
    """
    Process each sequence in the input file to detect and handle LTRs.
    """
    input_dir = os.path.dirname(input_file)
    input_name = os.path.basename(input_file)
    sep_seq_n = 0

    if output_dir is None:
        output_dir = os.path.join(input_dir, f"{input_name}_LTR_separate")
        output_file = f"{input_file}_LTR_separated.fa"
    else:
        output_dir = os.path.join(output_dir, f"{input_name}_LTR_separate")
        output_file = os.path.join(output_dir, f"{input_name}_LTR_separated.fa")

    os.makedirs(output_dir, exist_ok=True)
    sequences = list(SeqIO.parse(input_file, "fasta"))
    processed_sequences = []
    print("LTR separating is running ......")

    for record in sequences:
        # Process each sequence to detect LTR
        LTR_boundary = detect_ltr_for_sequence(record, output_dir)

        # Handle the sequence based on LTR detection
        if LTR_boundary:
            # Split the sequence at LTR boundary
            LTR_seq = record.seq[:LTR_boundary[1]]
            INT_seq = record.seq[LTR_boundary[1]-1:LTR_boundary[2]]

            # Update headers
            new_id_LTR = update_fasta_header(record.id, "_LTR")
            new_id_INT = update_fasta_header(record.id, "_INT")

            # Create new SeqRecord objects
            LTR_record = SeqRecord(Seq(LTR_seq), id=new_id_LTR, description="")
            INT_record = SeqRecord(Seq(INT_seq), id=new_id_INT, description="")

            processed_sequences.extend([LTR_record, INT_record])
            print(f"{record.id} contains 'LTR' and was written into separate sequence file.")
            sep_seq_n += 1
        else:
            # Keep sequence unchanged if LTR was not identified
            processed_sequences.append(record)
    print(f"\n{sep_seq_n} sequences were found to contain 'LTR' and written into separate sequence file."
          f"\nSeparation is finished.")

    # Write processed sequences to a new file
    with open(output_file, "w") as f:
        SeqIO.write(processed_sequences, f, "fasta")


def detect_ltr_for_sequence(record, output_dir):
    """
    Detect LTR in a single sequence using BLASTn.
    """
    record_len = record_len = len(record.seq)
    record_name = record.id.split("#")[0].replace("/", "_")
    sequence_folder = os.path.join(output_dir, record_name)
    os.makedirs(sequence_folder, exist_ok=True)
    sequence_file = os.path.join(sequence_folder, record_name)

    # Write the current sequence to a temporary file
    SeqIO.write(record, sequence_file, "fasta")

    # Prepare BLAST database and output paths
    database_file = os.path.join(sequence_folder, "temp_blast_database")
    makeblastdb_cmd = f"makeblastdb -in {sequence_file} -dbtype nucl -out {database_file}"

    try:
        subprocess.run(makeblastdb_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"\nmakeblastdb encountered an error for sequence {record_name} and returned error code {e.returncode}.\n")
        print(e.stdout)
        print(e.stderr)
        return None

    blast_cmd = f"blastn -query {sequence_file} -db {database_file} " \
                f"-outfmt \"6 qseqid qstart qend sstart send\" " \
                f"-evalue 0.05"  # Set a higher evalue for self-blast

    # Execute the command
    result = subprocess.run(blast_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if result.returncode != 0:
        print(f"An error occurred in BLAST: {result.stderr.decode('utf-8')}")
        return None

    blast_out = result.stdout.decode('utf-8').strip()
    if not blast_out:
        return None

    # Process BLAST output
    # Process BLAST output
    columns = ['qseqid', 'qstart', 'qend', 'sstart', 'send']
    df = pd.DataFrame([x.split('\t') for x in blast_out.split('\n')], columns=columns)

    if df.empty:
        return None

    # Convert columns to integer
    df[['qstart', 'qend', 'sstart', 'send']] = df[['qstart', 'qend', 'sstart', 'send']].astype(int)

    df_LTR = df[(df['qend'] - df['qstart'] >= 150) &
                (df['qstart'] != df['sstart']) &
                (df['qend'] != df['send']) &
                (df['sstart'] < df['send']) &
                (df['sstart'] > df['qend']) &
                (df['sstart'] > df['qstart'])].copy()

    if not df_LTR.empty:
        df_LTR["diff"] = df['send'] - df['qstart']
        df_LTR.reset_index(drop=True, inplace=True)

        # Find the row with the largest difference
        LTR_largest = df_LTR.iloc[df_LTR["diff"].idxmax()]

        # Check if the terminal repeat spans the majority of the query sequence. Because the query was extended,
        # assume the maximum redundant extension for left and right side at 2000.
        if abs(LTR_largest['send'] - LTR_largest['qstart']) >= (record_len - 200):
            # Because BLAST uses index starting from 1, modify the start position
            LTR_boundary = [LTR_largest['qstart'] - 1, LTR_largest['qend'], LTR_largest['sstart'], LTR_largest['send']]
            # print(record_len)
            # print(LTR_boundary)
            return LTR_boundary
        else:
            return None
    else:
        return None


def update_fasta_header(header, suffix):
    """
    Update FASTA header with suffix before '#', or add suffix at the end if '#' not present.
    """
    if '#' in header:
        parts = header.split('#')
        return f'{parts[0]}{suffix}#{parts[1]}'
    else:
        return f'{header}{suffix}'


if __name__ == '__main__':
    process_sequences()
