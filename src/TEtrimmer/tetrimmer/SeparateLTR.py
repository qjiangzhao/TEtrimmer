import logging
import os
import subprocess

import click
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


@click.command()
@click.option(
    '--input_file',
    '-i',
    required=True,
    type=str,
    help='Path to the input sequence FASTA file.',
)
@click.option(
    '--output_dir', '-o', default=None, type=str, help='Directory to save output files.'
)
def process_sequences(input_file, output_dir):
    """
    Process each sequence in the input file to detect and handle LTRs.
    """
    input_dir = os.path.dirname(input_file)
    input_name = os.path.basename(input_file)
    sep_seq_n = 0

    if output_dir is None:
        output_dir = os.path.join(input_dir, f'{input_name}_LTR_separate')
        output_file = f'{input_file}_LTR_separated.fa'
    else:
        output_dir = os.path.join(output_dir, f'{input_name}_LTR_separate')
        output_file = os.path.join(output_dir, f'{input_name}_LTR_separated.fa')

    os.makedirs(output_dir, exist_ok=True)
    sequences = list(SeqIO.parse(input_file, 'fasta'))
    processed_sequences = []
    logging.info('LTR separation is running ......')

    for record in sequences:
        # Process each sequence to detect LTR
        LTR_boundary = detect_ltr_for_sequence(record, output_dir)

        # Handle the sequence based on LTR detection
        if LTR_boundary:
            # Split the sequence at LTR boundary
            left_LTR_seq = record.seq[LTR_boundary[0] : LTR_boundary[1]]
            INT_seq = record.seq[LTR_boundary[1] : LTR_boundary[2]]
            right_LTR_seq = record.seq[LTR_boundary[2] : LTR_boundary[3]]

            # Count 'N' characters in each side LTR sequence
            left_LTR_N_count = left_LTR_seq.count('N')
            right_LTR_N_count = right_LTR_seq.count('N')

            # Choose the LTR sequence with fewer 'N's
            if left_LTR_N_count <= right_LTR_N_count:
                LTR_seq = left_LTR_seq
            else:
                LTR_seq = right_LTR_seq

            # Update headers
            new_id_LTR = update_fasta_header(record.id, '_LTR')
            new_id_INT = update_fasta_header(record.id, '_INT')

            # Create new SeqRecord objects
            LTR_record = SeqRecord(Seq(LTR_seq), id=new_id_LTR, description='')
            INT_record = SeqRecord(Seq(INT_seq), id=new_id_INT, description='')

            processed_sequences.extend([LTR_record, INT_record])
            logging.info(f"{record.id} contains 'LTR' and was separated.")
            sep_seq_n += 1
        else:
            # Keep sequence unchanged if LTR was not identified
            processed_sequences.append(record)
    logging.info(
        f"\n{sep_seq_n} sequences were found to contain 'LTR' and written into separate sequence file."
        f'\nSeparation is finished.'
    )

    # Write processed sequences to a new file
    with open(output_file, 'w') as f:
        SeqIO.write(processed_sequences, f, 'fasta')


def detect_ltr_for_sequence(record, output_dir):
    """
    Detect LTR in a single sequence using BLASTn.
    """
    record_len = record_len = len(record.seq)
    record_name = record.id.split('#')[0].replace('/', '_')
    sequence_folder = os.path.join(output_dir, record_name)
    os.makedirs(sequence_folder, exist_ok=True)
    sequence_file = os.path.join(sequence_folder, record_name)

    # Write the current sequence to a temporary file
    SeqIO.write(record, sequence_file, 'fasta')

    # Prepare BLAST database and output paths
    database_file = os.path.join(sequence_folder, 'temp_blast_database')
    makeblastdb_cmd = (
        f'makeblastdb -in {sequence_file} -dbtype nucl -out {database_file}'
    )

    try:
        subprocess.run(
            makeblastdb_cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        logging.error(
            f'\nmakeblastdb encountered an error for sequence {record_name} and returned error code {e.returncode}\n{e.stdout}\n{e.stderr}\n')
        return None

    blast_cmd = (
        f'blastn -query {sequence_file} -db {database_file} '
        f'-outfmt "6 qseqid qstart qend sstart send" '
        f'-evalue 0.05 -word_size 11 -gapopen 5 -gapextend 2 -reward 2 -penalty -3'
    )  # Set a higher evalue for self-blast

    try:
        # Execute the command
        result = subprocess.run(
            blast_cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        logging.error(
            f'\nblast encountered an error for sequence {record_name} and returned error code {e.returncode}\n{e.stdout}\n{e.stderr}\n')
        return None

    blast_out = result.stdout.strip()
    if not blast_out:
        return None

    # Process BLAST output
    # Process BLAST output
    columns = ['qseqid', 'qstart', 'qend', 'sstart', 'send']
    df = pd.DataFrame([x.split('\t') for x in blast_out.split('\n')], columns=columns)

    if df.empty:
        return None

    # Convert columns to integer
    df[['qstart', 'qend', 'sstart', 'send']] = df[
        ['qstart', 'qend', 'sstart', 'send']
    ].astype(int)

    df_LTR = df[
        (df['qend'] - df['qstart'] >= 150)
        & (df['qstart'] != df['sstart'])
        & (df['qend'] != df['send'])
        & (df['sstart'] < df['send'])
        & (df['sstart'] > df['qend'])
        & (df['sstart'] > df['qstart'])
    ].copy()

    if not df_LTR.empty:
        df_LTR['diff'] = df['send'] - df['qstart']
        df_LTR.reset_index(drop=True, inplace=True)

        # Find the row with the largest difference
        LTR_largest = df_LTR.iloc[df_LTR['diff'].idxmax()]

        # Check if the terminal repeat spans the majority of the query sequence. Because the query was extended
        if abs(LTR_largest['send'] - LTR_largest['qstart']) >= (record_len - 200):
            # Because BLAST uses index starting from 1, modify the start position
            LTR_boundary = [
                LTR_largest['qstart'] - 1,
                LTR_largest['qend'],
                LTR_largest['sstart'] - 1,
                LTR_largest['send'],
            ]
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
