# Standard library imports
import logging
import gzip
import os
import shutil
import subprocess
import time
import traceback

import click
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .boundarycrop import find_boundary_and_crop
# Local imports
from .functions import (blast, cd_hit_est, check_bed_uniqueness,
                        copy_files_with_start_pattern, decompress_gzip,
                        extract_fasta, fasta_file_to_dict,
                        handle_sequence_low_copy, handle_sequence_skipped,
                        modify_fasta_headers, multi_seq_dotplot,
                        parse_cd_hit_est_result, prcyan, prgre, process_lines,
                        remove_files_with_start_pattern, rename_cons_file,
                        rename_files_based_on_dict, repeatmasker,
                        repeatmasker_output_classify,
                        update_low_copy_cons_file)
from .MSAcluster import clean_and_cluster_MSA
from .orfdomain import PlotPfam, prepare_pfam_database
from .seqclass import SeqObject
from .TEaid import check_self_alignment


# Define a function to check the progress file, which will be used to continue analysis if program exited prematurely
def check_progress_file(progress_file_path):
    # Read the progress file into a pandas DataFrame
    df = pd.read_csv(progress_file_path)

    # Calculate skipped and low copy element number
    skipped_count = df[df['status'].str.strip().str.lower() == 'skipped'].shape[0]
    low_copy_count = df[
        df['low_copy'].astype(str).str.strip().str.lower() == 'true'
    ].shape[0]
    unknown_n = df[df['reclassified_type'].str.contains('Unknown', na=False)].shape[0]
    classifid_n = df[df['reclassified_type'].str.contains('/', na=False)].shape[0]

    # Calculate classified proportion
    if classifid_n != 0 and unknown_n != 0:
        classified_pro = classifid_n / (unknown_n + classifid_n)
    else:
        classified_pro = (
            0  # Set a default value (or any other value deemed appropriate)
        )

    # Get unique 'input_name' values
    local_completed_sequences = df['input_name'].unique().tolist()

    return local_completed_sequences, skipped_count, low_copy_count, classified_pro


def change_permissions_recursive(input_dir, mode):
    try:
        for dirpath, dirnames, filenames in os.walk(input_dir):
            os.chmod(dirpath, mode)
            for filename in filenames:
                os.chmod(os.path.join(dirpath, filename), mode)
    except PermissionError:
        click.echo(
            "TEtrimmer don't have right to change permissions. Pleas use sudo to run TEtrimmer"
        )
        return False
    return True


def calculate_genome_length(genome_file, outfile=None):
    """
    Calculate the length of each sequence in a genome file in FASTA format
    and write the lengths to an output file.

    :param genome_file: str, path to genome file in FASTA format
    :return: str, path to the output file containing sequence names and lengths
    """
    genome_lengths = {}
    with open(genome_file, 'r') as f:
        current_seq = None
        current_length = 0
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq is not None:
                    genome_lengths[current_seq] = current_length
                # If chromosome header contains empty spaces, only consider the content before the first space
                current_seq = line[1:].split(' ')[0]
                current_length = 0
            else:
                current_length += len(line)
        # Add length of last sequence to dictionary
        if current_seq is not None:
            genome_lengths[current_seq] = current_length

    # Write lengths to output file
    if not outfile:
        output_file = genome_file + '.length'
    else:
        output_file = outfile

    with open(output_file, 'w') as out:
        for seq_name, length in genome_lengths.items():
            out.write(f'{seq_name}\t{length}\n')

    return output_file


def check_database(genome_file, idx_dir=None, search_type='blast'):
    """
    Checks if the BLAST or MMseqs2 database and genome length file exist.
    If they do not exist, create them.

    :param genome_file: str, path to genome file (containing genome name)
    :param search_type: str, type of search tool to use ('blast' or 'mmseqs')
    """

    # Check if the genome file exists
    if not os.path.isfile(genome_file):
        raise FileNotFoundError(f'Genome file {genome_file} not found.')

    # Check if the genome file is gzipped
    is_gzipped = genome_file.endswith('.gz')
    if is_gzipped:
        decompressed_genome_file = decompress_gzip(genome_file)
    else:
        decompressed_genome_file = genome_file

    mmseqs_database_dir = None
    symlinked_genome_file = genome_file

    # If idx_dir is provided, set the database output directory to idx_dir/genome_index/genome_name
    if idx_dir:
        # Create the database directory if it does not exist
        database_dir = os.path.join(idx_dir, 'genome_index')
        # Create the database directory if it does not exist
        os.makedirs(database_dir, exist_ok=True)
        # Get basename and remove file extension if present
        database_name = os.path.splitext(os.path.basename(decompressed_genome_file))[0]
        # Create a symlink to the genome file in the database directory
        symlinked_genome_file = os.path.join(database_dir,os.path.basename(genome_file))
        if not os.path.exists(symlinked_genome_file):
            os.symlink(genome_file, symlinked_genome_file)
    else:
        # Set the database output directory to the same directory as the genome file
        database_dir = os.path.dirname(decompressed_genome_file)
        database_name = os.path.splitext(os.path.basename(decompressed_genome_file))[0]

    # Compose the output database file path
    output_database_path = os.path.join(database_dir, database_name)

    try:
        if search_type == 'blast':
            blast_database_file = output_database_path + '.nin'
            if not os.path.isfile(blast_database_file):
                logging.warning(f"Blast database doesn't exist at: {output_database_path}")

                try:
                    makeblastdb_cmd = (
                        f'makeblastdb -in {decompressed_genome_file} -dbtype nucl -out {output_database_path} '
                    )

                    logging.info(f'Creating BLAST database: {makeblastdb_cmd}')

                    subprocess.run(
                        makeblastdb_cmd,
                        shell=True,
                        check=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                    )

                except subprocess.CalledProcessError as e:
                    logging.error(f'makeblastdb failed with exit code {e.returncode}\n{e.stdout}\n{e.stderr}')
                    exit(1)

        elif search_type == 'mmseqs':
            mmseqs_database_dir = output_database_path + '_db'
            if not os.path.isdir(mmseqs_database_dir):
                print('\nMMseqs2 database does not exist. Creating MMseqs2 database...\n')

                try:
                    mmseqs_createdb_cmd = (
                        f'mmseqs createdb {decompressed_genome_file} {mmseqs_database_dir}'
                    )
                    subprocess.run(
                        mmseqs_createdb_cmd,
                        shell=True,
                        check=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                    )

                except subprocess.CalledProcessError as e:
                    logging.error(f'mmseqs createdb failed with exit code {e.returncode}\n{e.stdout}\n{e.stderr}')
                    exit(1)

                # Check if MMseqs2 database files have been created
                if os.path.exists(f'{mmseqs_database_dir}.dbtype'):
                    print(f'MMseqs2 database created successfully: {mmseqs_database_dir}')

                    # Create index for the MMseqs2 database
                    print('Creating index for MMseqs2 database...')
                    mmseqs_createindex_cmd = (
                        f'mmseqs createindex {mmseqs_database_dir} {mmseqs_database_dir}_tmp '
                        f'--search-type 3'
                    )
                    result = subprocess.run(
                        mmseqs_createindex_cmd,
                        shell=True,
                        check=True,
                        stderr=subprocess.PIPE,
                        text=True,
                    )
                    error_output = result.stderr
                    if error_output:
                        logging.error(f'Error creating MMseqs2 index: {error_output}\n')
                        exit(1)
                    else:
                        logging.info('MMseqs2 index created successfully.')
                else:
                    logging.error(
                        f'Error: MMseqs2 database files not found for {decompressed_genome_file}.')
                    exit(1)

        # Check if genome length files exist, otherwise create them in the same folder with genome file
        length_file = output_database_path + '.length'
        if not os.path.isfile(length_file):
            logging.info(f'File with genome lengths not found. Writing: {length_file}')
            calculate_genome_length(decompressed_genome_file, output=length_file)

        # Check if .fai index file exists, otherwise create it using samtools faidx
        # Index gzip file directly if available
        fai_file = os.path.join(database_dir,os.path.basename(genome_file)) + '.fai'

        if not os.path.isfile(fai_file):
            logging.info(
                f'Index file {fai_file} not found. Creating it using samtools faidx ...\n'
            )
            faidx_cmd = f'samtools faidx {genome_file} -o {fai_file}'

            try:
                subprocess.run(
                    faidx_cmd,
                    shell=True,
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )

            except FileNotFoundError:
                logging.error(
                    "'samtools' command not found. Please ensure 'samtools' is correctly installed."
                )
                exit(1)

            except subprocess.CalledProcessError as e:
                logging.error(f'\nsamtools faidx failed with error code {e.returncode}\n{e.stdout}\n{e.stderr}\n')
                exit(1)
    finally:
        # Remove the decompressed file if it was created
        if is_gzipped and os.path.isfile(decompressed_genome_file):
            os.remove(decompressed_genome_file)

    return (mmseqs_database_dir, output_database_path, length_file, fai_file, symlinked_genome_file)

# Print iterations progress
def printProgressBar(
    iteration,
    total,
    prefix='',
    suffix='',
    decimals=1,
    length=100,
    fill='█',
    printEnd='\r',
    final=False,
    ):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ('{0:.' + str(decimals) + 'f}').format(100 * (iteration / float(total)))
    if not final:
        percent = min(float(percent), 99)
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    # click.echo(f'\r{prefix} |{bar}| {iteration}/{total} = {percent}% {suffix}', nl=False)
    click.echo(f'{prefix} |{bar}| {iteration}/{total} = {percent}% {suffix}', nl=True)
    # print(f'{prefix} |{bar}| {iteration}/{total} = {percent}% {suffix}', end='\n', flush=True)

    # Print New Line on Complete
    if iteration == total:
        click.echo()


def separate_sequences(input_file, output_dir, continue_analysis=False):
    """
    Separates input file into single separate FASTA files and creates objects for each input sequence
    """
    os.makedirs(output_dir, exist_ok=True)
    seq_list = []

    def open_file(file_path):
        if file_path.endswith('.gz'):
            return gzip.open(file_path, 'rt')
        else:
            return open(file_path, 'r')

    if not continue_analysis:
        print(
            "TE Trimmer is modifying sequence names; any occurrence of '/', '-', ':', '...', '|' and empty spaces before '#' "
            "will be converted to '_'.\n"
            "You can find the original and modified names in the 'Sequence_name_mapping.txt' file in the output directory.\n"
        )
        # Initialize the name mapping file
        name_mapping_file = os.path.join(
            os.path.dirname(output_dir), 'Sequence_name_mapping.txt'
        )

        detected_pound = False
        with open_file(input_file) as fasta_file, open(
            name_mapping_file, 'w'
        ) as mapping_file:
            # Write header to the mapping file
            mapping_file.write('original_input_seq_name\tTEtrimmer_modified_seq_name\n')
            id_list = []
            # Required to add suffix 'fasta', as this pattern will be used for file deletion later
            for record in SeqIO.parse(fasta_file, 'fasta'):
                # Check if '#' is in seq.id. If '#' is found, the string before '#' becomes the seq_name and the string
                # after '#' is the seq_TE_type
                if len(record.id.split('#')) > 1:
                    detected_pound = True
                    sanitized_id = (
                        record.id.split('#')[0]
                        .replace('/', '_')
                        .replace(' ', '_')
                        .replace('-', '_')
                        .replace(':', '_')
                        .replace('...', '_')
                        .replace('|', '_')
                    )
                    te_type = record.id.split('#')[1]

                    # Normally SeqIO.parse only takes content before " " as record.id. Separate with " " to make
                    # the code more reliable
                    te_type = te_type.split(' ')[0]

                else:
                    sanitized_id = (
                        record.id.replace('/', '_')
                        .replace(' ', '_')
                        .replace('-', '_')
                        .replace(':', '_')
                        .replace('...', '_')
                    )
                    te_type = 'Unknown'
                    # modify header to add #Unknown
                    record.id = f'{record.id}#{te_type}'
                    record.description = record.id

                # double check if sanitized_id is unique. If not, modify sanitized_id
                if sanitized_id not in id_list:
                    id_list.append(sanitized_id)
                else:
                    # print(f"Duplicated seq_name {sanitized_id} during separate_sequences.")
                    id_list.append(sanitized_id)
                    count = id_list.count(sanitized_id)
                    sanitized_id = f'{sanitized_id}_n{count}'

                # Write original and modified names to the mapping file
                mapping_file.write(f'{record.id}\t{sanitized_id}\n')

                # Define output file name
                output_filename = os.path.join(output_dir, f'{sanitized_id}.fasta')
                seq_obj = SeqObject(
                    str(sanitized_id), str(output_filename), len(record.seq), te_type
                )

                # Store all input file information (object) to seq_list
                seq_list.append(seq_obj)

                # Convert sequence name to sanitized_id
                record.id = sanitized_id

                # Write single FASTA file using sanitized name
                # the record.id is now same as sanitized_id, which only contains string before #
                # the record.description = f"{record.id}#{te_type}"
                with open(output_filename, 'w') as output_file:
                    SeqIO.write(record, output_file, 'fasta')

            if detected_pound:
                print(
                    "TEtrimmer detected instances of '#' in your input FASTA sequence headers. The string before "
                    "'#' is denoted as the seq_name, and the string after '#' is denoted as the TE type.\n"
                )
        print('Finish to generate single sequence files.\n')

    elif continue_analysis:
        # If continue_analysis is 'True', generate seq_list based on single FASTA files
        for filename in os.listdir(output_dir):
            file = os.path.join(output_dir, filename)
            with open(file, 'r') as fasta_file:
                for record in SeqIO.parse(fasta_file, 'fasta'):
                    # Get sanitized_id from single FASTA file name
                    sanitized_id = os.path.splitext(filename)[0]

                    # single FASTA file name is the same as sanitized_id and record.id
                    te_type = record.description.split('#')[-1]
                    te_type = te_type.split(' ')[0]
                    seq_obj = SeqObject(
                        str(sanitized_id), str(file), len(record.seq), te_type
                    )
                    seq_list.append(seq_obj)
        print(
            '\nFinished to read single sequence files generated by previous analysis.\n'
        )

    single_fasta_n = len(seq_list)

    return seq_list, single_fasta_n


def repeatmasker_classification(
    final_unknown_con_file,
    final_classified_con_file,
    classification_dir,
    num_threads,
    progress_file,
    final_con_file,
    proof_curation_dir,
    perfect_proof,
    good_proof,
    intermediate_proof,
    need_check_proof,
    low_copy_dir,
    hmm,
    hmm_dir,
):
    if os.path.exists(final_unknown_con_file) and os.path.exists(
        final_classified_con_file
    ):
        temp_repeatmasker_dir = os.path.join(
            classification_dir, 'temp_repeatmasker_classification'
        )
        reclassified_recording_path = os.path.join(
            temp_repeatmasker_dir, 'Reclassified_recoring.txt'
        )
        os.makedirs(temp_repeatmasker_dir, exist_ok=True)
        classification_out = repeatmasker(
            final_unknown_con_file,
            final_classified_con_file,
            temp_repeatmasker_dir,
            thread=num_threads,
            classify=True,
        )

        if classification_out:
            repeatmasker_out = os.path.join(
                temp_repeatmasker_dir, 'temp_TEtrimmer_unknown_consensus.fasta.out'
            )
            reclassified_dict = repeatmasker_output_classify(
                repeatmasker_out, progress_file, min_iden=70, min_len=80, min_cov=0.5
            )
            if reclassified_dict:
                click.echo(
                    f'\n{len(reclassified_dict)} TE elements were re-classified by the '
                    f'final classification module.'
                )

                # Update final consensus file
                rename_cons_file(final_con_file, reclassified_dict)
                rename_files_based_on_dict(proof_curation_dir, reclassified_dict)
                rename_files_based_on_dict(perfect_proof, reclassified_dict)
                rename_files_based_on_dict(good_proof, reclassified_dict)
                rename_files_based_on_dict(intermediate_proof, reclassified_dict)
                rename_files_based_on_dict(need_check_proof, reclassified_dict)
                rename_files_based_on_dict(
                    low_copy_dir, reclassified_dict, seq_name=True
                )
                if hmm:
                    rename_files_based_on_dict(hmm_dir, reclassified_dict)

                # Write reclassified ID into a file
                # Open the file for writing
                with open(reclassified_recording_path, 'w') as file:
                    # Iterate through the dictionary and write each key-value pair to the file
                    for key, value in reclassified_dict.items():
                        file.write(f'{key}\t{value}\n')

            else:
                click.echo(
                    '0 TE elements were re-classified by the final classification module.'
                )

    else:
        prcyan('\nThe final classification module failed.')
        prgre(
            '\nThis does not affect the final TE consensus sequences You can choose to ignore this error.\n'
        )


def merge_cons(
    classification_dir,
    final_con_file,
    progress_file,
    cd_hit_est_final_merged,
    num_threads,
):
    # Do first round of CD-HIT-EST
    cd_hit_merge_output_round1 = os.path.join(
        classification_dir, 'TEtrimmer_consensus_merged_round1.fasta'
    )
    cd_hit_merge_output_round1_clstr = f'{cd_hit_merge_output_round1}.clstr'

    # Round 1 merge only requires that the alignment coverage for the shorter sequence is greater than 0.9
    # and the similarity is greater than 0.9
    cd_hit_est(
        final_con_file,
        cd_hit_merge_output_round1,
        identity_thr=0.9,
        aL=0,
        aS=0.9,
        s=0,
        thread=num_threads,
    )

    # Read progress file
    progress_df = pd.read_csv(progress_file)

    # Create a dictionary with sequence names as keys
    sequence_info = {}
    for index, row in progress_df.iterrows():
        sequence_name = row['consensus_name']
        evaluation = (
            row['evaluation'] if pd.notna(row['evaluation']) else 'Unknown'
        )  # Default value for NaN
        te_type = (
            row['reclassified_type']
            if pd.notna(row['reclassified_type'])
            else 'Unknown'
        )
        length = (
            row['cons_length'] if pd.notna(row['cons_length']) else 0
        )  # Default value for NaN
        sequence_info[sequence_name] = {
            'evaluation': evaluation,
            'type': te_type,
            'length': length,
        }

    # Parse cd-hit-est result, clusters is a dictionary, the key is cluster number, value is a list
    # contain all sequence names in this cluster
    clusters, detailed_clusters = parse_cd_hit_est_result(
        cd_hit_merge_output_round1_clstr
    )

    # Check if sequences scored "Perfect" and "Good" are included in the cluster and choose the longest sequence
    best_sequences = []  # Define list to store "Perfect" or "Good" sequence names

    # Define list to store sequence in clusters that do not contain "Perfect" or "Good" sequences
    sequence_for_round2 = []
    for cluster_name, sequences in clusters.items():
        perfect_sequences = []
        good_sequences = []
        best_seq = None

        for seq in sequences:
            if seq in sequence_info:
                evaluation = sequence_info[seq]['evaluation']
                length = sequence_info[seq]['length']
                if evaluation == 'Perfect':
                    perfect_sequences.append((seq, length))
                elif evaluation == 'Good':
                    good_sequences.append((seq, length))
        # Choose the longest "Perfect" sequence, if have Perfect
        if perfect_sequences:
            best_seq = max(perfect_sequences, key=lambda x: x[1])[0]
        # If no "Perfect", choose the longest "Good" sequence
        elif good_sequences:
            best_seq = max(good_sequences, key=lambda x: x[1])[0]

        if best_seq:
            best_sequences.append(best_seq)
        else:
            sequence_for_round2.extend(
                clusters[cluster_name]
            )  # extend() will create a flat list

    # Read the original consensus file
    consensus_sequences = SeqIO.parse(final_con_file, 'fasta')

    # Define temporary file to store "Perfect" and "Good" sequences
    temp_consensus_round1 = os.path.join(
        classification_dir, 'temp_consensus_round1.fasta'
    )

    # Define temporary file to store remaining sequences for second round of CD-HIT-EST
    temp_consensus_round2_input = os.path.join(
        classification_dir, 'temp_consensus_round2_input.fasta'
    )

    # Write sequences to files
    with open(temp_consensus_round1, 'w') as high_quality_file, open(
        temp_consensus_round2_input, 'w'
    ) as round2_file:
        for seq_record in consensus_sequences:
            # Sequence names in best_sequences and sequence_for_round2 do not contain classification
            seq_id = seq_record.id.split('#')[0]

            if seq_id in best_sequences:
                SeqIO.write(seq_record, high_quality_file, 'fasta')
            elif seq_id in sequence_for_round2:
                SeqIO.write(seq_record, round2_file, 'fasta')

    # Do second round of CD-HIT-EST based on temp_consensus_round2_input
    cd_hit_merge_output_round2 = os.path.join(
        classification_dir, 'TEtrimmer_consensus_merged_round2.fasta'
    )

    # Round 2 merge requires that the alignment coverage for the long and short sequence are both greater than 0.8
    # and the similarity is greater than 0.85
    cd_hit_est(
        temp_consensus_round2_input,
        cd_hit_merge_output_round2,
        identity_thr=0.85,
        aL=0.8,
        aS=0.8,
        s=0.8,
        thread=num_threads,
    )

    # Combine the two files into a merged file
    with open(temp_consensus_round1, 'r') as file1, open(
        cd_hit_merge_output_round2, 'r'
    ) as file2, open(cd_hit_est_final_merged, 'w') as combined_file:
        # Write contents of the first file
        for line in file1:
            combined_file.write(line)

        # Write contents of the second file
        for line in file2:
            combined_file.write(line)

    # Find sequence names that are not included inside in cd_hit_est_final_merged file
    # Parse the sequences in the original and merged files
    original_sequences = SeqIO.parse(final_con_file, 'fasta')
    merged_sequences = SeqIO.parse(cd_hit_est_final_merged, 'fasta')

    # Extract sequence IDs from both files and store to set
    original_ids = {seq_record.id.split('#')[0] for seq_record in original_sequences}
    merged_ids = {seq_record.id.split('#')[0] for seq_record in merged_sequences}

    # Find the difference between the two sets to identify sequence names not included in the merged file
    missing_ids = original_ids - merged_ids

    """
    # Based on missing_ids delete files in proof curation folder and HMM folder
    for missing_id in missing_ids:

        # if not, set evaluation_level to "Need_check". The "get" method will return the default value
        # when the key does not exist.
        evaluation_level = sequence_info.get(missing_id, {"evaluation": "Need_check"})["evaluation"]

        # Add '#' to the end of missing_id, this can avoid to delete 140 when the id is 14
        missing_id = f"{missing_id}#"
        if evaluation_level == "Perfect":
            remove_files_with_start_pattern(perfect_proof, missing_id, if_seq_name=False)
        elif evaluation_level == "Good":
            remove_files_with_start_pattern(good_proof, missing_id, if_seq_name=False)
        elif evaluation_level == "Reco_check":
            remove_files_with_start_pattern(intermediate_proof, missing_id, if_seq_name=False)
        elif evaluation_level == "Need_check":
            remove_files_with_start_pattern(need_check_proof, missing_id, if_seq_name=False)
        else:
            remove_files_with_start_pattern(low_copy_dir, missing_id, if_seq_name=False)

        if hmm:
            remove_files_with_start_pattern(hmm_dir, missing_ids)
    """
    click.echo('\nFinished to remove sequence duplications.\n')
    return sequence_info


def cluster_proof_anno_file(
    multi_dotplot_dir,
    final_con_file_no_low_copy,
    continue_analysis,
    cluster_proof_anno_dir,
    num_threads,
    sequence_info,
    perfect_proof,
    good_proof,
    intermediate_proof,
    need_check_proof,
):
    # Load fast file to a dictionary, key is record.id, value is record project
    # When separate_name is true, the key of the dictionary will be the sequence name separated by '#'
    final_con_file_no_low_copy_dict = fasta_file_to_dict(
        final_con_file_no_low_copy, separate_name=True
    )

    # Clean cluster_proof_anno_dir when --continue_analysis is on.
    if continue_analysis:
        # When the start pattern isn't given, all files inside the folder will be removed
        remove_files_with_start_pattern(cluster_proof_anno_dir)
        remove_files_with_start_pattern(multi_dotplot_dir)

    # Do CD-HIT-EST for final consensus file without low copy elements
    final_con_file_no_low_copy_cd_out = f'{final_con_file_no_low_copy}_cd.fa'
    final_con_file_no_low_copy_clstr = f'{final_con_file_no_low_copy_cd_out}.clstr'

    # Round 1 merge only requires that the alignment coverage for the shorter sequence is greater than 0.9
    # and the similarity is greater than 0.9
    cd_hit_est(
        final_con_file_no_low_copy,
        final_con_file_no_low_copy_cd_out,
        identity_thr=0.9,
        aL=0,
        aS=0.9,
        s=0,
        thread=num_threads,
    )
    clusters_proof_anno, detailed_clusters_proof_anno = parse_cd_hit_est_result(
        final_con_file_no_low_copy_clstr
    )

    for (
        cluster_name_proof_anno,
        seq_info_proof_anno,
    ) in detailed_clusters_proof_anno.items():
        # Create cluster folder
        cluster_folder = os.path.join(cluster_proof_anno_dir, cluster_name_proof_anno)
        os.makedirs(cluster_folder, exist_ok=True)

        seq_info_proof_anno_len = len(seq_info_proof_anno)

        cluster_record_list = []
        for i in range(seq_info_proof_anno_len):
            try:
                seq_length_proof_anno = seq_info_proof_anno[i][0]
            except Exception:
                seq_length_proof_anno = None

            seq_name_proof_anno = seq_info_proof_anno[i][1]

            try:
                seq_per_proof_anno = seq_info_proof_anno[i][2]
            except Exception:
                seq_per_proof_anno = None
            try:
                seq_direction_proof_anno = seq_info_proof_anno[i][3]
            except Exception:
                seq_direction_proof_anno = None

            # Copy sequence files into cluster folder
            # if not, set evaluation_level to "Need_check". The "get" method will return the default value
            # when the key does not exist.
            evaluation_level = sequence_info.get(
                seq_name_proof_anno, {'evaluation': 'Need_check'}
            )['evaluation']

            # Add '#' to the end of seq_name_proof_anno, this can avoid to delete 140 when the id is 14
            seq_name_proof_anno_m = f'{seq_name_proof_anno}#'
            if evaluation_level == 'Perfect':
                copy_files_with_start_pattern(
                    perfect_proof,
                    seq_name_proof_anno_m,
                    cluster_folder,
                    seq_length_proof_anno,
                    seq_per_proof_anno,
                    evaluation_level,
                )
            elif evaluation_level == 'Good':
                copy_files_with_start_pattern(
                    good_proof,
                    seq_name_proof_anno_m,
                    cluster_folder,
                    seq_length_proof_anno,
                    seq_per_proof_anno,
                    evaluation_level,
                )
            elif evaluation_level == 'Reco_check':
                copy_files_with_start_pattern(
                    intermediate_proof,
                    seq_name_proof_anno_m,
                    cluster_folder,
                    seq_length_proof_anno,
                    seq_per_proof_anno,
                    evaluation_level,
                )
            elif evaluation_level == 'Need_check':
                copy_files_with_start_pattern(
                    need_check_proof,
                    seq_name_proof_anno_m,
                    cluster_folder,
                    seq_length_proof_anno,
                    seq_per_proof_anno,
                    evaluation_level,
                )

            # Plot multiple sequence dotplot when more than one sequence are included inside one cluster
            if seq_info_proof_anno_len > 1:
                # Extract record from dictionary
                cluster_record = final_con_file_no_low_copy_dict.get(
                    seq_name_proof_anno
                )

                # When the sequence direction is negative, reverse complement it
                if cluster_record is not None and seq_direction_proof_anno == '-':
                    rev_comp_cluster_record_seq = (
                        cluster_record.seq.reverse_complement()
                    )
                    cluster_record = SeqRecord(
                        rev_comp_cluster_record_seq,
                        id=cluster_record.id,
                        description='',
                    )
                if cluster_record is not None:
                    cluster_record_list.append(cluster_record)

        if len(cluster_record_list) > 1:
            # Define and write cluster fasta file a
            cluster_fasta = os.path.join(
                multi_dotplot_dir, f'{cluster_name_proof_anno}.fa'
            )
            SeqIO.write(cluster_record_list, cluster_fasta, 'fasta')

            # Do multiple sequence dotplot
            multi_dotplot_pdf = multi_seq_dotplot(
                cluster_fasta, multi_dotplot_dir, cluster_name_proof_anno
            )

            # Move muti_dotplot_pdf to proof curation cluster folder
            if os.path.isfile(multi_dotplot_pdf):
                shutil.copy(multi_dotplot_pdf, cluster_folder)


#####################################################################################################
# Code block: Define analyze_sequence function
#####################################################################################################


def analyze_sequence_helper(params):
    return analyze_sequence(*params)


def analyze_sequence(
    seq_obj,
    genome_file,
    MSA_dir,
    min_blast_len,
    min_seq_num,
    max_msa_lines,
    top_msa_lines,
    max_cluster_num,
    cons_thr,
    ext_thr,
    ex_step,
    classification_dir,
    max_extension,
    gap_thr,
    gap_nul_thr,
    crop_end_thr,
    crop_end_win,
    crop_end_gap_thr,
    crop_end_gap_win,
    start_patterns,
    end_patterns,
    output_dir,
    pfam_dir,
    mini_orf,
    single_fasta_n,
    hmm,
    hmm_dir,
    check_extension_win,
    debug,
    progress_file,
    classify_unknown,
    classify_all,
    final_con_file,
    final_con_file_no_low_copy,
    final_unknown_con_file,
    final_classified_con_file,
    low_copy_dir,
    fast_mode,
    error_files,
    plot_skip,
    skipped_dir,
    plot_query,
    engine,
    proof_curation_dir,
):
    #####################################################################################################
    # Code block: Set different elongation number for different elements and do BLAST search
    #####################################################################################################

    try:
        # Get query fasta file path
        seq_name = seq_obj.get_seq_name()
        seq_type = seq_obj.get_old_TE_type()
        seq_file = seq_obj.get_input_fasta()  # Return complete file path

        # Since DNA element are significantly shorter than LTR and LINE elements, adjust default parameters
        if 'DNA' in seq_type:
            ex_step = 500
            max_extension = 7000
            min_blast_len = 150
            crop_end_gap_win = 100
            check_extension_win = 50

        # The average length of SINE elements is around 500 bp, adjust default parameters
        if 'SINE' in seq_type:
            ex_step = 200
            max_extension = 1400
            min_blast_len = 80
            crop_end_gap_win = 50
            check_extension_win = 50

        if 'Helitron' in seq_type:
            ex_step = 500
            max_extension = 7000
            min_blast_len = 150
            crop_end_gap_win = 100
            check_extension_win = 50

        if 'MITE' in seq_type:
            ex_step = 100
            max_extension = 500
            min_blast_len = 50
            crop_end_gap_win = 40
            check_extension_win = 50

        # run BLAST search for each FASTA file and return a BED file absolute path
        bed_out_file_dup, blast_hits_count, blast_out_file = blast(
            seq_file,
            genome_file,
            MSA_dir,
            min_length=min_blast_len,
            task='blastn',
            seq_obj=seq_obj,
            search_type=engine,
        )

    except Exception as e:
        with open(error_files, 'a') as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f'Error while running blast for sequence: {seq_name}\n')
            f.write(tb_content + '\n\n')
        prcyan(
            f'Error while running blast for sequence: {seq_name}. Main Error: {str(e)}. \n'
            f'Trace back content: {tb_content}\n'
        )
        return

    #####################################################################################################
    # Code block: Perform ORF and PFAM prediction for input sequences
    #####################################################################################################
    input_orf_domain_plot = None
    try:
        input_orf_pfam_obj = PlotPfam(
            seq_file,
            MSA_dir,
            pfam_database_dir=pfam_dir,
            mini_orf=mini_orf,
            after_tetrimmer=False,
        )

        # "run_getorf()" function will return 'True' if any ORF was detected. Otherwise, it will return 'False'.
        if input_orf_pfam_obj.run_getorf():
            # "run_pfam_scan()" will return 'True' if any PFAM domains were found. Otherwise, it will return 'False'.
            pfam_scan_result = input_orf_pfam_obj.run_pfam_scan()
            input_orf_domain_plot = input_orf_pfam_obj.orf_domain_plot()

    except Exception:
        input_orf_domain_plot = None
        with open(error_files, 'a') as f:
            tb_content = traceback.format_exc()
            f.write(
                f'Error when doing ORF and PFAM prediction for input sequence {seq_name}\n'
            )
            f.write(tb_content + '\n\n')
        # prcyan(f"Error while performing input sequence ORF and PFAM predictions: {seq_name}. Main Error: {str(e)}. \n"
        # f"Trace back content: {tb_content}\n")

    #####################################################################################################
    # Code block: Check BLAST hit number
    #####################################################################################################

    try:
        # Check if BLAST hit number is exactly 0. If so, skip this sequence
        if blast_hits_count == 0:
            click.echo(f'\n{seq_name} is skipped due to blast hit number is 0\n')
            handle_sequence_skipped(
                seq_obj,
                progress_file,
                debug,
                MSA_dir,
                classification_dir,
                skip_proof_dir=skipped_dir,
            )
            return

        # Check if BLAST hit number is smaller than "min_seq_num"; do not include "min_seq_num"
        elif blast_hits_count != 0 and blast_hits_count < min_seq_num:
            check_low_copy, blast_full_length_n, found_match, TE_aid_plot = (
                check_self_alignment(
                    seq_obj,
                    seq_file,
                    MSA_dir,
                    genome_file,
                    blast_hits_count,
                    blast_out_file,
                    plot_skip=plot_skip,
                )
            )

            if check_low_copy is True:
                # Update terminal repeat and BLAST full length number, remove low-copy intermediate files
                handle_sequence_low_copy(
                    seq_obj,
                    progress_file,
                    debug,
                    MSA_dir,
                    classification_dir,
                    found_match=found_match,
                    blast_full_length_n=blast_full_length_n,
                    te_aid_plot=TE_aid_plot,
                    orf_plot=input_orf_domain_plot,
                    low_copy_dir=low_copy_dir,
                )

                # Integrate low-copy element sequences into consensus file
                update_low_copy_cons_file(
                    seq_obj,
                    final_con_file,
                    final_unknown_con_file,
                    final_classified_con_file,
                    low_copy_dir,
                    TE_aid_plot,
                )
            else:
                click.echo(
                    f'\n{seq_name} was skipped because the BLAST hit number is smaller than {min_seq_num} '
                    f'and check_low_copy is {check_low_copy}.\n'
                )

                # handle_sequence_skipped will update skipped status
                handle_sequence_skipped(
                    seq_obj,
                    progress_file,
                    debug,
                    MSA_dir,
                    classification_dir,
                    plot_skip=plot_skip,
                    te_aid_plot=TE_aid_plot,
                    skip_proof_dir=skipped_dir,
                    orf_plot=input_orf_domain_plot,
                )

            return  # if BLAST hit number is smaller than 10, code will execute next FASTA file

    except Exception as e:
        # Add sequence to skip if check low-copy module returns errors
        handle_sequence_skipped(
            seq_obj,
            progress_file,
            debug,
            MSA_dir,
            classification_dir,
            skip_proof_dir=skipped_dir,
        )
        with open(error_files, 'a') as f:
            # Return the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(
                f'\nError while checking low-copy status for sequence: {seq_name}\n'
            )
            f.write(tb_content + '\n\n')
        prcyan(
            f'\nError while checking low-copy status for sequence: {seq_name}. Error: {str(e)}\n'
        )
        prgre(
            '\nThe low-copy TE check module is an optional additional analysis. You may ignore this error, as it '
            'will not affect the final result significantly.\n'
        )
        return

    #####################################################################################################
    # Code block: Separate MSA based on sequence relatedness
    #####################################################################################################

    try:
        # Remove duplicated lines
        bed_out_file = check_bed_uniqueness(MSA_dir, bed_out_file_dup)

        # Test bed_out_file line number and extract the longest lines for process_lines() function.
        # The threshold represents the maximum number of lines to keep for MSA.
        # top_longest_lines_count means the number of sequences with top length.
        # For example, if threshold = 100, top_longest_lines_count = 50, then 50 sequences will be
        # randomly selected from the remaining sequences.
        # top_msa_lines has to be equal to or smaller than max_msa_lines.
        bed_out_filter_file = process_lines(
            bed_out_file,
            MSA_dir,
            threshold=max_msa_lines,
            top_longest_lines_count=top_msa_lines,
        )

        # Extract FASTA from bed_out_filter_file
        # Return fasta_out_flank_file absolute path of FASTA file
        # It is imperative to group the sequences in the first round of MSA; extend for both ends is 0.
        fasta_out_flank_file, bed_out_flank_file = extract_fasta(
            bed_out_filter_file,
            genome_file,
            MSA_dir,
            left_ex=0,
            right_ex=0,
            nameonly=True,
        )

        # Return 'False' if cluster number is 0 (sequence number in each cluster is smaller than 10).
        # Otherwise, return the subset BED and alignment files.
        cluster_MSA_result = clean_and_cluster_MSA(
            fasta_out_flank_file,
            bed_out_filter_file,
            MSA_dir,
            clean_column_threshold=0.02,
            min_length_num=min_seq_num,
            cluster_num=max_cluster_num,
            cluster_col_thr=100,
            fast_mode=fast_mode,
        )
    except Exception as e:
        with open(error_files, 'a') as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f'Error while grouping MSA: {seq_name}\n')
            f.write(tb_content + '\n\n')
        prcyan(
            f'\nError while grouping MSA for sequence: {seq_name}. Error: {str(e)}\n'
        )
        prcyan('\n' + tb_content + '\n')
        return

    #####################################################################################################
    # Code block: Perform find_boundary_and_crop on clustered MSA if necessary
    #####################################################################################################

    try:
        # cluster_false means too few sequences were found in clusters from MSA (all_cluster_size < 10); TE Trimmer will skip this sequence.
        if cluster_MSA_result is False:
            check_low_copy, blast_full_length_n, found_match, TE_aid_plot = (
                check_self_alignment(
                    seq_obj,
                    seq_file,
                    MSA_dir,
                    genome_file,
                    blast_hits_count,
                    blast_out_file,
                    plot_skip=plot_skip,
                )
            )

            if check_low_copy is True:
                # Update terminal repeat and BLAST full-length number, remove low-copy intermediate files
                handle_sequence_low_copy(
                    seq_obj,
                    progress_file,
                    debug,
                    MSA_dir,
                    classification_dir,
                    found_match=found_match,
                    blast_full_length_n=blast_full_length_n,
                    te_aid_plot=TE_aid_plot,
                    orf_plot=input_orf_domain_plot,
                    low_copy_dir=low_copy_dir,
                )
                update_low_copy_cons_file(
                    seq_obj,
                    final_con_file,
                    final_unknown_con_file,
                    final_classified_con_file,
                    low_copy_dir,
                    TE_aid_plot,
                )

            else:
                click.echo(
                    f'\n{seq_name} was skipped because the sequence number in each cluster was smaller '
                    f'than {min_seq_num} and check_low_copy is {check_low_copy}.\n'
                )
                handle_sequence_skipped(
                    seq_obj,
                    progress_file,
                    debug,
                    MSA_dir,
                    classification_dir,
                    plot_skip=plot_skip,
                    te_aid_plot=TE_aid_plot,
                    skip_proof_dir=skipped_dir,
                    orf_plot=input_orf_domain_plot,
                )
            return
        else:
            cluster_bed_files_list, fasta_out_flank_mafft_gap_rm = cluster_MSA_result

            # modify fasta_out_flank_mafft_gap_rm fasta header based on the bed file, this can allow the
            # extension function in the final GUI.
            # For example: change 1(+) to scaffold_1:23256-24757(+)
            fasta_out_flank_mafft_gap_rm_nm = modify_fasta_headers(
                bed_out_filter_file, fasta_out_flank_mafft_gap_rm
            )

            # cluster_pattern_alignment_list has the same index as cluster_bed_files_list
            all_inner_skipped = True
            for i in range(len(cluster_bed_files_list)):
                try:
                    find_boundary_result = find_boundary_and_crop(
                        cluster_bed_files_list[i],
                        genome_file,
                        MSA_dir,
                        pfam_dir,
                        seq_obj,
                        hmm,
                        classify_all,
                        classify_unknown,
                        error_files,
                        plot_query,
                        classification_dir,
                        final_con_file,
                        final_con_file_no_low_copy,
                        proof_curation_dir,
                        hmm_dir,
                        cons_threshold=cons_thr,
                        ext_threshold=ext_thr,
                        ex_step_size=ex_step,
                        max_extension=max_extension,
                        gap_threshold=gap_thr,
                        gap_nul_thr=gap_nul_thr,
                        crop_end_thr=crop_end_thr,
                        crop_end_win=crop_end_win,
                        crop_end_gap_thr=crop_end_gap_thr,
                        crop_end_gap_win=crop_end_gap_win,
                        start_patterns=start_patterns,
                        end_patterns=end_patterns,
                        mini_orf=mini_orf,
                        define_boundary_win=check_extension_win,
                        fast_mode=fast_mode,
                        engine=engine,
                        input_orf_pfam=input_orf_domain_plot,
                        debug=debug,
                        cluster_msa=fasta_out_flank_mafft_gap_rm_nm,
                    )
                except Exception:
                    return

                if not find_boundary_result:
                    continue
                elif find_boundary_result:
                    all_inner_skipped = False

            # Check the flag after the loop. If all inner clusters were skipped, write the progress file.
            if all_inner_skipped:
                check_low_copy, blast_full_length_n, found_match, TE_aid_plot = (
                    check_self_alignment(
                        seq_obj,
                        seq_file,
                        MSA_dir,
                        genome_file,
                        blast_hits_count,
                        blast_out_file,
                        plot_skip=plot_skip,
                    )
                )

                if check_low_copy is True:
                    # Update terminal repeat and BLAST full-length number, remove low-copy intermediate files.
                    handle_sequence_low_copy(
                        seq_obj,
                        progress_file,
                        debug,
                        MSA_dir,
                        classification_dir,
                        found_match=found_match,
                        blast_full_length_n=blast_full_length_n,
                        te_aid_plot=TE_aid_plot,
                        orf_plot=input_orf_domain_plot,
                        low_copy_dir=low_copy_dir,
                    )
                    update_low_copy_cons_file(
                        seq_obj,
                        final_con_file,
                        final_unknown_con_file,
                        final_classified_con_file,
                        low_copy_dir,
                        TE_aid_plot,
                    )
                else:
                    handle_sequence_skipped(
                        seq_obj,
                        progress_file,
                        debug,
                        MSA_dir,
                        classification_dir,
                        plot_skip=plot_skip,
                        te_aid_plot=TE_aid_plot,
                        skip_proof_dir=skipped_dir,
                        orf_plot=input_orf_domain_plot,
                    )
                    click.echo(
                        f'\n{seq_name} was skipped because sequence is too short and check_low_copy is {check_low_copy}.\n'
                    )
                return

    except Exception as e:
        with open(error_files, 'a') as f:
            # Return the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(
                f'\nError during boundary finding and cropping for sequence: {seq_name}\n'
            )
            f.write(tb_content + '\n\n')
        prcyan(
            f'\nError during boundary finding and cropping for sequence: {seq_name}. Error: {str(e)}\n'
        )
        prcyan(tb_content + '\n')
        return

    # After all processing is done, change status to 'process' and write the file name to the progress file
    seq_obj.update_status('processed', progress_file)

    # If analysis of this sequence has been completed, remove all files contain sequence name
    if not debug:
        remove_files_with_start_pattern(MSA_dir, f'{seq_name}.fasta')

    # Read and count sequences from Finished_sequence_name.txt
    completed_sequence, skipped_count, low_copy_count, classified_pro = (
        check_progress_file(progress_file)
    )

    # Calculate the total count
    processed_count = len(completed_sequence)

    # Calculate number of sequences not processed by TE Trimmer
    rest_sequence = single_fasta_n - processed_count

    printProgressBar(
        processed_count,
        single_fasta_n,
        prefix='Progress:',
        suffix='Complete',
        length=50,
    )


def create_dir(
    continue_analysis, hmm, pfam_dir, output_dir, input_file, genome_file, plot_skip
):
    if not os.path.isfile(input_file):
        prcyan(
            f'The FASTA file {input_file} does not exist. Please check the input file path!'
        )
        raise FileNotFoundError
    input_file = os.path.abspath(input_file)  # get absolute path for input

    # check if genome file exist
    if not os.path.isfile(genome_file):
        prcyan(
            f'The genome FASTA file {genome_file} does not exist. Please check the genome file path!'
        )
        raise FileNotFoundError
    genome_file = os.path.abspath(genome_file)  # get absolute path for genome

    # bin_py_path contains all classes and BASH code
    # so.path.abspath(__file__) will return the current executable Python file
    bin_py_path = os.path.dirname(os.path.abspath(__file__))

    # Check if output path exists; otherwise create it
    os.makedirs(output_dir, exist_ok=True)
    output_dir = os.path.abspath(output_dir)  # get absolute path

    # Check if output directory is empty when --continue_analysis is 'False'
    if os.listdir(output_dir) and not continue_analysis:
        """
        prcyan(f"\nWARNING: The output directory {output_dir} is not empty. Please empty the output directory or "
               f"choose another empty directory.")
        prgre("\nNOTE: TE Trimmer can create output directory if it does not exist.")
        """
        # If the current folder is not empty, create a new folder with current time stamp
        current_time = time.strftime('%Y%m%d_%H%M%S')
        new_output_dir = os.path.join(output_dir, f'TEtrimmer_output_{current_time}')
        os.makedirs(new_output_dir, exist_ok=True)
        output_dir = new_output_dir
        prgre(
            f'\nThe given output directory is not empty. Results will be stored into folder: \n'
            f'{output_dir}\n'
        )

    # Create a new folder for single FASTA sequences
    single_file_dir = os.path.join(output_dir, 'Single_fasta_files')
    os.makedirs(single_file_dir, exist_ok=True)

    # Create a new folder for MSA
    MSA_dir = os.path.join(output_dir, 'Multiple_sequence_alignment')
    os.makedirs(MSA_dir, exist_ok=True)

    # Make a new folder for classification
    classification_dir = os.path.join(output_dir, 'Classification_and_deduplication')
    os.makedirs(classification_dir, exist_ok=True)

    # Create a new folder for HMM files
    if hmm:
        hmm_dir = os.path.join(output_dir, 'HMM_files')
        os.makedirs(hmm_dir, exist_ok=True)
    else:
        hmm_dir = ''

    # Define proof_curation folder path
    proof_curation_dir = os.path.join(output_dir, 'TEtrimmer_for_proof_curation')
    os.makedirs(proof_curation_dir, exist_ok=True)

    # Define clustered proof curation folder inside proof_curation_dir
    cluster_proof_anno_dir = os.path.join(
        proof_curation_dir, 'Clustered_proof_curation'
    )
    os.makedirs(cluster_proof_anno_dir, exist_ok=True)

    # Define skipped folder if required
    if plot_skip:
        skipped_dir = os.path.join(proof_curation_dir, 'TE_skipped')
        os.makedirs(skipped_dir, exist_ok=True)
    else:
        skipped_dir = None

    # Define low-copy folder
    low_copy_dir = os.path.join(proof_curation_dir, 'TE_low_copy')
    os.makedirs(low_copy_dir, exist_ok=True)

    # Define annotation evaluation folders
    perfect_proof = os.path.join(proof_curation_dir, 'Annotations_perfect')
    good_proof = os.path.join(proof_curation_dir, 'Annotations_good')
    intermediate_proof = os.path.join(
        proof_curation_dir, 'Annotations_check_recommended'
    )
    need_check_proof = os.path.join(proof_curation_dir, 'Annotations_check_required')
    os.makedirs(perfect_proof, exist_ok=True)
    os.makedirs(good_proof, exist_ok=True)
    os.makedirs(intermediate_proof, exist_ok=True)
    os.makedirs(need_check_proof, exist_ok=True)

    # Define the progress_file, finished sequence IDs will be stored here
    progress_file = os.path.join(output_dir, 'summary.txt')

    # Check and create progress_file if it does not exist yet
    if not os.path.exists(progress_file):
        with open(progress_file, 'a') as f:
            f.write(
                'input_name,consensus_name,blast_hit_n,cons_MSA_seq_n,cons_full_blast_n,input_length,cons_length,'
                'input_TE_type,reclassified_type,terminal_repeat,low_copy,evaluation,status\n'
            )

    # Define error files to store non-mandatory function errors including RepeatClassified classification,
    # RepeatMasker classification, PFAM scanning, MUSCLE alignment
    error_files = os.path.join(MSA_dir, 'error_file.txt')

    # If a PFAM database was not provided, create PFAM database in the TE Trimmer software folder, so the
    # database can be downloaded here.
    # If the pfam_dir was provided but the database cannot be found there, TE Trimmer will download the PFAM database
    # into the provided directory and generate the index file.
    if pfam_dir is None:
        pfam_dir = os.path.join(os.path.dirname(bin_py_path), 'pfam_database')
    try:
        os.makedirs(pfam_dir, exist_ok=True)
        if_pfam = prepare_pfam_database(pfam_dir)

        if not if_pfam:  # Check if if_pfam is 'False'
            raise Exception

    except Exception:
        with open(error_files, 'a') as f:
            # Return the traceback content as a string
            tb_content = traceback.format_exc()
            f.write('PFAM database building error\n')
            f.write(tb_content + '\n\n')

        prgre(
            'Note: Cannot download PFAM database from internet, please use a local PFAM database.\n'
            'For example: --pfam_dir <your_PFAM_directory>\n'
            'Your PFAM directory should contain: \n'
            'Pfam-A.hmm\n'
            'Pfam-A.hmm.h3f\n'
            'Pfam-A.hmm.h3m\n'
            'Pfam-A.hmm.dat\n'
            'Pfam-A.hmm.h3i\n'
            'Pfam-A.hmm.h3p\n\n'
            'PFAM predictions will be used to determine the direction of TEs. The database is therefore mandatory.\n\n'
            'You can download <Pfam-A.hmm.gz> and <Pfam-A.hmm.dat.gz> from '
            'https://www.ebi.ac.uk/interpro/download/pfam/\n'
            'Afterwards, do: \n'
            'gzip -d Pfam-A.hmm.gz\n'
            'gzip -d Pfam-A.hmm.dat.gz\n'
            'hmmpress Pfam-A.hmm\n\n'
        )
        return

    # Define consensus files. Temporary (temp) files will be used for the final RepeatMasker classification.
    final_con_file = os.path.join(output_dir, 'TEtrimmer_consensus.fasta')
    final_unknown_con_file = os.path.join(
        classification_dir, 'temp_TEtrimmer_unknown_consensus.fasta'
    )
    final_classified_con_file = os.path.join(
        classification_dir, 'temp_TEtrimmer_classified_consensus.fasta'
    )

    # Define consensus file without low copy
    final_con_file_no_low_copy = os.path.join(
        classification_dir, 'TEtrimmer_consensus_no_low_copy.fasta'
    )

    return (
        bin_py_path,
        output_dir,
        single_file_dir,
        MSA_dir,
        classification_dir,
        hmm_dir,
        proof_curation_dir,
        low_copy_dir,
        perfect_proof,
        good_proof,
        intermediate_proof,
        need_check_proof,
        progress_file,
        pfam_dir,
        final_con_file,
        final_con_file_no_low_copy,
        final_unknown_con_file,
        final_classified_con_file,
        error_files,
        input_file,
        genome_file,
        skipped_dir,
        cluster_proof_anno_dir,
    )
