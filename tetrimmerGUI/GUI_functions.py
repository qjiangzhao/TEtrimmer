import gzip
import logging
import os
import random
import shutil
import subprocess
import sys
import tarfile
import threading
import traceback
import urllib.request

import pandas as pd
import requests
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

special_character = {
    '/': '__',
    '\\': '__',
    ' ': '_',
    '|': '_',
    '<': '_',
    ':': '_',
    '"': '_',
    '*': '_',
    '?': '_',
}

def check_cmd_in_path(cmd):
    cmd_path = shutil.which(cmd)
    if cmd_path:
        #logging.info(f"{cmd} is available at: {cmd_path}")
        return True
    else:
        #logging.error(f"{cmd} is not available in the system's PATH")
        # Raise error if tool is not available
        #raise FileNotFoundError
        return False

def decompress_gzip(file_path):
    decompressed_file = file_path.rstrip('.gz')
    with gzip.open(file_path, 'rt') as f_in, open(decompressed_file, 'w') as f_out:
        shutil.copyfileobj(f_in, f_out)
    return decompressed_file

def get_original_file_path():
    """Get the path of the original script file, whether running in development or as a PyInstaller bundle."""
    if hasattr(sys, '_MEIPASS'):
        # Path to the temporary folder where PyInstaller extracts your package
        original_file_path = sys._MEIPASS
    else:
        # If the application is running in a normal Python environment, use __file__
        original_file_path = os.path.dirname(os.path.abspath(__file__))
    return original_file_path


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
                sanitized_id = f'{sanitized_id}_{count}'

            # Define output file name
            output_filename = os.path.join(output_dir, f'{sanitized_id}.fa')

            # Convert sequence name to sanitized_id
            record.id = sanitized_id
            record.description = ''

            # Write single FASTA file using sanitized name
            # the record.id is now same as sanitized_id
            with open(output_filename, 'w') as output_file:
                SeqIO.write(record, output_file, 'fasta')


def check_database(
        genome_file,
        output_dir=None,
        os_type='Darwin',
        use_system_blast=False
):
    #script_dir = get_original_file_path()
    #blast_dir = os.path.join(script_dir, 'blast')

    script_dir = get_original_file_path()

    blast_dir = os.path.join(script_dir, "blast")
    # Define database path
    genome_n = os.path.basename(genome_file)
    if output_dir is None:
        database_path = genome_file
    else:
        database_path = os.path.join(output_dir, genome_n)

    database_file = f'{database_path}.nin'

    if not os.path.isfile(database_file):
        try:

            # Check if the makeblastdb command is available in PATH
            # Use the user system makeblastdb when it exists, otherwise use the build-in version
            if use_system_blast and check_cmd_in_path('makeblastdb'):
                # Make makeblastdb command to create a BLAST database
                makeblastdb_cmd = f'makeblastdb -in {genome_file} -dbtype nucl -out {database_path}'
            else:

                if os_type == "Linux":
                    makeblastdb_cmd = f"{os.path.join(blast_dir, 'makeblastdb_linux')} -in {genome_file} -dbtype nucl -out {database_path}"
                elif os_type == "Darwin":
                    makeblastdb_cmd = f"{os.path.join(blast_dir, 'makeblastdb_mac')} -in {genome_file} -dbtype nucl -out {database_path}"
                else:
                    makeblastdb_cmd = f"{os.path.join(blast_dir, 'makeblastdb.exe')} -in {genome_file} -dbtype nucl -out {database_path}"

            subprocess.run(
                makeblastdb_cmd,
                shell=True,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            return database_path

        except FileNotFoundError:
            logging.error(
                f"'makeblastdb' command not found. Please ensure 'makeblastdb' is installed correctly.\n"
                f'{traceback.format_exc()}'
            )
            return 'makeblastdb_not_found'

        except subprocess.CalledProcessError as e:
            logging.error(
                f'makeblastdb failed with exit code {e.returncode} \n'
                f'{traceback.format_exc()}'
            )
            logging.error(e.stderr)
            return 'makeblastdb_got_error'
    else:
        return database_path


def blast(
    input_file,
    blast_database,
    blast_out_dir,
    e_value=1e-40,
    bed_file=False,
    self_blast=False,
    os_type='Darwin',
    num_threads=1,
    use_system_blast=False
):
    script_dir = get_original_file_path()
    blast_dir = os.path.join(script_dir, 'blast')

    # Define blast out without header
    blast_out_file = os.path.join(
        blast_out_dir, f'{os.path.basename(input_file)}_no_header.b'
    )
    # Define blast out with header
    blast_out_file_header = os.path.join(
        blast_out_dir, f'{os.path.basename(input_file)}_blast.txt'
    )

    # Check if the blastn command is available in PATH
    # Use the user system blastn when it exists, otherwise use the build-in version
    if use_system_blast and check_cmd_in_path('blastn'):
        blastn = "blastn"
    else:
        if os_type == "Linux":
            blastn = os.path.join(blast_dir, "blastn_linux")
        elif os_type == "Darwin":
            blastn = os.path.join(blast_dir, "blastn_mac")
        else:
            blastn = os.path.join(blast_dir, "blastn.exe")

    # -outfmt 6 use "\t" as deliminator. -outfmt 10 use "," as deliminator
    blast_cmd = [
        str(blastn),
        "-query", input_file,
        "-db", blast_database,
        "-max_target_seqs", str(10000),
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand evalue bitscore",
        "-evalue", str(e_value),
        "-out", blast_out_file,
        "-num_threads", str(num_threads)
    ]

    if self_blast:
        # Add more parameters for self-blast
        blast_cmd.extend([
            "-word_size", "11",
            "-gapopen", "5",
            "-gapextend", "2",
            "-reward", "2",
            "-penalty", "-3"
        ])

    try:
        # Run the BLAST command
        subprocess.run(blast_cmd, check=True, capture_output=True, text=True)

    except FileNotFoundError:
        logging.error(
            f"'blastn' command not found. Please ensure 'blastn' is installed correctly.\n"
            f'{traceback.format_exc()}'
        )
        return 'blastn_not_found', False

    except subprocess.CalledProcessError as e:
        logging.error(f'An error occurred during BLAST: \n {traceback.format_exc()}')
        logging.error(f'\nBLAST failed with error code {e.returncode}')
        # Print the error generated from the BLAST itself
        logging.error(e.stderr)
        return 'blastn_got_error', False

    # Check if the blast hit number is 0
    if os.path.isfile(blast_out_file) and os.path.getsize(blast_out_file) == 0:
        return 'blast_n_zero', False

    # TE-Aid analysis don't need bed file
    if not bed_file:
        # TE-Aid used different headers to deal with the blast output file.
        # Define the header for TE-Aid
        col_n = ["V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"]

        # Read the BLAST output file into a pandas DataFrame
        try:
            # Define blast output file for TE-Aid
            blast_out_teaid_file = f'{blast_out_file}_TEAid.b'

            # Add TE-Aid style header to the blast output
            blast_data_teaid_df = pd.read_csv(
                blast_out_file, delimiter='\t', header=None, names=col_n
            )

            # Save the DataFrame back to the file with the header, use "," as deliminator for TE-Aid
            blast_data_teaid_df.to_csv(blast_out_teaid_file, sep=',', index=False)

            return blast_out_teaid_file, False

        except Exception:
            logging.error(
                f'An error occurred while processing the BLAST output for TE-Aid : \n {traceback.format_exc()}'
            )
            return False, False
    # "Blast" button need bed file to allow to extract sequence from the genome
    else:
        try:
            # Add header to blast out
            # Define the header
            header = [
                'qseqid',
                'sseqid',
                'pident',
                'length',
                'mismatch',
                'gapopen',
                'qstart',
                'qend',
                'sstart',
                'send',
                'sstrand',
                'evalue',
                'bitscore',
            ]

            # Read the BLAST output into a DataFrame and add header to the blast output
            blast_out_file_df = pd.read_csv(
                blast_out_file, delimiter='\t', header=None, names=header
            )

            # Save the DataFrame back to blast_out_file_header file, this file will be shown in the GUI
            blast_out_file_df.to_csv(blast_out_file_header, sep='\t', index=False)
            blast_out_bed_file = f'{blast_out_file}.bed'

            blast_data_bed_df = pd.read_csv(blast_out_file, delimiter='\t', header=None)

            # Define a function to determine the strand and adjust columns accordingly
            def process_row(row):
                if 'plus' in row[10]:
                    return [row[1], row[8], row[9], row[0], 0, '+']
                else:
                    return [row[1], row[9], row[8], row[0], 0, '-']

            # Apply the function to each row and create a new DataFrame
            blast_data_bed_df_processed = blast_data_bed_df.apply(
                process_row, axis=1, result_type='expand'
            )

            # Drop duplicate rows
            blast_data_bed_df_processed_unique = (
                blast_data_bed_df_processed.drop_duplicates()
            )

            # Save bed file
            blast_data_bed_df_processed_unique.to_csv(
                blast_out_bed_file, sep='\t', index=False, header=False
            )

            return blast_out_bed_file, blast_out_file_header

        except Exception:
            logging.error(
                f'An error occurred while processing the BLAST output for BED file : \n {traceback.format_exc()}'
            )
            return False, False


def select_top_longest_lines(lines, n):
    return sorted(lines, key=lambda line: int(line[2]) - int(line[1]), reverse=True)[:n]


def select_random_lines(n, remaining_lines):
    random.shuffle(remaining_lines)
    return remaining_lines[:n]


def process_bed_lines(
    input_file, output_dir, max_lines=100, top_longest_lines_count=70
):
    """
    Select desired line numbers from bed file.

    :param output_dir: str, the absolute path of output folder directory.
    :param max_lines: num default=100, the maximum line number for bed file to keep.
    :param top_longest_lines_count: num default=100 smaller or equal to max_lines.
    When the bed file line number is excess than the max_lines, sort bed file by sequence length and choose
    "top_longest_lines_count" lines. When this number is smaller than the max_lines, randomly choose the rest
    number of lines from the bed file.

    :return: the absolute selected bed file path.
    """
    # check if top_longest_lines_count is equal or smaller than max_lines.
    if top_longest_lines_count > max_lines:
        top_longest_lines_count = max_lines

    with open(input_file, 'r') as file:
        lines = [line.strip().split('\t') for line in file]

    bed_out_filter_file = os.path.join(
        output_dir, f'{os.path.basename(input_file)}_f.bed'
    )

    if len(lines) > max_lines:
        top_longest_lines = select_top_longest_lines(lines, top_longest_lines_count)

        # Eliminate selected lines
        remaining_lines = [line for line in lines if line not in top_longest_lines]

        # Randomly choose lines
        random_lines = select_random_lines(
            max_lines - top_longest_lines_count, remaining_lines
        )
        selected_lines = top_longest_lines + random_lines

        # Write the selected lines to the output file.
        with open(bed_out_filter_file, 'w') as file:
            for line in selected_lines:
                file.write('\t'.join(line) + '\n')
    else:
        # Copy the original file to the output directory with the new name to keep the file name consistency.
        shutil.copy(input_file, bed_out_filter_file)

    return bed_out_filter_file


# Generate bed file based on fasta header，this is used for the extend function
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
def extend_bed_regions(
    bed_file, left_extension, right_extension, chrom_sizes, output_bed
):
    with open(bed_file, 'r') as infile, open(output_bed, 'w') as outfile:
        for line in infile:
            parts = line.strip().split()
            chrom, start, end, strand = parts[0], int(parts[1]), int(parts[2]), parts[5]

            # Adjust extensions based on strand
            if strand == '+':
                new_start = max(0, start - left_extension)
                new_end = min(chrom_sizes[chrom], end + right_extension)
            else:  # For anti sense strand
                new_start = max(
                    0, start - right_extension
                )  # Extend "start" less, as it's the "end"
                new_end = min(
                    chrom_sizes[chrom], end + left_extension
                )  # Extend "end" more, as it's the "start"

            # Write the extended region to the output BED file
            outfile.write(f'{chrom}\t{new_start}\t{new_end}\tTEtrimmer\t0\t{strand}\n')
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
                'strand': parts[5],
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
                    sequence = (
                        sequence.reverse_complement()
                    )  # Reverse complement if on negative strand
                seq_record = SeqRecord(
                    sequence, id=f'{chrom}:{start}-{end}({strand})', description=''
                )
                SeqIO.write(seq_record, out_fasta, 'fasta')
    return output_fasta


def align_sequences(input_file, output_file):
    """
    Aligns FASTA sequences using MAFFT

    :param input_file: str, input file path
    :param output_dir: str, output file directory
    :return: absolute path for multiple sequence alignment file
    """
    input_file_n = os.path.basename(input_file)

    # Construct the command as a list of strings
    mafft_cmd = ['mafft', '--quiet', '--nuc', '--retree', '1', input_file]

    try:
        # Execute the command
        result = subprocess.run(
            mafft_cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

    except FileNotFoundError as e:
        logging.error(
            "'mafft' command not found. Please ensure 'mafft' is installed correctly."
        )
        raise Exception from e

    except subprocess.CalledProcessError as e:
        error_message = e.stderr
        logging.error(f'\nMAFFT failed for {input_file_n} with error code {e.returncode}\n{e.stdout}\n{e.stderr}\n')

        if (
            'Killed' in error_message
            and 'disttbfast' in error_message
            and 'memopt' in error_message
        ):
            logging.error(
                'It seems insufficient RAM was available for MAFFT multiple sequence alignment. Please assign more '
                'RAM or reduce the thread number to solve this problem.\n'
            )
        raise Exception from e

    # Write the output to the file
    with open(output_file, 'w') as f:
        f.write(result.stdout)

    return output_file


def printProgressBar(
    iteration, total, prefix='', suffix='', decimals=1, length=50, fill='█', final=False
):
    """
    Custom terminal progress bar with MB display
    """
    percent = ('{0:.' + str(decimals) + 'f}').format(100 * (iteration / float(total)))
    if not final:
        percent = min(float(percent), 100)
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)

    # Convert bytes to MB for display
    iteration_mb = iteration / (1024 * 1024)
    total_mb = total / (1024 * 1024)

    print(
        f'\r{prefix} |{bar}| {iteration_mb:.2f}/{total_mb:.2f} MB = {percent}% {suffix}',
        end='',
        flush=True,
    )

    # Print a new line on completion
    if iteration == total:
        print()


def show_progress(block_num, block_size, total_size):
    if total_size > 0:
        iteration = block_num * block_size
        # Ensure we don't exceed total size
        iteration = min(iteration, total_size)
        printProgressBar(
            iteration,
            total_size,
            prefix='Downloading CDD:',
            suffix='Complete',
            length=50,
        )
    else:
        print('Total size unknown. Downloading without progress bar.')


def check_and_download(directory, check_pattern, filename, url):
    """
    Function to check if file exists, otherwise download and unzip it
    :param directory: str, directory the file to be checked should be in
    :param filename: str, file name used to store downloaded file
    :param url: str, url address used to download file if the file cannot be found in the given folder
    :return: boolean, 'True' means file was not found but was successfully downloaded. 'False' if file was found
    """

    check_pattern_file_path = os.path.join(directory, check_pattern)

    # Check if the cdd database is already downloaded but not unzipped
    if os.path.isfile(os.path.join(directory, filename)) and not os.path.isfile(check_pattern_file_path):
        logging.info(
            'CDD database is found but not unzipped. Unzipping......'
        )

        # Unzipping using tarfile
        with tarfile.open(os.path.join(directory, filename), 'r:gz') as tar:
            tar.extractall(path=directory)

        logging.info(
            f'\n{filename} is unzipped. CDD database was stored in \n'
            f'{directory}\n'
        )
        # Return True to indicate that the cdd database is found
        return True

    # If cdd database was not found, download it
    if not os.path.isfile(check_pattern_file_path) and not os.path.isfile(os.path.join(directory, filename)):
        logging.warning(
            'CDD database not found. Downloading... This might take some time. Please be patient.\n'
        )

        # Check if the URL is valid
        try:
            response = requests.get(url, stream=True)
            response.raise_for_status()  # Raise an exception if the GET request was unsuccessful

        except (
            requests.exceptions.HTTPError,
            requests.exceptions.ConnectionError,
        ) as e:
            logging.error(
                f'Failed to reach the server at {url} for downloading cdd database. {e}\n'
            )
            logging.info(
                'You can download the CDD database by yourself and refer to your database by '
                'the option --cdd_dir.'
            )
            return False

        try:
            gz_file_path = os.path.join(
                directory, filename
            )  # Provide a defined name for the file that will be downloaded
            urllib.request.urlretrieve(url, gz_file_path, show_progress)

            logging.info('\n CDD database is downloaded. Unzipping......')

            # Unzipping using tarfile
            with tarfile.open(gz_file_path, 'r:gz') as tar:
                tar.extractall(path=directory)

            # Delete gz file after extraction
            # os.remove(gz_file_path)

            logging.info(
                f'\n{filename} is downloaded and unzipped. CDD database was stored in \n'
                f'{directory}\n'
            )

        except Exception:
            logging.error('TEtrimmer failed to unpack the downloaded cdd database.')
            return False

    # Check if download was successful
    if os.path.isfile(check_pattern_file_path):
        logging.info('CDD database is found.')

        return True

    else:
        logging.error(
            'The CDD database cannot be downloaded by TEtrimmerGUI. Please check your internet connection '
            "or download CDD database manually and use '--cdd_dir' to indicate your cdd database path. "
            'After downloading, you have to unzip it by yourself. '
            "\n CDD database is only used to detect TE protein domains and doesn't affect other functions."
        )
        return False


def check_cdd_index_files(directory):
    """
    Check if the cdd index files exist in the provided directory.
    :param directory: str, the directory to search for the cdd index files.
    :return: boolean, 'True' if all cdd index files exist, 'False' otherwise.
    """
    # Define the list of cdd index files
    cdd_files = [
        'cdd_profile.24.freq',
        'cdd_profile.24.aux',
        'cdd_profile.24.rps',
        'cdd_profile.24.loo',
        'cdd_profile.24.ptf',
        'cdd_profile.24.pot',
        'cdd_profile.24.pdb',
        'cdd_profile.24.pos',
    ]

    # Check if all index files exist, store missing files in a list
    missing_files = []
    for file in cdd_files:
        if not os.path.isfile(os.path.join(directory, file)):
            missing_files.append(file)

    # Return True if all files exist, False otherwise
    if not missing_files:
        return True
    else:
        return False


def prepare_cdd_database(cdd_database_dir, os_type='Darwin', use_system_blast = False):
    # The downloaded database zipped name is cdd_tetrimmer.tar.gz
    # The name of cdd_database_dir is cdd_database, it is stored under TEtrimmerGUI path
    # Inside cdd_tetrimmer there should have a file called Cdd.cn, generate index file based that.
    # The index file is called Cdd_profile

    # cdd_database_dir is TEtrimmerGUI/cdd_database

    global prepared_cdd_g

    cdd_url = 'https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/cdd.tar.gz'
    cdd_pn_path = os.path.join(cdd_database_dir, 'Cdd.pn')

    try:
        # Check if CDD database downloads and unzipped in target dir
        if check_and_download(cdd_database_dir, r'Cdd.pn', r'Cdd.tar.gz', cdd_url):
            # Create index file for cdd database
            script_dir = get_original_file_path()
            blast_dir = os.path.join(script_dir, "blast")

            # Use the user system makeprofiledb when it exists, otherwise use the build-in version
            if use_system_blast and check_cmd_in_path('makeprofiledb'):
                makeprofiledb = "makeprofiledb"
            else:
                if os_type == "Linux":
                    makeprofiledb = os.path.join(blast_dir, "makeprofiledb_linux")
                elif os_type == "Darwin":
                    makeprofiledb = os.path.join(blast_dir, "makeprofiledb_mac")
                else:
                    makeprofiledb = os.path.join(blast_dir, "makeprofiledb.exe")

            makeprofiledb_cmd = [
                str(makeprofiledb),
                "-in", cdd_pn_path,
                "-out", os.path.join(cdd_database_dir, "cdd_profile"),
                "-threshold", str(9.82),
                "-scale", str(100.0),
                "-dbtype", "rps"]
            try:
                logging.info(
                    'CDD database index file not found. Making it ...... This can take around 10 mins.'
                )

                logging.info(' '.join(makeprofiledb_cmd))

                # Run the makeprofiledb command from the cdd_database_dir (with cwd).
                subprocess.run(
                    makeprofiledb_cmd,
                    check=True,
                    cwd=cdd_database_dir,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )

                if check_cdd_index_files(cdd_database_dir):
                    logging.info('CDD profile database index file generated.')
                else:
                    logging.error('CDD profile database index file generation failed.')
                    raise FileNotFoundError

                # Change global variable prepared_cdd_g when it is generated successfully
                prepared_cdd_g = os.path.join(cdd_database_dir, 'cdd_profile')
                return os.path.join(cdd_database_dir, 'cdd_profile')

            except FileNotFoundError:
                logging.error("'makeprofiledb' command not found.")
                return False

            except subprocess.CalledProcessError as e:
                logging.error(
                    f'\nCDD index file generation failed with error code {e.returncode}\n{e.stdout}\n{e.stderr}')
                return False
        else:
            return False

    except Exception as e:
        logging.error(f'Prepare CDD database failed with error {e}')
        return False


def rpstblastn(
    input_file,
    rpsblast_database,
    rpsblast_out_dir,
    e_value=0.01,
    os_type='Darwin',
    num_threads=20,
    use_system_blast = False
):
    script_dir = get_original_file_path()
    blast_dir = os.path.join(script_dir, 'blast')

    # Define output file
    rpstblastn_out_file = os.path.join(
        rpsblast_out_dir, f'{os.path.basename(input_file)}_pcc.txt'
    )

    # Check if the rpstblastn command is available in PATH
    # Use the user system rpstblastn when it exists, otherwise use the build-in version
    if use_system_blast and check_cmd_in_path('rpstblastn'):
        rpstblastn = "rpstblastn"
    else:
        if os_type == "Linux":
            rpstblastn = os.path.join(blast_dir, "rpstblastn_linux")
        elif os_type == "Darwin":
            rpstblastn = os.path.join(blast_dir, "rpstblastn_mac")
        else:
            rpstblastn = os.path.join(blast_dir, "rpstblastn.exe")

    # -outfmt 6 use "\t" as deliminator.
    rpsblast_cmd = [
        str(rpstblastn),
        "-query", input_file,
        "-db", rpsblast_database,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue score stitle",
        "-evalue", str(e_value),
        "-out", rpstblastn_out_file,
        "-num_threads", str(num_threads)
    ]

    try:
        # Print the rpstblastn command
        logging.info(' '.join(rpsblast_cmd))
        # Run the rpstblastn command
        subprocess.run(rpsblast_cmd, check=True, capture_output=True, text=True)


        # Check if the blast hit number is 0
        if (
            os.path.isfile(rpstblastn_out_file)
            and os.path.getsize(rpstblastn_out_file) == 0
        ):
            logging.warning('rpstblastn hit number is 0.')
            return 'rpstblastn_n_zero'

        return rpstblastn_out_file

    except FileNotFoundError:
        logging.error(f"'rpsblastn' command not found.\n{traceback.format_exc()}")
        return False

    except subprocess.CalledProcessError as e:
        logging.error(f'An error occurred during rpstblastn: \n {traceback.format_exc()}')
        logging.error(f'\nrpstblastn failed with error code {e.returncode}')
        # Print the error generated from the BLAST itself
        logging.error(e.stderr)
        return False


def run_func_in_thread(func, *args, **kwargs):
    """
    Runs a given function in a new thread with the provided arguments.

    Parameters:
        func (callable): The function to run in a thread.
        *args: Positional arguments to pass to the function.
        **kwargs: Keyword arguments to pass to the function.

    Returns:
        threading.Thread: The thread object running the function.
    """
    # Create and start a new thread for the function
    thread = threading.Thread(target=func, args=args, kwargs=kwargs)
    thread.start()
    return thread  # Optionally return the thread for monitoring


class CustomFormatter(logging.Formatter):
    """
    Custom logging formatter with optional color and per-level formatting.
    """

    # ANSI escape codes for terminal colors
    grey = "\x1b[38;21m"
    black = "\x1b[30m"
    blue = "\x1b[38;5;39m"
    orange = "\x1b[38;5;214m"
    red = "\x1b[38;5;196m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"

    def __init__(self, use_color=True):
        super().__init__()
        self.use_color = use_color
        self.datefmt = "%H:%M:%S"  # Short timestamp format

        # Two format styles
        self.simple_fmt = "%(asctime)s | %(levelname)s | %(message)s"
        self.detailed_fmt = "%(asctime)s | %(levelname)s | %(module)s | %(funcName)s | %(lineno)d | %(message)s"

        # Color and plain format mappings
        self.color_formats = {
            logging.DEBUG: self.grey + self.detailed_fmt + self.reset,
            logging.INFO: self.black + self.simple_fmt + self.reset,
            logging.WARNING: self.orange + self.simple_fmt + self.reset,
            logging.ERROR: self.red + self.detailed_fmt + self.reset,
            logging.CRITICAL: self.bold_red + self.detailed_fmt + self.reset,
        }

        self.plain_formats = {
            logging.DEBUG: self.detailed_fmt,
            logging.INFO: self.simple_fmt,
            logging.WARNING: self.simple_fmt,
            logging.ERROR: self.detailed_fmt,
            logging.CRITICAL: self.detailed_fmt,
        }

    def format(self, record):
        formats = self.color_formats if self.use_color else self.plain_formats
        log_fmt = formats.get(record.levelno)
        formatter = logging.Formatter(log_fmt, datefmt=self.datefmt)
        return formatter.format(record)


def init_logging(loglevel="DEBUG", logfile=None):
    """
    Initialize logging with colored console output and optional file logging.

    Parameters:
    - loglevel: str (e.g., "DEBUG", "INFO")
    - logfile: str or None (if provided, logs will also be written to file)
    """
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {loglevel}")

    # Console handler (with color)
    handler_console = logging.StreamHandler(sys.stderr)
    handler_console.setFormatter(CustomFormatter(use_color=True))

    handlers = [handler_console]

    # Optional file handler (no color)
    if logfile is not None:
        handler_file = logging.FileHandler(logfile)
        handler_file.setFormatter(CustomFormatter(use_color=False))
        handlers.append(handler_file)

    logging.basicConfig(level=numeric_level, handlers=handlers)
