import shutil
import subprocess
import os
import os.path
import click
import shutil
import random
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment, AlignInfo
from Bio.motifs import Motif
import pandas as pd
import pandas.errors
import numpy as np
from PyPDF2 import PdfMerger, PdfFileReader, PdfFileWriter

import warnings
from Bio import BiopythonDeprecationWarning

# Suppress all deprecation warnings
warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)


# Check if file name contains the string LTR
def is_LTR(input_file):
    input_file_name = os.path.basename(input_file)
    return 'LTR' in input_file_name


# Generate consensus sequence
def con_generater_no_file(input_file, threshold=0.8, ambiguous="N"):
    # Read input file
    alignment = AlignIO.read(input_file, "fasta")

    # Generate a summary of the alignment
    summary_align = AlignInfo.SummaryInfo(alignment)

    # Calculate consensus sequence with the specified threshold
    consensus = summary_align.dumb_consensus(threshold=threshold, ambiguous=ambiguous).upper()

    return consensus


def con_generater(input_file, output_dir, threshold=0.8, ambiguous="N"):
    # Read input file
    alignment = AlignIO.read(input_file, "fasta")

    # Generate a summary of the alignment
    summary_align = AlignInfo.SummaryInfo(alignment)

    # Calculate consensus sequence with the specified threshold
    consensus = summary_align.dumb_consensus(threshold=threshold, ambiguous=ambiguous)

    # Create SeqRecord for consensus sequence
    consensus_record = SeqRecord(consensus, id=os.path.basename(input_file), description="")

    # Write consensus sequence to a FASTA file
    output_file = os.path.join(output_dir, f"{os.path.basename(input_file)}_co.fa")
    with open(output_file, "w") as file:
        SeqIO.write(consensus_record, file, "fasta")

    return output_file


# Check if start and end are matchable with the given pattern
def check_start_end(consensus_seq, start, end, start_patterns, end_patterns):
    start_matched = end_matched = True

    if start_patterns:
        for start_pattern in start_patterns:
            if consensus_seq[start: start + len(start_pattern)] == list(start_pattern):
                break
        else:
            start_matched = False

    if end_patterns:
        for end_pattern in end_patterns:
            if consensus_seq[end - len(end_pattern): end] == list(end_pattern):
                break
        else:
            end_matched = False

    return start_matched, end_matched


def check_and_update(consensus_seq, start, end, start_patterns, end_patterns):
    """
    Check if start and end of the MSA are equal to the given patterns
    """
    # Ensure all patterns are uppercase
    start_patterns = [pattern.upper() for pattern in start_patterns.split(',')] if start_patterns else None  
    end_patterns = [pattern.upper() for pattern in end_patterns.split(',')] if end_patterns else None

    start_matched, end_matched = check_start_end(consensus_seq, start, end, start_patterns, end_patterns)

    if not start_matched or not end_matched:
        # If exact position matching fails, try sliding window approach
        if not start_matched and start_patterns:
            for start_pattern in start_patterns:
                start_window = consensus_seq[max(0, start - 15): start + 15]
                for i in range(len(start_window) - len(start_pattern) + 1):
                    if start_window[i: i + len(start_pattern)] == list(start_pattern):
                        start = max(0, start - 15) + i
                        break
                break

        if not end_matched and end_patterns:
            for end_pattern in end_patterns:
                end_window = consensus_seq[max(0, end - 15): end + 15]
                for i in reversed(range(len(end_window) - len(end_pattern) + 1)):
                    if end_window[i: i + len(end_pattern)] == list(end_pattern):
                        end = max(0, end - 15) + i + len(end_pattern)
                        break
                break

    # Check again after updating positions
    start_matched, end_matched = check_start_end(consensus_seq, start, end, start_patterns, end_patterns)
    return start_matched, end_matched, start, end


def prcyan(text):
    click.echo(click.style(text, fg='cyan'))


def prgre(text):
    click.echo(click.style(text, fg='green'))


def fasta_file_to_dict(input_file, separate_name=False):
    sequences = {}
    for record in SeqIO.parse(input_file, "fasta"):
        if separate_name:
            sequences[record.id.split('#')[0]] = record
        else:
            sequences[record.id] = record
    return sequences


def blast(seq_file, genome_file, output_dir, min_length=150, search_type="blast", task="blastn", seq_obj=None):
    """
    Runs BLAST calling a specified task type and saves the results as a BED file.

    :param seq_file: str, path to input FASTA file
    :param genome_file: str, path to genome FASTA file
    :param output_dir: str, prefix for output files
    :param min_length: int, minimum alignment length. Default: 150
    :param task: str, BLAST task type ("blastn", "dc-megablast", etc.). Default: "blastn"
    :param seq_obj: object (optional), sequence object to update with BLAST hits
    :return: tuple, output BED file name and BLAST hit count for each sequence
    """
    input_file = seq_file
    input_file_n = os.path.basename(input_file)
    blast_hits_count = 0
    bed_out_file = None
    # define blast outfile
    blast_out_file = os.path.join(output_dir, f"{os.path.basename(input_file)}.b")

    if search_type == "blast":
        # Modify the blast command to include the specified task
        blast_cmd = (f"blastn -task {task} -query {input_file} -db {genome_file} "
                     f"-outfmt \"6 qseqid sseqid pident length mismatch qstart qend sstart send sstrand evalue qcovhsp\" "
                     f"-evalue 1e-40 -qcov_hsp_perc 20 | "
                     f"awk -v ml={min_length} 'BEGIN{{OFS=\"\\t\"}} $4 > ml {{print $0}}' >> {blast_out_file}")
        try:
            subprocess.run(blast_cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        except FileNotFoundError:
            prcyan("'blastn' command not found. Please ensure 'blastn' is installed correctly.")
            raise Exception

        except subprocess.CalledProcessError as e:
            prcyan(f"\nBLAST failed for {input_file_n} with error code {e.returncode}")
            prcyan(e.stderr)
            raise Exception

    elif search_type == "mmseqs":
        # MMseqs2 search
        mmseqs_cmd = (f"mmseqs easy-search {seq_file} {genome_file}_db {blast_out_file} {blast_out_file}_tmp --search-type 3 "
                      f"--min-seq-id 0.6 --format-output \"query,target,pident,alnlen,mismatch,qstart,qend,tstart,tend,evalue,qcov\" "
                      f"--cov-mode 4 -c 0.5 --e-profile 1e-40 --threads 1 --min-aln-len {min_length}")
        print(f"Running analysis with MMseqs2. \nWARNING: MMseqs2 is less accurate than 'blastn' and may return faulty results. We strongly recommend to use 'blastn' instead.")
        result = subprocess.run(mmseqs_cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        error_output = result.stderr.decode("utf-8")

        if error_output:
            print(f"Error running MMseqs2: {error_output}")

    # Check the number of BLAST hits
    with open(blast_out_file) as blast_file:
        for _ in blast_file:
            blast_hits_count += 1

    if blast_hits_count > 0:
        # Define BED outfile
        # add $4 alignment length
        bed_out_file = os.path.join(output_dir, f"{os.path.basename(input_file)}.b.bed")

        if search_type == "blast":
            bed_cmd = (f"awk 'BEGIN{{OFS=\"\\t\"; counter=0}} !/^#/ {{counter+=1; "
                       f"if ($10~/plus/){{print $2, $8, $9, counter, $3, \"+\", $4, $1}} "
                       f"else {{print $2, $9, $8, counter, $3, \"-\", $4, $1}}}}' < {blast_out_file} > {bed_out_file}")
        elif search_type == "mmseqs":
            bed_cmd = (f"awk 'BEGIN{{OFS=\"\\t\"; counter=0}} !/^#/ {{counter+=1; "
                       f"if ($7>$6){{print $2, $8, $9, counter, $3, \"+\", $4, $1}} "
                       f"else {{print $2, $8, $9, counter, $3, \"-\", $4, $1}}}}' < {blast_out_file} > {bed_out_file}")
        try:
            subprocess.run(bed_cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        except subprocess.CalledProcessError as e:
            prcyan(f"\nBLAST to BED file conversion failed for {input_file_n} with error code {e.returncode}")
            prcyan(e.stderr)
            raise Exception

        # Append BLAST hit number to sequence object when seq_obj is provided
        if seq_obj is not None:
            seq_obj.update_blast_hit_n(blast_hits_count)

    return bed_out_file, blast_hits_count, blast_out_file


def check_bed_uniqueness(output_dir, bed_file):
    """
    Check if the BED file contains duplicated lines, and if so, delete any duplicates

    :param output_dir: str, output directory
    :param bed_file: str, path for input BED file
    :return: absolute path of modified BED file
    """
    with open(bed_file, "r") as f:
        lines = f.readlines()
    entries = [line.strip().split("\t") for line in lines]
    unique_keys = set()

    unique_entries = []
    for entry in entries:
        key = "".join(entry[:3])
        if key not in unique_keys:
            unique_keys.add(key)
            unique_entries.append(entry)

    bed_out_file = os.path.join(output_dir, f"{os.path.basename(bed_file)}_u")

    with open(bed_out_file, "w") as file:
        for entry in unique_entries:
            line = "\t".join(entry)
            file.write(f"{line}\n")
    return bed_out_file

def read_bed_file(input_file):
    with open(input_file, 'r') as file:
        return [line.strip().split('\t') for line in file]


def remove_duplicates(lines):
    unique_lines = []
    line_set = set()
    for line in lines:
        line_key = tuple(line[:3])
        if line_key not in line_set:
            line_set.add(line_key)
            unique_lines.append(line)

    return unique_lines


def select_top_longest_lines(lines, n):
    return sorted(lines, key=lambda line: int(line[2]) - int(line[1]), reverse=True)[:n]


def select_random_lines(n, remaining_lines):
    random.shuffle(remaining_lines)
    return remaining_lines[:n]


def process_lines(input_file, output_dir, threshold=100, top_longest_lines_count=100):
    """
    Select desired line numbers from bed file.

    :param output_dir: str, the absolute path of output folder directory.
    :param threshold: num default=100, the maximum line number for bed file to keep.
    :param top_longest_lines_count: num default=100 smaller or equal to threshold.
    When the bed file line number is excess than the threshold, sort bed file by sequence length and choose
    "top_longest_lines_count" lines. When this number is smaller than the threshold, randomly choose the rest
    number of lines from the bed file.

    :return: the absolute selected bed file path.
    """
    # check if top_longest_lines_count is equal or smaller than threshold.
    if top_longest_lines_count > threshold:
        raise ValueError("top_longest_lines must be equal to or smaller than threshold.")
    lines = read_bed_file(input_file)
    lines = remove_duplicates(lines)
    bed_out_filter_file = os.path.join(output_dir, f"{os.path.basename(input_file)}f.bed")

    if len(lines) > threshold:
        top_longest_lines = select_top_longest_lines(lines, top_longest_lines_count)
        # Eliminate selected lines
        remaining_lines = [line for line in lines if line not in top_longest_lines]
        # Randomly choose lines
        random_lines = select_random_lines(threshold - top_longest_lines_count, remaining_lines)
        selected_lines = top_longest_lines + random_lines
        # Write the selected lines to the output file.
        with open(bed_out_filter_file, 'w') as file:
            for line in selected_lines:
                file.write("\t".join(line) + "\n")
    else:
        # Copy the original file to the output directory with the new name to keep the file name consistency.
        shutil.copy(input_file, bed_out_filter_file)

    return bed_out_filter_file

# Calculate average sequence length in the bed file
def bed_ave_sequence_len(bed_content, start_rank, end_rank):
    """Compute average length for regions within a specified rank range in BED file."""

    # Extracting lengths from the BED file
    lengths = []
    for line in bed_content.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) >= 3:
            start, end = int(parts[1]), int(parts[2])
            lengths.append(end - start)

    # Sorting lengths in descending order
    lengths.sort(reverse=True)

    # Checking if there are enough entries for the specified range
    if len(lengths) < end_rank:
        return "Not enough entries in BED file."

    # Extracting lengths of regions within the specified rank range
    selected_lengths = lengths[start_rank - 1:end_rank]

    # Calculating average
    average = sum(selected_lengths) / len(selected_lengths)

    return average


def extract_fasta(input_file, genome_file, output_dir, left_ex, right_ex, nameonly=False):
    """
    Extracts FASTA sequence from the reference genome using BEDtools.

    :param genome_file: str, path to genome FASTA file
    :param output_dir: str, prefix for output files
    :param left_ex: int, number of bases to extend the start position of each feature
    :param right_ex: int, number of bases to extend the end position of each feature
    :param nameonly: boolean, only use BED file name for the FASTA header. Default: False
    :return: absolute path of FASTA file derived from BED file
    """
    input_file_n = os.path.basename(input_file)

    bed_out_flank_file_dup = os.path.join(output_dir, f"{os.path.basename(input_file)}_{left_ex}_{right_ex}.bed")

    fasta_out_flank_file = os.path.join(output_dir, f"{os.path.basename(input_file)}_{left_ex}_{right_ex}.fa")

    fasta_out_flank_file_nucleotide_clean = os.path.join(output_dir,
                                                         f"{os.path.basename(input_file)}_{left_ex}_{right_ex}_bcl.fa")

    if nameonly:
        bed_out_flank_file_dup = os.path.join(output_dir, f"{os.path.basename(input_file)}_{left_ex}_{right_ex}_n.bed")

        fasta_out_flank_file = os.path.join(output_dir, f"{os.path.basename(input_file)}_{left_ex}_{right_ex}_n.fa")

        fasta_out_flank_file_nucleotide_clean = os.path.join(output_dir,
                                                             f"{os.path.basename(input_file)}_{left_ex}_{right_ex}_bcln.fa")

    bed_cmd = f"bedtools slop -s -i {input_file} -g {genome_file}.length -l {left_ex} -r {right_ex} > {bed_out_flank_file_dup}"

    try:
        subprocess.run(bed_cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    except FileNotFoundError:
        prcyan("'bedtools slop' command not found. Please ensure 'bedtools' is installed correctly.")
        raise Exception

    except subprocess.CalledProcessError as e:
        prcyan(f"\nbedtools slop failed for {input_file_n} with error code {e.returncode}")
        prcyan(e.stderr)
        raise Exception

    bed_out_flank_file = check_bed_uniqueness(output_dir, bed_out_flank_file_dup)

    if nameonly:
        fasta_cmd = f"bedtools getfasta -s -nameOnly -fi {genome_file} -fo {fasta_out_flank_file} -bed {bed_out_flank_file}"
    else:
        fasta_cmd = f"bedtools getfasta -s -fi {genome_file} -fo {fasta_out_flank_file} -bed {bed_out_flank_file}"

    try:
        subprocess.run(fasta_cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    except FileNotFoundError:
        prcyan("'bedtools getfasta' command not found. Please ensure 'bedtools' is installed correctly.")
        raise Exception

    except subprocess.CalledProcessError as e:
        prcyan(f"\nbedtools getfasta failed for {input_file_n} with error code {e.returncode}")
        prcyan(e.stderr)
        raise Exception

    # Use awk to remove letters other than: A G C T a g c t
    fasta_nucleotide_clean = f"awk '/^>/ {{print}} !/^>/ {{gsub(/[^AGCTagct]/, \"\"); print}}' {fasta_out_flank_file}" \
                             f" > {fasta_out_flank_file_nucleotide_clean}"
    try:
        subprocess.run(fasta_nucleotide_clean, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    except subprocess.CalledProcessError as e:
        prcyan(e.stderr)
        raise Exception

    return fasta_out_flank_file_nucleotide_clean, bed_out_flank_file_dup


def align_sequences(input_file, output_dir):
    """
    Aligns FASTA sequences using MAFFT

    :param input_file: str, input file path
    :param output_dir: str, output file directory
    :return: absolute path for multiple sequence alignment file
    """
    input_file_n = os.path.basename(input_file)
    fasta_out_flank_mafft_file = os.path.join(output_dir, f"{os.path.basename(input_file)}_aln.fa")

    # Read the top 5 sequences and determine if any are longer than 10000 bp
    long_sequences = False
    with open(input_file, "r") as f:
        for i, record in enumerate(SeqIO.parse(f, "fasta")):
            if len(record.seq) > 10000:
                long_sequences = True
                break
            if i >= 10:  # Only check the top 10 sequences
                break

    # Construct the command as a list of strings
    mafft_cmd = ["mafft", "--quiet", "--nuc", "--retree", "1", input_file]

    # If any of the top 5 sequences are longer than 10000, add --memsave to save memory
    if long_sequences:
        mafft_cmd.insert(1, "--memsave")  # Insert after 'mafft'

    try:
        # Execute the command
        result = subprocess.run(mafft_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    except FileNotFoundError:
        prcyan("'mafft' command not found. Please ensure 'mafft' is installed correctly.")
        raise Exception

    except subprocess.CalledProcessError as e:
        error_message = e.stderr
        prcyan(f"MAFFT failed for {input_file_n} with error code {e.returncode}\n")
        prcyan(e.stderr)

        if "Killed" in error_message and "disttbfast" in error_message and "memopt" in error_message:
            prgre(f"It seems insufficient RAM was available for MAFFT multiple sequence alignment. Please assign more "
                  f"RAM or reduce the thread number to solve this problem.\n")
        raise Exception

    # Write the output to the file
    with open(fasta_out_flank_mafft_file, 'w') as f:
        f.write(result.stdout)

    return fasta_out_flank_mafft_file


def muscle_align(input_file, output_dir, ite_times=4):
    output_file = os.path.join(output_dir, f"{os.path.basename(input_file)}_maln.fa")
    muscle_cmd = ["muscle", "-maxiters", str(ite_times), "-in", input_file, "-out", output_file]

    # Execute MUSCLE
    result = subprocess.run(muscle_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if result.returncode != 0:
        return False
        # raise Exception(f"MUSCLE command failed with error: {os.path.basename(input_file)}\n{result.stderr.decode('utf-8')}")
    # Convert the sequences in the output file to lowercase
    sequences = list(SeqIO.parse(output_file, "fasta"))
    for seq in sequences:
        seq.seq = seq.seq.lower()
    SeqIO.write(sequences, output_file, "fasta")

    return output_file


def calc_proportion(input_file):
    # Read input file
    alignment = AlignIO.read(input_file, "fasta")
    max_props = []

    # Loop through each column of the alignment
    for i in range(alignment.get_alignment_length()):
        counts = {"a": 0, "c": 0, "g": 0, "t": 0}
        for record in alignment:
            nucleotide = record.seq[i].lower()
            if nucleotide in counts:
                counts[nucleotide] += 1

        # Calculate the proportion of each nucleotide at each site
        total = sum(counts.values())

        # Check if total is zero
        if total == 0:
            max_props.append(0)
        else:
            proportions = {nucleotide: count / total for nucleotide, count in counts.items()}
            max_props.append(max(proportions.values()))  # Append the maximum proportion

    return np.array(max_props)  # Convert the list to a NumPy array and return it


# Function to help define the threshold for crop end by similarity
def define_crop_end_simi_thr(input_file, window_size=40, max_steps=100):
    array = calc_proportion(input_file)
    array_length = len(array)
    start_sum = 0
    end_sum = 0
    steps = 0

    # Determine the range for sliding windows
    if array_length < 2 * max_steps:
        iteration_range = range(array_length // 2)
    else:
        iteration_range = range(max_steps)

    # Calculate sum for start and end sliding windows
    for i in iteration_range:
        # Start window
        start_window = array[i:i + window_size]
        start_sum += np.sum(start_window)

        # End window
        end_window = array[-(i + window_size): -i if i != 0 else None]
        end_sum += np.sum(end_window)
        steps += 1
    mean_start_sum = start_sum / steps
    mean_end_sum = end_sum / steps
    # Calculate mean of the sum of proportions
    # mean_sum = (start_sum + end_sum) / (steps * 2)
    return mean_start_sum, mean_end_sum


def calc_conservation(col):
    """
    Calculate the conservation of a column as the fraction of the most common nucleotide.
    :param col: one single column content for multiple sequence alignment
    :return: most abundant nucleotide proportion for this column (excluding gaps)
    """
    # MAFFT will convert all nucleotide to lowercase
    nucleotide_counts = {'a': 0, 'c': 0, 'g': 0, 't': 0}

    for nucleotide in col:
        if nucleotide in nucleotide_counts:
            nucleotide_counts[nucleotide] += 1

    total_nucleotides = sum(nucleotide_counts.values())
    max_count = max(nucleotide_counts.values())

    return max_count / total_nucleotides


def generate_hmm_from_msa(input_msa_file, output_hmm_file, error_file):
    """
    Generate HMM profiles using hmmbuild on a multiple sequence alignment file.

    :param input_msa_file: path to the multiple sequence alignment file (in FASTA or Stockholm format)
    :param output_hmm_file: path to the output HMM file
    """

    # Construct the command as a list
    cmd = ["hmmbuild", "--dna", output_hmm_file, input_msa_file]

    try:
        # Execute the command
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    except FileNotFoundError:
        prcyan("'hmmbuild' command not found. Please ensure 'hmmbuild' is installed correctly.")
        prgre("This hmm error will not affect the final TE consensus library.")

    except subprocess.CalledProcessError as e:
        with open(error_file, 'a') as f:
            f.write(f"\nhmm file generation failed for {os.path.basename(output_hmm_file)} with error code {e.returncode}")
            f.write('\n' + e.stderr)
        prcyan(f"\nhmm file generation failed for {os.path.basename(output_hmm_file)} with error code {e.returncode}")
        prgre("\nThis hmm error will not affect the final TE consensus library, you can ignore it."
              "\nFor traceback text, please refer to 'error_file.txt' in the 'Multiple_sequence_alignment' folder.\n")


def reverse_complement_seq_file(input_file, output_file):
    """
    Takes an MSA FASTA file and writes the reverse complement of the sequences to a new file.
    """
    reverse_complemented_records = []

    for record in SeqIO.parse(input_file, "fasta"):
        # Create a new record with the reverse complemented sequence
        rev_comp_seq = record.seq.reverse_complement()
        rev_comp_record = SeqRecord(rev_comp_seq, id=record.id, description="")
        reverse_complemented_records.append(rev_comp_record)

    # Write the reverse-complemented sequences to the output file
    SeqIO.write(reverse_complemented_records, output_file, "fasta")
    return output_file


def remove_short_seq(input_file, output_file, threshold=0.1):
    alignment = AlignIO.read(input_file, "fasta")

    # Get MSA length
    len = alignment.get_alignment_length()

    # Initialize an empty list to store the filtered sequences
    filtered_alignment = []

    # Loop through each sequence in the alignment
    for record in alignment:
        # Get the length of the sequence, excluding gaps
        seq_length = len(record.seq.ungap("-"))

        # Check if the sequence length is greater than or equal to the length threshold
        if seq_length >= len * threshold:
            filtered_alignment.append(record)

    # Write the filtered list record to a file
    AlignIO.write(filtered_alignment, output_file)


def remove_gaps(input_file, output_dir, threshold=0.8, min_nucleotide=5):
    """
    Remove gaps if the gap percentage is above threshold. Remove columns if nucleotide number is less than 5.
    :param input_file: str, input file path
    :param output_dir: str, output file directory
    :param threshold: int (0-1), columns with gap percentage greater than <threshold> will be removed. Default: 0.8
    :param min_nucleotide: int (>0), columns with less than <min_nucleotide> nucleotides will be removed. Default: 5
    :return: absolut path for the gap-purged multiple sequence alignment file
    """
    keep_list = []
    fasta_out_flank_mafft_gap_filter_file = os.path.join(output_dir,
                                                         f"{os.path.basename(input_file)}_g.fa")
    MSA_mafft = AlignIO.read(input_file, "fasta")

    column_mapping = {}  # Stores the mapping of column indices of filtered MSA to original MSA

    for col_idx in range(MSA_mafft.get_alignment_length()):
        col = MSA_mafft[:, col_idx]
        gap_count = col.count('-')
        gap_fraction = gap_count / len(col)

        # Check if the count of nucleotides in the column is 5 or greater
        nt_count = len(col) - gap_count
        if nt_count < min_nucleotide:
            continue

        # If the gap fraction is less than the threshold, add the column to the filtered MSA
        if gap_fraction <= threshold:
            keep_list.append(col_idx)
            # Be careful when using len for indexing, because len starts with 1
            column_mapping[
                len(keep_list) - 1] = col_idx  # Store the mapping of original MSA index to filtered MSA index

    # The index in Python does not include the second value. As a consequence, end_posit from the DefineBoundary()
    # function contains one more index than the final [len(keep_list) -1]. 
    # Adding one more value to correct index
    column_mapping[len(keep_list)] = column_mapping[len(keep_list) - 1] + 1
    # Keep the columns if they meet requirements
    MSA_mafft_filtered = MSA_mafft[:, keep_list[0]:keep_list[0] + 1]
    for i in keep_list[1:]:
        MSA_mafft_filtered += MSA_mafft[:, i:i + 1]

    # Write the filtered MSA to the output file
    with open(fasta_out_flank_mafft_gap_filter_file, 'w') as f:
        AlignIO.write(MSA_mafft_filtered, f, 'fasta')
    return fasta_out_flank_mafft_gap_filter_file, column_mapping


def remove_gaps_block(input_file, output_dir, threshold=0.8, conservation_threshold=0.5, min_nucleotide=5):
    """
    Remove gap blocks flanked by conserved regions of the MSA. Remove columns if nucleotide number is less than 5.
    :param input_file: str, input file path
    :param output_dir: str, output file directory
    :param threshold: int (0-1) columns with gap percentage greater than <threshold> will be removed. Default: 0.8
    :param conservation_threshold: int (0-1), after connected gaps have been classified as a gap block, the nucleotide
    conservation of two flanking non-gap columns at both ends of the gap block will be calculated. If the proportion of
    the most abundant nucleotide is greater than <conservation_threshold>, this gap block will be removed from the MSA. 
    Default: 0.5
    :param min_nucleotide: int (>0), columns with less than <min_nucleotide> nucleotides will be removed. Default: 5
    :return: absolute path to gap block-purged MSA FASTA file
    """
    # The keep_list is a list of tuples containing start and end position for each gap block
    gap_blocks = []
    fasta_out_flank_mafft_gap_filter_file = os.path.join(output_dir,
                                                         f"{os.path.basename(input_file)}_gb.fa")
    MSA_mafft = AlignIO.read(input_file, "fasta")

    # Initialize the start of the block as 'None'
    block_start = None

    for col_idx in range(MSA_mafft.get_alignment_length()):
        col = MSA_mafft[:, col_idx]
        gap_count = col.count('-')
        gap_fraction = gap_count / len(col)
        # Check the count of nucleotides in the column
        nt_count = len(col) - gap_count

        # If the column meets the gap requirement and is not within an already called block, start a new block
        if (gap_fraction > threshold or nt_count < min_nucleotide) and block_start is None:
            block_start = col_idx
        # If the column does not meet the gap requirement and is within a gap block, close the block
        elif gap_fraction <= threshold and nt_count >= min_nucleotide and block_start is not None:
            gap_blocks.append((block_start, col_idx))
            block_start = None

    # If the block extends to the end of the MSA, close the block
    if block_start is not None:
        gap_blocks.append((block_start, MSA_mafft.get_alignment_length()))

    # Filter blocks by conservation score
    delete_blocks = []

    if gap_blocks:
        for block in gap_blocks:
            col_before = MSA_mafft[:, block[0] - 1] if block[0] > 0 else None
            col_after = MSA_mafft[:, block[1]] if block[1] < MSA_mafft.get_alignment_length() else None
            if col_before is not None and calc_conservation(col_before) > conservation_threshold and \
                    col_after is not None and calc_conservation(col_after) > conservation_threshold:
                delete_blocks.append(block)

    # If no gap block was detected, return the original MSA
    else:
        # No blocks to delete, return the original MSA. Write into file with converted file name
        with open(fasta_out_flank_mafft_gap_filter_file, 'w') as f:
            AlignIO.write(MSA_mafft, f, 'fasta')
        return fasta_out_flank_mafft_gap_filter_file

    if delete_blocks:

        # Identify the sequence ranges to keep, which are in between the blocks to delete
        keep_ranges = [(0, delete_blocks[0][0])]
        for i in range(len(delete_blocks) - 1):
            keep_ranges.append((delete_blocks[i][1], delete_blocks[i + 1][0]))
        keep_ranges.append((delete_blocks[-1][1], MSA_mafft.get_alignment_length()))

        # Concatenate the columns in each keep_range object to create the final filtered MSA
        MSA_mafft_filtered = MSA_mafft[:, keep_ranges[0][0]:keep_ranges[0][1]]
        for keep_range in keep_ranges[1:]:
            MSA_mafft_filtered += MSA_mafft[:, keep_range[0]:keep_range[1]]

        # Write the filtered MSA to the output file
        with open(fasta_out_flank_mafft_gap_filter_file, 'w') as f:
            AlignIO.write(MSA_mafft_filtered, f, 'fasta')

        return fasta_out_flank_mafft_gap_filter_file
    else:
        # No blocks to delete, return the original MSA. Write into file with converted file name
        with open(fasta_out_flank_mafft_gap_filter_file, 'w') as f:
            AlignIO.write(MSA_mafft, f, 'fasta')
        return fasta_out_flank_mafft_gap_filter_file


def remove_gaps_with_similarity_check(input_file, output_dir, gap_threshold=0.8, simi_check_gap_thre=0.4,
                                      similarity_threshold=0.7, min_nucleotide=5, return_map=False):
    
    # This function replaces remove_gaps_block
    """
    Remove gaps flanked by conserved sequences in MSA. 
    Gaps are removed if: gap proportion is greater than 80%, less than 5 nucleotides are found in the column, 
    or if the gap proportion is between 40-80% and the nucleotide proportion is below 70% excluding gaps. 
    :param input_file: str, input file path
    :param output_dir: str, output file directory
    :param gap_threshold: int (0-1), columns with gap proportion greater than <gap_threshold> will be removed. Default: 0.8
    :param simi_check_gap_thre: int (0-1), if gap proportion greater or equal to <simi_check_gap_thre>, nucleotide
    proportions will be calculated for this column (excluding gaps). Default: 0.4
    :param similarity_threshold: int (0-1), if the gap proportion is between <simi_check_gap_thre> and <gap_threshold>
    and the proportion of the most abundant nucleotide is less than <similarity_threshold>, remove this column. Default: 0.7
    :param min_nucleotide: int (>0), columns with less than <min_nucleotide> nucleotides will be removed. Default: 5
    :return: absolute path to gap-purged MSA FASTA file
    """
    keep_list = []
    fasta_out_flank_mafft_gap_filter_file = os.path.join(output_dir,
                                                         f"{os.path.basename(input_file)}_gs.fa")
    MSA_mafft = AlignIO.read(input_file, "fasta")

    column_mapping = {}  # Stores the mapping of column indices of filtered MSA to original MSA

    #if len(MSA_mafft) < 5:
        #raise ValueError("Number of sequences is less than 5. Cannot remove gaps.")

    for col_idx in range(MSA_mafft.get_alignment_length()):
        col = MSA_mafft[:, col_idx]
        gap_count = col.count('-')
        gap_fraction = gap_count / len(col)

        # Check if the number of nucleotides in a column is less than 5
        nt_count = len(col) - gap_count
        if nt_count < min_nucleotide:
            continue

        # If the gap fraction is less than the gap_threshold, mark the column for further analysis
        if gap_fraction <= gap_threshold:
            # If the gap fraction is between 0.4 and 0.8, calculate proportions of remaining nucleotides
            if simi_check_gap_thre > gap_fraction:
                keep_list.append(col_idx)
                column_mapping[len(keep_list) - 1] = col_idx
            elif simi_check_gap_thre <= gap_fraction:
                nt_fraction = calc_conservation(col)

                # If the proportion of the most abundant nucleotide is less than similarity_threshold, skip this column
                if nt_fraction < similarity_threshold:
                    continue
                else:
                    # If the column passes all checks, add it to the keep_list
                    keep_list.append(col_idx)
                    # Store the mapping of original MSA index to filtered MSA index
                    column_mapping[len(keep_list) - 1] = col_idx

    # The index in Python does not include the second value. As a consequence, end_posit from the DefineBoundary()
    # function contains one more index than the final [len(keep_list) -1]. 
    # Adding one more value to correct index
    if len(keep_list) < 50:
        if return_map:
            return False, False
        else:
            return False
    column_mapping[len(keep_list)] = column_mapping[len(keep_list) - 1] + 1

    # Keep the columns
    MSA_mafft_filtered = MSA_mafft[:, keep_list[0]:keep_list[0] + 1]
    for i in keep_list[1:]:
        MSA_mafft_filtered += MSA_mafft[:, i:i + 1]

    # Write the filtered MSA to the output file
    with open(fasta_out_flank_mafft_gap_filter_file, 'w') as f:
        AlignIO.write(MSA_mafft_filtered, f, 'fasta')

    if return_map:
        return fasta_out_flank_mafft_gap_filter_file, column_mapping
    else:
        return fasta_out_flank_mafft_gap_filter_file


def remove_gaps_block_with_similarity_check(input_file, output_dir, gap_threshold=0.8,
                                            simi_check_gap_thre=0.4, similarity_threshold=0.7,
                                            conservation_threshold=0.6, min_nucleotide=5):
    """
    Remove gaps that are only from conserved regions. Remove columns when nucleotide number is less than 5.
    When gap percentage is equal or bigger than "simi_check_gap_thre". It will calculate most abundant nucleotide
    proportion for this column. If it is less than "similarity_threshold" this column will be removed
    :param input_file: str, input file path
    :param output_dir: str, output file directory
    :param gap_threshold: int (0-1), columns with gap proportion greater than <gap_threshold> will be removed. Default: 0.8
    :param simi_check_gap_thre: int (0-1), if gap proportion greater or equal to <simi_check_gap_thre>, nucleotide
    proportions will be calculated for this column (excluding gaps). Default: 0.4
    :param similarity_threshold: int (0-1), if the gap proportion is between <simi_check_gap_thre> and <gap_threshold>
    and the proportion of the most abundant nucleotide is less than <similarity_threshold>, remove this column. Default: 0.7
    :param conservation_threshold: int (0-1), after connected gaps have been classified as a gap block, the nucleotide
    conservation of two flanking non-gap columns at both ends of the gap block will be calculated. If the proportion of
    the most abundant nucleotide is greater than <conservation_threshold>, this gap block will be removed from the MSA. 
    Default: 0.6
    :param min_nucleotide: int (>0), columns with less than <min_nucleotide> nucleotides will be removed. Default: 5
    :return: absolute path to gap-purged MSA FASTA file
    """
    # Identify blocks of gaps
    gap_blocks = []

    # Output file name
    fasta_out_flank_mafft_gap_filter_file = os.path.join(output_dir, f"{os.path.basename(input_file)}_gbs.fa")

    # Load the MSA file
    MSA_mafft = AlignIO.read(input_file, "fasta")

    # Raise an error if the number of sequences is less than 5
    if len(MSA_mafft) < 5:
        raise ValueError("Number of sequences in MSA is less than 5. Cannot remove gaps.")

    # Define the starting point of a block
    block_start = None

    # Go through each column in the alignment
    for col_idx in range(MSA_mafft.get_alignment_length()):

        # Define metrics for the column
        col = MSA_mafft[:, col_idx]
        gap_count = col.count('-')
        gap_fraction = gap_count / len(col)
        nt_fraction = calc_conservation(col)
        nt_count = len(col) - gap_count

        # Check if the column meets the criteria to open a new block
        if ((simi_check_gap_thre <= gap_fraction <= gap_threshold and nt_fraction <= similarity_threshold) or
            gap_fraction > gap_threshold or nt_count < min_nucleotide) and block_start is None:
            block_start = col_idx

        # Check if the column meets the criteria to close a block
        elif simi_check_gap_thre <= gap_fraction <= gap_threshold and nt_fraction > \
                similarity_threshold and nt_count >= min_nucleotide and block_start is not None:
            gap_blocks.append((block_start, col_idx))
            block_start = None

        # Check if the block should be closed due to gap proportion below simi_check_gap_thre in a column
        elif gap_fraction <= simi_check_gap_thre and nt_count >= min_nucleotide and block_start is not None:
            gap_blocks.append((block_start, col_idx))
            block_start = None

    # If a block was opened but not closed, close it at the last column of the MSA
    if block_start is not None:
        gap_blocks.append((block_start, MSA_mafft.get_alignment_length()))

    # Identify the blocks that should be deleted based on the conservation of the surrounding columns
    delete_blocks = []

    # Check if gap_blocks is empty
    if gap_blocks:
        for block in gap_blocks:
            # Calculate the nucleotide divergence of flanking columns of the gap block
            col_before = MSA_mafft[:, block[0] - 1] if block[0] > 0 else None
            col_before_2 = MSA_mafft[:, block[0] - 2] if block[0] > 1 else None
            col_after = MSA_mafft[:, block[1]] if block[1] < MSA_mafft.get_alignment_length() else None
            col_after_2 = MSA_mafft[:, block[1] + 1] if block[1] < MSA_mafft.get_alignment_length() - 1 else None
            # If the flanking columns contain conserved nucleotides, mark the block for deletion
            if col_before is not None and (calc_conservation(col_before) > conservation_threshold or
                                           calc_conservation(
                                               col_before_2) > conservation_threshold) and col_after is not None and \
                    (calc_conservation(col_after) > conservation_threshold or
                     calc_conservation(col_after_2) > conservation_threshold):
                delete_blocks.append(block)
    else:
        # No blocks to delete, return the original MSA. Write into file with converted file name
        with open(fasta_out_flank_mafft_gap_filter_file, 'w') as f:
            AlignIO.write(MSA_mafft, f, 'fasta')
        return fasta_out_flank_mafft_gap_filter_file

    # Check if "delete_blocks" is empty
    if delete_blocks:

        # Define the ranges of the MSA that should be kept
        keep_ranges = [(0, delete_blocks[0][0])]
        for i in range(len(delete_blocks) - 1):
            keep_ranges.append((delete_blocks[i][1], delete_blocks[i + 1][0]))
        keep_ranges.append((delete_blocks[-1][1], MSA_mafft.get_alignment_length()))

        # Create a new MSA by concatenating the kept ranges
        MSA_mafft_filtered = MSA_mafft[:, keep_ranges[0][0]:keep_ranges[0][1]]
        for keep_range in keep_ranges[1:]:
            MSA_mafft_filtered += MSA_mafft[:, keep_range[0]:keep_range[1]]

        # Write the filtered MSA to a file
        with open(fasta_out_flank_mafft_gap_filter_file, 'w') as f:
            AlignIO.write(MSA_mafft_filtered, f, 'fasta')

        # Return the path to the output file
        return fasta_out_flank_mafft_gap_filter_file

    else:
        # No blocks to delete, return the original MSA. Write file again to convert file name
        with open(fasta_out_flank_mafft_gap_filter_file, 'w') as f:
            AlignIO.write(MSA_mafft, f, 'fasta')
        return fasta_out_flank_mafft_gap_filter_file


def select_gaps_block_with_similarity_check(input_file,
                                            simi_check_gap_thre=0.4, similarity_threshold=0.8,
                                            conservation_threshold=0.6, min_nucleotide=10):
    """
    :return: A list containing column numbers used for clustering
    """
    # Identify blocks of gaps
    gap_blocks = []

    # Load the MSA file
    MSA_mafft = AlignIO.read(input_file, "fasta")

    # Raise an error if the number of sequences in the MSA is less than 5
    if len(MSA_mafft) < 5:
        raise ValueError("Number of sequences in MSA is less than 5. Cannot remove gaps.")

    # Define the starting point of a block
    block_start = None

    # Variables to keep track of the gap count in the previous column and the condition met count
    prev_gap_count = None
    condition_met_count = 0

    # Go through each column in the alignment. Start from 0
    for col_idx in range(MSA_mafft.get_alignment_length()):
        # Define metrics for the column
        col = MSA_mafft[:, col_idx]
        gap_count = col.count('-')
        gap_fraction = gap_count / len(col)
        nt_fraction = calc_conservation(col)
        nt_count = len(col) - gap_count

        # Check if the column meets the criteria to open a new block
        if simi_check_gap_thre <= gap_fraction and nt_fraction >= similarity_threshold and \
                nt_count >= min_nucleotide and block_start is None:
            block_start = col_idx
            condition_met_count = 0  # Reset the count for a new block

        # Additional check for gap variance
        elif prev_gap_count is not None and (
                gap_count < 0.8 * prev_gap_count or gap_count > 1.2 * prev_gap_count) and block_start is not None:
            gap_blocks.append((block_start, col_idx, condition_met_count))
            block_start = None

        # Close the gap block if nucleotide similarity is smaller than the threshold
        elif (
                simi_check_gap_thre <= gap_fraction and nt_fraction < similarity_threshold and block_start is not None) or \
                (nt_count < min_nucleotide and block_start is not None):
            # Check the gap count of the next column if it exists
            next_gap_count = None
            if col_idx + 1 < MSA_mafft.get_alignment_length():
                next_col = MSA_mafft[:, col_idx + 1]
                next_gap_count = next_col.count('-')

            # If gap count is not similar to both previous and next columns, close the block
            if prev_gap_count is not None and next_gap_count is not None:
                if not (0.8 * prev_gap_count <= gap_count <= 1.2 * prev_gap_count and
                        0.8 * next_gap_count <= gap_count <= 1.2 * next_gap_count):
                    condition_met_count += 1

        # Stop the gap block if the gap number is below threshold simi_check_gap_thre
        elif gap_fraction < simi_check_gap_thre and block_start is not None:
            gap_blocks.append((block_start, col_idx, condition_met_count))
            block_start = None

        # Update the gap count of the previous column
        prev_gap_count = gap_count

    # Process the identified gap blocks
    keep_blocks = []

    if gap_blocks:
        for block in gap_blocks:
            start, end, block_condition_count = block  # Now also unpacking the block_condition_count

            # Calculate the maximum allowable condition_met_count based on the gap block length
            block_length = end - start
            if block_length < 10:
                max_condition_count = 1
            else:
                max_condition_count = min(0.1 * block_length, 50)

            # Check if the condition_met_count exceeds the maximum allowable value for the gap block
            if block_condition_count > max_condition_count:
                continue

            if start and end and block_length > 10:
                if start > 0:
                    col_before = [MSA_mafft[i, start - 1] for i in range(len(MSA_mafft[:, start])) if
                                  MSA_mafft[i, start] == '-']
                    gap_fraction_before = col_before.count('-') / len(col_before)
                else:
                    col_before = None
                    gap_fraction_before = 0

                if end < MSA_mafft.get_alignment_length():
                    col_after = [MSA_mafft[i, end] for i in range(len(MSA_mafft[:, end - 1])) if
                                 MSA_mafft[i, end - 1] == '-']
                    gap_fraction_after = col_after.count('-') / len(col_after)
                else:
                    col_after = None
                    gap_fraction_after = 0

                if gap_fraction_before <= 0.3 and gap_fraction_after <= 0.3 and \
                        col_before is not None and calc_conservation(col_before) > conservation_threshold and \
                        col_after is not None and calc_conservation(col_after) > conservation_threshold:
                    keep_blocks.append((start, end))

    if keep_blocks:
        flat_keep_blocks = [item for start, end in keep_blocks for item in range(start, end)]
        return flat_keep_blocks
    else:
        return False


# Select MSA columns according to the provided start and end position
def select_star_to_end(input_file, output_dir, start, end):
    alignment = AlignIO.read(input_file, "fasta")
    MSA_len = alignment.get_alignment_length()
    # Ensure the start position is within the sequence range
    start = max(start, 0)
    # Ensure the end position is within the sequence range
    end = min(end, MSA_len)
    # Select the window columns from the alignment
    select_alignment = alignment[:, start:end]
    # Create MultipleSeqAlignment object with the select alignment
    select_alignment_object = MultipleSeqAlignment(select_alignment)
    # Define the output file
    output_file = os.path.join(output_dir, f"{os.path.basename(input_file)}_se.fa")
    AlignIO.write(select_alignment_object, output_file, "fasta")

    return output_file


def select_start_end_and_join(input_file, output_dir, start, end, window_size=50):
    """
    Select start and end columns of MSA
    :param input_file: str, absolute input file path
    :param output_dir: str, output directory
    :param start: int, start point, columns to the right of the start point will be selected
    :param end: int, end point, columns to the left of the end point will be selected
    :param window_size: int, number of columns beginning from <start> and <end> to be selected. Default: 50
    :return: Selected MSA object (no file path)
    """
    alignment = AlignIO.read(input_file, "fasta")
    sequence_length = end - start

    new_alignment = []  # Define a list to store new alignment file

    for record in alignment:
        sequence = record.seq
        selected_sequence = str(sequence[start:start + window_size]) + '----------' \
                            + str(sequence[end - window_size:end])
        new_record = SeqRecord(Seq(selected_sequence), id=record.id, description="")
        new_alignment.append(new_record)

    new_alignment = MultipleSeqAlignment(new_alignment)
    output_file = os.path.join(output_dir, f"{os.path.basename(input_file)}_plot.fa")
    AlignIO.write(new_alignment, output_file, "fasta")

    return new_alignment, sequence_length


def select_window_columns(input_file, output_dir, start_point, direction, window_size=50):
    """
    Select column block for MSA
    :param input_file: str, absolute input file path
    :param output_dir: str, output directory
    :param start_point: int, start point to select column block
    :param direction: str "right" or "left", direction to select column block
    :param window_size: int, number of columns beginning from <start_point> to be selected. Default: 50
    :return: Select MSA object (not file path)
    """
    # Read the MSA file using Biopython's AlignIO
    alignment = AlignIO.read(input_file, "fasta")

    # Get the total number of columns in the alignment
    total_columns = alignment.get_alignment_length()

    # Determine the start and end column indices based on the direction and start point
    if direction.lower() == 'left':
        start_col = max(0, start_point - window_size)
        end_col = start_point
    elif direction.lower() == 'right':
        start_col = start_point
        end_col = min(start_point + window_size, total_columns)  # get as many columns as possible
    else:
        raise ValueError("Invalid direction. Please enter 'left' or 'right'.")

    # Select the window columns from the alignment
    selected_alignment = alignment[:, start_col:end_col]

    # Convert selected_alignment to a MultipleSeqAlignment object
    selected_alignment = MultipleSeqAlignment(selected_alignment)

    output_file = os.path.join(output_dir, f"{os.path.basename(input_file)}_{direction}_{window_size}.fa")
    AlignIO.write(selected_alignment, output_file, "fasta")
    # Return the new alignment
    return selected_alignment


"""
The asterisk '*' before alignments in the definition of the function allows the function to accept any number 
of positional arguments. These arguments will be gathered into a tuple called alignments. This is particularly 
useful if you do not know beforehand how many arguments will be passed to the function.
"""


def concatenate_alignments(*alignments, input_file_name, output_dir):
    """
    Concatenate MSA files containing the same number of sequences
    :param alignments: str, alignment objects, will be concatenated in given order
    :param input_file_name: str, used for generating output file name
    :param output_dir: str, the absolute output directory
    :return: the absolute output file path, the start and end points of the concatenated MSA
    """
    # Check if at least two alignments were provided
    if len(alignments) < 2:
        raise ValueError("At least two alignments must be provided.")

    # Check if all alignments have the same number of sequences
    num_sequences = len(alignments[0])  # Give sequence number of MSA
    if not all(len(alignment) == num_sequences for alignment in alignments):
        raise ValueError("All alignments must have the same number of sequences.")

    # Concatenate the alignments
    concatenated_alignment = alignments[0]
    for alignment in alignments[1:]:
        concatenated_alignment += alignment

    # alignments[0].get_alignment_length() will return 50, but python starts with 0
    # The last position of alignments[0] is alignments[0].get_alignment_length() - 1
    concat_start = alignments[0].get_alignment_length()
    concat_end = concatenated_alignment.get_alignment_length() - alignments[-1].get_alignment_length() - 1

    # Write to a file
    output_file = os.path.join(output_dir, f"{os.path.basename(input_file_name)}_me.fa")
    AlignIO.write(concatenated_alignment, output_file, "fasta")

    # Return the concatenated alignment
    return output_file, concat_start, concat_end


def cd_hit_est(input_file, output_file, identity_thr=0.8, aL=0.9, aS=0.9, s=0.9, thread=10):
    """
    -l int, length of throw_away_sequences. Default: 10
    -d int, length of description in .clstr file; if set to 0, it takes the FASTA defline and stops at the
    first space. Default: 20
    -s float, length difference cutoff; if set to 0.9, the shorter sequences need to be at least 90% in length of the
    representative sequence of the cluster. Default: 0.0
    -aL	float, alignment coverage for the longer sequence; if set to 0.9, the alignment must cover 90% of the sequence.
    Default: 0.0
    -aS	float, alignment coverage for the shorter sequence; if set to 0.9, the alignment must cover 90% of the sequence.
    Default: 0.0
    Note: -s only considers length, but -aL and -aS considers alignment
    Note: -sc sort cluster according sequence numbers in the cluster
    """
    command = [
        "cd-hit-est",
        "-i", input_file,
        "-o", output_file,
        "-c", str(identity_thr),
        "-aL", str(aL),
        "-aS", str(aS),
        "-M", "0",
        "-T", str(thread),
        "-l", "30",
        "-d", "0",
        "-s", str(s),
        "-sc", "1"
    ]

    try:
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    except FileNotFoundError:
        prcyan("'\ncd-hit-est' command not found. Please ensure 'cd-hit-est' is installed correctly.\n")
        raise Exception

    except subprocess.CalledProcessError as e:
        prcyan(f"\ncd-hit-est failed for {os.path.basename(input_file)} with error code {e.returncode}")
        prcyan(e.stderr)
        raise Exception


def parse_cd_hit_est_result(input_file):
    clusters = {}  # The key is cluster name, the values are sequence names in this cluster (list)
    detailed_clusters = {}  # Key: cluster name, Value: list of tuples (sequence length, sequence name, percentage)
    current_cluster = []
    current_detailed_cluster = []
    with open(input_file, "r") as f:
        for line in f:
            # cd-hit-est introduces empty spaces in cluster headers like ">Cluster 53", which can cause
            # errors in downstream analysis. The following code fixes the issue.
            if line.startswith(">Cluster"):
                if current_cluster:  # If the current_cluster is not empty
                    # clusters is a dictionary, the key is cluster_name
                    clusters[cluster_name] = current_cluster
                    detailed_clusters[cluster_name] = current_detailed_cluster
                # Remove the empty space and ">" in the cluster name
                cluster_name = line.strip().replace(" ", "").replace(">", "")
                current_cluster = []
                current_detailed_cluster = []
            else:
                """
                This is an example of part of cd-hit-est output 
                >Cluster 198
                0	5965nt, >scaffold_23_3352139..3358595#scaffold_23:3352139..3358595... at +/99.55%
                1	5966nt, >scaffold_23_1868463..1880051_01#scaffold_23:1868463..1880051... at +/99.53%
                2	5972nt, >scaffold_25_1053723..1064724#scaffold_25:1053723..1064724... at +/99.45%
                3	5962nt, >scaffold_29_441223..447627#scaffold_29:441223..447627... at +/99.53%
                4	5977nt, >scaffold_32_604726..610694#scaffold_32:604726..610694... at +/99.20%
                5	5969nt, >scaffold_8_5853353..5859876#scaffold_8:5853353..5859876... at +/99.51%
                6	5969nt, >scaffold_8_5853353..5859876_01#scaffold_8:5853353..5859876... at +/99.53%
                7	5970nt, >TE_00000565_INT#TE_00000565_INT... at +/99.51%
                8	5982nt, >TE_00000545_LTR#TE_00000545_LTR... *
                """
                seq_name = line.split(">")[1].split("...")[0].split("#")[0].strip()
                try:
                    seq_length = line.split('nt,')[0].split('\t')[1].strip()
                except Exception:
                    seq_length = None
                if '... at' in line:
                    seq_per_and_direction = line.split('... at')[1].strip()
                    seq_per = seq_per_and_direction.split('/')[1].strip()
                    seq_direction = seq_per_and_direction.split('/')[0].strip()

                    if seq_direction != '+' and seq_direction != '-':
                        seq_direction = None

                elif '... *' in line:
                    seq_per = "standard"
                    seq_direction = '+'
                else:
                    seq_per = None
                    seq_direction = None
                current_cluster.append(seq_name)
                current_detailed_cluster.append((seq_length, seq_name, seq_per, seq_direction))
        if current_cluster:
            clusters[cluster_name] = current_cluster
            detailed_clusters[cluster_name] = current_detailed_cluster
    return clusters, detailed_clusters


def repeatmasker(genome_file, library_file, output_dir, thread=1, classify=False):
    """
    Run RepeatMasker with the provided parameters.
    """

    # Construct the RepeatMasker command
    if classify:
        command = ["RepeatMasker",
                   genome_file,
                   "-lib", library_file,
                   "-s",  # Slow search; 0-5% more sensitive, 2-3 times slower than default
                   "-dir", output_dir,
                   "-pa", str(thread)
                   ]
    else:
        command = ["RepeatMasker",
                   genome_file,
                   "-lib", library_file,
                   "-pa", str(thread),
                   "-dir", output_dir,
                   "-s",  # Slow search; 0-5% more sensitive, 2-3 times slower than default
                   "-gff",  # Creates an additional Gene Feature Finding format output
                   "-xm",  # Creates an additional output file in cross_match format (for parsing)
                   "-a",  # Writes alignments in .align output file
                   ]
    try:
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True

    except FileNotFoundError:
        prcyan("'RepeatMasker' command not found. Please ensure 'RepeatMasker' is installed correctly.")
        raise Exception

    except subprocess.CalledProcessError as e:
        if classify:
            prcyan(f"\nRepeatMasker failed during final classification step with error code {e.returncode}")
            prcyan(e.stderr)
            prgre("This will not affect the final result. Only the classification of TE may not be correct.")
            raise Exception
        else:
            prcyan(f"\nRepeatMasker failed during final whole-genome TE annotation with error code {e.returncode}")
            prcyan(e.stderr)
            prgre("This does not affect the final TE consensus library. You can perform the final genome-wide TE"
                  " annotation by yourself with RepeatMasker.")
            raise Exception


def repeatmasker_output_classify(repeatmasker_out, progress_file, min_iden=70, min_len=80, min_cov=0.8):
    # Read RepeatMasker output file (.out) into a DataFrame
    # The regex '\s+' matches one or more whitespace characters
    # error_bad_lines=False to skip errors
    try:
        df = pd.read_csv(repeatmasker_out, sep='\s+', header=None, skiprows=3, usecols=range(15))
    except pandas.errors.EmptyDataError:
        return False
    except pd.errors.ParserError:
        df = pd.read_csv(repeatmasker_out, sep='\s+', header=None, skiprows=3, error_bad_lines=False, usecols=range(15))

    # Rename columns for easier referencing
    df.columns = [
        'score', 'perc_div', 'perc_del', 'perc_ins', 'query_name',
        'query_start', 'query_end', 'query_left', 'strand',
        'repeat_name', 'repeat_class', 'repeat_start',
        'repeat_end', 'repeat_left', 'ID'
    ]

    # Filter rows based on query identity
    df = df[df['perc_div'] <= (100 - min_iden)]

    # Calculate coverage length and add to a new column cov_len
    df['cov_len'] = abs(df['query_start'] - df['query_end'])

    # Select dataframe columns
    df_filter = df[["query_name", "repeat_name", "repeat_class", "cov_len"]]

    # Group by columns and calculate sum of 'cov_len'
    grouped_df = df_filter.groupby(['query_name', 'repeat_name', 'repeat_class'])['cov_len'].sum().reset_index()
    grouped_df_filter = grouped_df[grouped_df['cov_len'] >= min_len]

    # Group by 'repeat_name' and get the index of the row with the maximum 'cov_len'
    idx = grouped_df_filter.groupby('query_name')['cov_len'].idxmax()

    # Use these indices to filter the DataFrame
    max_cov_len_df = grouped_df_filter.loc[idx]

    # Convert max_cov_len_df to a dictionary for faster look-up
    max_cov_dict = max_cov_len_df.set_index('query_name')['cov_len'].to_dict()

    # Read progress file
    progress_df = pd.read_csv(progress_file)

    # Dictionary to store consensus_name and its reclassified_type
    reclassified_dict = {}

    # Iterate over each row in progress_df
    for index, row in progress_df.iterrows():
        consensus_name = row['consensus_name']
        cons_length = row['cons_length']
        cons_type = row['reclassified_type']

        # Check if consensus_name exists in max_cov_dict and compute the ratio
        if consensus_name in max_cov_dict:
            cov_len = max_cov_dict[consensus_name]
            ratio = int(cov_len) / int(cons_length)

            # Check if ratio meets the threshold
            if ratio >= min_cov and "Unknown" in cons_type:
                # Find the corresponding repeat_class for this consensus_name
                repeat_class = max_cov_len_df[max_cov_len_df['query_name'] == consensus_name]['repeat_class'].iloc[0]

                # Modify the reclassified_type column in progress_df
                progress_df.at[index, 'reclassified_type'] = repeat_class

                # Update the dictionary with the new reclassified_type
                reclassified_dict[consensus_name] = repeat_class

    # Save the modified progress_df back to the original file
    progress_df.to_csv(progress_file, index=False, na_rep='NaN')

    return reclassified_dict


def rename_cons_file(consensus_file, reclassified_dict):
    # Define a temporary file for the updated content
    temp_file = consensus_file + ".tmp"

    with open(consensus_file, 'r') as infile, open(temp_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # Extract the sequence name (up to the '#' character)
                seq_name = line.split('#')[0][1:].strip()  # Remove '>' and split at '#', then keep the first part

                # Look up the sequence name in the reclassified_dict and modify the header if found
                if seq_name in reclassified_dict:
                    new_type = reclassified_dict[seq_name]
                    line = f">{seq_name}#{new_type}\n"

            outfile.write(line)

    # Use shutil.move to replace the original file with the temporary file
    shutil.move(temp_file, consensus_file)


def rename_files_based_on_dict(directory, reclassified_dict, seq_name=False):
    # List all files in the directory
    files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]

    # For each file, check and rename if necessary
    for filename in files:
        # Extract the part before '#'
        consensus_name = filename.split('#')[0]
        for key, value in reclassified_dict.items():
            if consensus_name == key:
                # Extract the part after # and before the first .
                te_type = filename.split('#')[1].split('.')[0].replace('__', '/')

                # If file name does not match the value in the dictionary, rename the file
                if te_type != value:
                    new_te_type = value.replace('/', '__')
                    new_filename = filename.replace(te_type.replace('/', '__'), new_te_type)
                    old_file_path = os.path.join(directory, filename)
                    new_file_path = os.path.join(directory, new_filename)

                    # If seq_name is 'True', rename the sequence header in the FASTA file
                    if seq_name:
                        record = SeqIO.read(old_file_path, 'fasta')
                        record.id = f"{seq_name}#{value}"
                        record.description = ""  # Clear the description to avoid duplication
                        SeqIO.write(record, new_file_path, 'fasta')
                        # Delete the old file
                        os.remove(old_file_path)
                    else:
                        # Use shutil.move to rename the file
                        shutil.move(old_file_path, new_file_path)


def remove_files_with_start_pattern(input_dir, start_pattern=None):
    # Remove files and folders whose name starts with given pattern
    # When start pattern is missing, remove everything in this folder
    for filename in os.listdir(input_dir):
        file_path = os.path.join(input_dir, filename)

        if start_pattern is not None and not filename.startswith(start_pattern):
            continue

        if os.path.isfile(file_path):
            os.remove(file_path)
        elif os.path.isdir(file_path):
            try:
                shutil.rmtree(file_path)
            except Exception as e:
                # File deletion does not affect the final consensus sequence. Skip if an error occurs
                pass


def copy_files_with_start_pattern(input_dir, start_pattern, output_dir, seq_len=None, seq_pre=None, evaluation=None):
    # List all files in the input directory
    for file in os.listdir(input_dir):
        # Check if the file name starts with the given pattern
        if file.startswith(start_pattern):
            # Construct full file path
            full_file_path = os.path.join(input_dir, file)

            # Split the file name into the name and extension
            name, extension = os.path.splitext(file)

            # Define destination file
            if extension == ".pdf":
                if evaluation is not None:
                    name = f"{name}_{evaluation}"
                if seq_len is not None:
                    name = f"{name}_{seq_len}_bp"
                if seq_pre is not None:
                    name = f"{name}_{seq_pre}"
                output_file = os.path.join(output_dir, f"{name}.pdf")
            else:
                output_file = os.path.join(output_dir, file)
            # Copy file to the output directory
            shutil.copy(full_file_path, output_file)


# Define a function to handle sequence skipping and removal of files
def handle_sequence_low_copy(seq_obj, progress_file, debug, MSA_dir, classification_dir,
                             found_match=None, blast_full_length_n=None, te_aid_plot=None, orf_plot=None,
                             low_copy_dir=None):
    seq_name = seq_obj.get_seq_name()
    te_type = seq_obj.get_old_TE_type()
    te_type_modified = te_type.replace("/", "__")
    try:
        if found_match is not None and blast_full_length_n is not None:
            seq_obj.set_old_terminal_repeat(found_match)
            seq_obj.set_old_blast_full_n(blast_full_length_n)
            seq_obj.update_status("processed", progress_file)

            if (te_aid_plot is not None or orf_plot is not None) and low_copy_dir is not None:
                # Merge TE Aid and ORF plots
                merge_pdfs(low_copy_dir, f"{seq_name}#{te_type_modified}", te_aid_plot, orf_plot)
        if not debug:
            remove_files_with_start_pattern(MSA_dir, f"{seq_name}.fasta")
            remove_files_with_start_pattern(classification_dir, f"{seq_name}.fasta")
    except Exception as e:
        click.echo(f"\nAn error occurred while handling low-copy sequence {seq_name}:\n {e}\n")


def handle_sequence_skipped(seq_obj, progress_file, debug, MSA_dir, classification_dir, plot_skip=True,
                            te_aid_plot=None, orf_plot=None, skip_proof_dir=None):
    seq_name = seq_obj.get_seq_name()
    te_type = seq_obj.get_old_TE_type()
    te_type_modified = te_type.replace("/", "__")
    input_fasta = seq_obj.get_input_fasta()

    try:
        seq_obj.update_status("skipped", progress_file)
        if plot_skip and (te_aid_plot is not None or orf_plot is not None) and skip_proof_dir is not None:
            merge_pdfs(skip_proof_dir, f"{seq_name}#{te_type_modified}", te_aid_plot, orf_plot)

            # Copy skipped input sequence into skip_proof_dir
            skip_fasta_file = os.path.join(skip_proof_dir, f"{seq_name}#{te_type_modified}.fa")
            shutil.copy(input_fasta, skip_fasta_file)
        if not debug:
            remove_files_with_start_pattern(MSA_dir, f"{seq_name}.fasta")
            remove_files_with_start_pattern(classification_dir, f"{seq_name}.fasta")
    except Exception as e:
        click.echo(f"\nAn error occurred while handling skipped sequence {seq_name}:\n {e}\n")


def update_cons_file(updated_type, unknown_concensus_file, consensus_file):
    if os.path.exists(unknown_concensus_file):
        with open(unknown_concensus_file, 'r') as fasta_file:
            for record in SeqIO.parse(fasta_file, 'fasta'):
                header = record.id
                sequence = str(record.seq)
                if header in updated_type:
                    te_type = updated_type[header]
                else:
                    te_type = "Unknown"
                with open(consensus_file, 'a') as f:
                    f.write(">" + header + "#" + te_type + "\n" + sequence + "\n")


# if the sequence object seq_obj contains a low-copy TE, append to consensus_file and final_unknown_con_file used for classification
def update_low_copy_cons_file(seq_obj, consensus_file, final_unknown_con_file, final_classified_con_file, proof_dir,
                              te_aid_pdf):
    seq_name = seq_obj.get_seq_name()
    te_type = seq_obj.get_old_TE_type()
    te_type_modified = te_type.replace("/", "__")
    input_fasta = seq_obj.get_input_fasta()

    record = SeqIO.read(input_fasta, "fasta")
    sequence = str(record.seq)

    if "Unknown" in te_type:
        with open(final_unknown_con_file, "a") as f:  # 'a' mode for appending
            f.write(">" + seq_name + "\n" + sequence + "\n")
    else:
        with open(final_classified_con_file, "a") as f:
            f.write(">" + seq_name + "#" + te_type + "\n" + sequence + "\n")

        # Write all consensus sequences to final_cons_file.
    with open(consensus_file, "a") as f:
        f.write(">" + seq_name + "#" + te_type + "\n" + sequence + "\n")

    low_copy_single_fasta_file = os.path.join(proof_dir, f"{seq_name}#{te_type_modified}.fa")
    #low_copy_te_aid_pdf_file = os.path.join(proof_dir, f"{seq_name}#{te_type_modified}_TE_Aid.pdf")

    shutil.copy(input_fasta, low_copy_single_fasta_file)

    # Sometimes TE-Aid cannot be plotted properly due to low input sequence quality. Only move plot if it exists.
    #if os.path.exists(te_aid_pdf) and os.path.getsize(te_aid_pdf) > 0:
        #shutil.copy(te_aid_pdf, low_copy_te_aid_pdf_file)


# Classify single FASTA file
def classify_single(consensus_fasta):
    """
    Run RepeatClassifier with the provided parameters.
    """

    # Store the current working directory
    original_dir = os.getcwd()

    # Change the working directory to the directory of the consensus_fasta
    os.chdir(os.path.dirname(consensus_fasta))

    # Define RepeatClassifier command, the output file will be stored in the same directory as consensus_fasta
    command = ["RepeatClassifier", "-consensi", consensus_fasta]

    try:
        # Run RepeatClassifier using subprocess
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    except FileNotFoundError:
        prcyan("'RepeatClassifier' command not found. Please ensure 'RepeatModeler' is installed correctly.")
        return False

    except subprocess.CalledProcessError as e:
        prcyan(f"RepeatClassifier failed for {os.path.basename(consensus_fasta)} with error code {e.returncode}")
        prcyan(e.stderr)
        prgre("This only affects classification but not the consensus sequence. "
              "You can run 'RepeatClassifier -consensi <your_consensus_file>' manually.")
        # Change the working directory back to the original directory
        os.chdir(original_dir)
        return False

    # Change the working directory back to the original directory
    os.chdir(original_dir)

    classified_file = f'{consensus_fasta}.classified'

    # Get the first record of file with classified consensus sequences
    record = next(SeqIO.parse(classified_file, "fasta"))

    # seq_name = record.id.split("#")[0]
    seq_TE_type = record.id.split("#")[-1]

    return seq_TE_type


def check_terminal_repeat(input_file, output_dir, teaid_blast_out=None, TIR_adj=2000, LTR_adj=3000):
    """
    output_dir: used to store self-BLAST database
    """
    # Read input file and get sequence length
    record = SeqIO.read(input_file, "fasta")
    record_len = len(record.seq)
    found_match = False
    LTR_boundary = None
    TIR_boundary = None

    # If TE-Aid output was not provided, the file does not exist or is empty, do self-BLAST
    if teaid_blast_out is None or not file_exists_and_not_empty(teaid_blast_out):
        os.makedirs(output_dir, exist_ok=True)
        database_file = os.path.join(output_dir, "Tem_blast_database")
        makeblastdb_cmd = f"makeblastdb -in {input_file} -dbtype nucl -out {database_file}"

        # If an error was encountered, abort and return
        try:
            subprocess.run(makeblastdb_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError:
            return None, None, False

        # Check if database was built, otherwise return
        if not file_exists_and_not_empty(f"{database_file}.nhr") or not \
                file_exists_and_not_empty(f"{database_file}.nin") or not file_exists_and_not_empty(f"{database_file}.nsq"):
            return None, None, False

        blast_cmd = f"blastn -query {input_file} -db {database_file} " \
                    f"-outfmt \"6 qseqid qstart qend sstart send \" " \
                    f"-evalue 0.05"  # Set more stringent e-value for self-BLAST

        try:
            result = subprocess.run(blast_cmd, shell=True, check=True, text=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            blast_out = result.stdout

            # If self-BLAST result is empty, return
            if blast_out.strip() == "":
                return None, None, False
        except subprocess.CalledProcessError as e:
            prcyan(f"\nTerminal repeat detection failed for {os.path.basename(input_file)} with error code {e.returncode}")
            prcyan('\n' + e.stderr + '\n')
            raise Exception

        # Define BLAST output file
        blast_out_file = os.path.join(output_dir, "tem_blast_out.txt")

        with open(blast_out_file, 'w') as f:
            # Add a header to BLAST result
            f.write("qseqid\tqstart\tqend\tsstart\tsend\n")
            f.write(blast_out)
        df = pd.read_csv(blast_out_file, sep="\t", header=None, skiprows=1)
    else:
        blast_out_file = teaid_blast_out

        # TE-Aid self-BLAST output default format is separated by white space
        df = pd.read_csv(blast_out_file, sep='\s+', header=None, skiprows=1)

    # Return 'None' if the self-BLAST result is empty. This can happen if sequence contains too many ambiguous letters
    # like 'N' or 'X'.
    if df.empty:
        return None, None, False

    df_LTR = df[(df.iloc[:, 2] - df.iloc[:, 1] >= 150) &
            (df.iloc[:, 1] != df.iloc[:, 3]) &
            (df.iloc[:, 2] != df.iloc[:, 4]) &
            (df.iloc[:, 3] < df.iloc[:, 4]) &
            (df.iloc[:, 3] > df.iloc[:, 1])].copy()

    if not df_LTR.empty:
        df_LTR["5"] = df.iloc[:, 4] - df.iloc[:, 1]
        df_LTR.reset_index(drop=True, inplace=True)

        # Find the row with the largest difference
        LTR_largest = df_LTR.iloc[df_LTR["5"].idxmax()]

        # Check if the terminal repeat spans the majority of the query sequence. Because self-BLAST is done after extension,
        # the maximum redundant extension for left and right side are both 2000.
        if abs(LTR_largest[4] - LTR_largest[1]) >= (record_len - LTR_adj):
            # Because blast use index start from 1, modify the start position
            LTR_boundary = [LTR_largest[1] - 1, LTR_largest[4]]
            found_match = "LTR"
        else:
            LTR_boundary = None

        # Return directly if LTR was found
        return LTR_boundary, TIR_boundary, found_match

    df_TIR = df[(df.iloc[:, 2] - df.iloc[:, 1] >= 50) &
            (df.iloc[:, 1] != df.iloc[:, 3]) &
            (df.iloc[:, 2] != df.iloc[:, 4]) &
            (df.iloc[:, 3] > df.iloc[:, 4]) &
            (df.iloc[:, 4] > df.iloc[:, 1])].copy()
    
    if not df_TIR.empty:
        df_TIR["5"] = df.iloc[:, 3] - df.iloc[:, 1]
        df_TIR.reset_index(drop=True, inplace=True)

        # Same like LTR, check the terminal repeat-spanning region
        TIR_largest = df_TIR.iloc[df_TIR["5"].idxmax()]
        if abs(TIR_largest[3] - TIR_largest[1]) >= (record_len - TIR_adj):
            TIR_boundary = [TIR_largest[1] - 1, TIR_largest[3]]
            found_match = "TIR"
        else:
            TIR_boundary = None
    else:
        TIR_boundary = None

    return LTR_boundary, TIR_boundary, found_match


def filter_out_big_gap_seq(input, output=None, gap_threshold=1):
    if output is None:
        alignment = input
    else:
        # Read the alignment from the input file
        alignment = AlignIO.read(input, "fasta")
    alignment_len = alignment.get_alignment_length()

    # Filter based on gap fraction
    gap_alignment_filter_list = []
    for record in alignment:
        gap_count = record.seq.count("-")
        gap_fraction = gap_count / alignment_len

        if gap_fraction < gap_threshold:
            gap_alignment_filter_list.append(record)

    # Create a new MultipleSeqAlignment object with the filtered records
    filtered_alignment = MultipleSeqAlignment(gap_alignment_filter_list)

    if output is None:
        return filtered_alignment
    else:
        # Write the filtered alignment to the output file
        AlignIO.write(filtered_alignment, output, "fasta")


def file_exists_and_not_empty(file_path):
    # Check if the file exists
    if os.path.isfile(file_path):
        # Check if the file is not empty
        if os.path.getsize(file_path) > 0:
            return True
        else:
            return False
    else:
        return False


def merge_pdfs(output_dir, output_file_n, *pdfs):
    """
    Merge PDF files into one single PDF file. 
    The order in which file paths are provided in <*pdfs> defines the order of plots in the final single PDF file.

    """
    merger = PdfMerger()
    valid_pdf_count = 0  # Counter to keep track of valid PDFs added

    # Iterate over the list of file paths
    for pdf in pdfs:
        # Check if the file exists and is not empty before appending
        if pdf is not None and os.path.exists(pdf) and os.path.getsize(pdf) > 0:
            # Append PDF files
            merger.append(pdf)
            valid_pdf_count += 1

    if valid_pdf_count > 0:
        merged_pdf_path = os.path.join(output_dir, os.path.join(output_dir, f"{output_file_n}_me.pdf"))
        merger.write(merged_pdf_path)
        merger.close()
        return merged_pdf_path

    if valid_pdf_count == 0:
        merger.close()
        return False


def dotplot(sequence1, sequence2, output_dir):

    # Define based on TEtrimmer analysis file name
    n_after_tetrimmer = os.path.basename(sequence1)

    # Define the output filenames
    pdf_out = os.path.join(output_dir, f"{n_after_tetrimmer}.ps.pdf")

    # Define command for dotmatcher
    dotmatcher_command = [
        "dotmatcher",
        "-asequence", sequence2,
        "-bsequence", sequence1,
        "-windowsize", "25",
        "-threshold", "50",
        "-gtitle", "Dotplot",
        "-gxtitle", "After TEtrimmer treatment",
        "-gytitle", "Before TEtrimmer treatment",
        "-graph", "ps",
        "-goutfile", sequence1
    ]

    try:
        subprocess.run(dotmatcher_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    except FileNotFoundError:
        prcyan("\n'dotmatcher' command not found. Please ensure 'emboss' is installed correctly.")
        prgre("\ndotmatcher does not affect the final consensus sequence. You can choose to ignore this error.\n")
        return None

    except subprocess.CalledProcessError as e:
        prcyan(f"\n'dotmatcher' failed for {n_after_tetrimmer} with error code {e.returncode}")
        prgre("\ndotmatcher does not affect the final consensus sequence. You can choose to ignore this error.")
        return None

    # Define command to convert ps to pdf
    ps2pdf_command = [
        "ps2pdf",
        f"{sequence1}.ps",
        pdf_out
    ]

    try:
        subprocess.run(ps2pdf_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    except FileNotFoundError:
        prcyan("\n'ps2pdf' command not found. Please install it with 'sudo apt-get install ghostscript'")
        prgre("ps2pdf does not affect the final consensus sequence. You can choose to ignore this error.")
        return None

    except subprocess.CalledProcessError as e:
        prcyan(f"\n'ps2pdf' failed for {n_after_tetrimmer} with error code {e.returncode}")
        prgre("\nps2pdf does not affect the final consensus sequence. You can choose to ignore this error.\n")
        return None

    return pdf_out


def multi_seq_dotplot(input_file, output_dir, title):

    input_name = os.path.basename(input_file)

    # Define the output filenames
    ps_out = os.path.join(output_dir, f"{title}.ps")
    pdf_out = os.path.join(output_dir, f"0TEtrimmer_{title}_multiple_sequence_dotplot.pdf")

    # Define command for dotmatcher
    polydot_command = [
        "polydot",
        str(input_file),
        "-wordsize", str(40),
        "-gtitle", str(title),
        "-gdirectory", str(output_dir),
        "-goutfile", str(title),
        "-graph", "ps"
    ]

    try:
        subprocess.run(polydot_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    except FileNotFoundError:
        prcyan("'\npolydot' command not found. Please ensure 'emboss' is correctly installed.")
        prgre("\npolydot won't affect the final consensus sequence. You can choose to ignore this error\n")
        raise Exception

    except subprocess.CalledProcessError as e:
        prcyan(f"\npolydot failed for {input_name} with error code {e.returncode}")
        prgre("\npolydot won't affect the final consensus sequence. You can choose to ignore this error")
        raise Exception

    # Define command to convert ps to pdf
    ps2pdf_command = [
        "ps2pdf",
        str(ps_out),
        str(pdf_out)
    ]

    try:
        subprocess.run(ps2pdf_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    except FileNotFoundError:
        prcyan("'\nps2pdf' command not found. Please install it by 'sudo apt-get install ghostscript'")
        prgre("ps2pdf won't affect the final consensus file. You can choose to ignore it.")
        raise Exception

    except subprocess.CalledProcessError as e:
        prcyan(f"\nps2pdf failed for {input_name} with error code {e.returncode}")
        prgre("\nps2pdf won't affect the final consensus sequence. You can choose to ignore this error\n")
        raise Exception

    return pdf_out


def scale_single_page_pdf(input_pdf_path, output_pdf_path, scale_ratio):
    pdf_reader = PdfFileReader(input_pdf_path)
    pdf_writer = PdfFileWriter()

    page = pdf_reader.getPage(0)  # Get the first (and only) page
    page.scale_by(scale_ratio)  # Scale the page
    pdf_writer.addPage(page)

    with open(output_pdf_path, 'wb') as out_file:
        pdf_writer.write(out_file)

    return output_pdf_path


# Function used to find and return poly A end position
def find_poly_a_end_position(input_file, min_length=10):

    try:
        # Read input file and get sequence length
        record = SeqIO.read(input_file, "fasta")
        seq_length = len(record)
        seq_mid_position = seq_length // 2

        reverse_sequence = record[::-1].upper()
        poly_a_length = 0
        interruption_allowed = True  # Allow for one non 'A' interruption

        for i, nucleotide in enumerate(reverse_sequence):
            # Only check the second half of the sequence
            if i < seq_mid_position:
                if nucleotide == 'A':
                    poly_a_length += 1
                elif interruption_allowed and 5 <= poly_a_length < min_length:
                    # Allow for one interruption if we have encountered at least 5 'A's
                    interruption_allowed = False  # Use up the allowance for interruption
                    poly_a_length += 1  # Count the interruption as an 'A'
                else:
                    if poly_a_length >= min_length:
                        # Correctly calculate the end position in the original orientation
                        return seq_length - i - 1 + poly_a_length
                    poly_a_length = 0  # Reset counter if sequence is interrupted
                    interruption_allowed = True  # Reset the allowance for an interruption
        return None  # No poly A sequence found or does not meet the minimum length requirement
    except Exception as e:
        prcyan(f"\nPoly A detection failed for {os.path.basename(input_file)} with error:\n{e}")
        prcyan('\n' + 'This is will not affect the result too much, you can choose ignore this error' + '\n')
        return None


def calculate_con_coverage_num(consensus_file, blast_file):

    # Use Biopython to read the consensus sequence from the FASTA file
    record = SeqIO.read(consensus_file, "fasta")
    consensus_seq = record.seq
    cons_len = len(consensus_seq)

    # Initialize a list with zeros for each position in the consensus sequence
    coverage = [0] * cons_len

    # Process BLAST hits from the file
    with open(blast_file, 'r') as f:
        for line in f:
            parts = line.split()
            start, end = int(parts[6]), int(parts[7])
            # Adjust for Python's 0-based indexing and ensure the positions are within the consensus length
            start, end = max(1, start), min(end, cons_len)
            for position in range(start - 1, end):  # Adjusting start index to 0-based for Python
                coverage[position] += 1

    return coverage
