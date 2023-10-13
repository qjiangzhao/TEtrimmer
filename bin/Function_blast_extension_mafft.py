import shutil
import subprocess
import os
import click
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment, AlignInfo
import pandas as pd
import pandas.errors


def calculate_genome_length(genome_file):
    """
    Calculate the length of each sequence in a genome file in FASTA format
    and write the lengths to an output file.

    :param genome_file: str, path to genome file in FASTA format
    :return: str, path to the output file containing sequence names and lengths
    """
    genome_lengths = {}
    with open(genome_file, "r") as f:
        current_seq = None
        current_length = 0
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_seq is not None:
                    genome_lengths[current_seq] = current_length
                # If chromosome header contain space, only consider the content before the first space
                current_seq = line[1:].split(" ")[0]
                current_length = 0
            else:
                current_length += len(line)
        # Add length of last sequence to dictionary
        if current_seq is not None:
            genome_lengths[current_seq] = current_length

    # Write lengths to output file
    output_file = genome_file + ".length"
    with open(output_file, "w") as out:
        for seq_name, length in genome_lengths.items():
            out.write(f"{seq_name}\t{length}\n")

    return output_file


def check_database(genome_file):
    """
    Checks if the blast database and genome length file exist.
    If they don't exist, create them.

    :param genome_file: str, path to genome file (contain genome name)
    """

    blast_database_file = genome_file + ".nin"
    if not os.path.isfile(blast_database_file):
        print(f"\nBlast database doesn't exist. Running makeblastdb\n")

        makeblastdb_cmd = f"makeblastdb -in {genome_file} -dbtype nucl -out {genome_file} "
        result = subprocess.run(makeblastdb_cmd, shell=True, check=True, stderr=subprocess.PIPE)
        error_output = result.stderr.decode("utf-8")

        if error_output:
            print(f"Error check_database {error_output}\n")

    # Check if genome length files exists, otherwise create it at the same folder with genome file
    length_file = genome_file + ".length"
    if not os.path.isfile(length_file):
        print(f"\nFile with genome lengths not found. Making it now\n")
        calculate_genome_length(genome_file)

    # Check if .fai index file exists, otherwise create it using bedtools
    fai_file = genome_file + ".fai"
    if not os.path.isfile(fai_file):
        print(f"\nIndex file {fai_file} not found. Creating it by samtools faidx\n")
        faidx_cmd = f"samtools faidx {genome_file}"
        result_fai = subprocess.run(faidx_cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        error_output_fai = result_fai.stderr.decode("utf-8")

        if error_output_fai:
            print(f"Error check_database {error_output_fai}\n")


def blast(seq_obj, seq_file, genome_file, output_dir, min_length=150):
    """
    Runs blastn and saves the results in a bed file.

    :param input_file: str, path to input fasta file
    :param genome_file: str, path to genome file
    :param output_dir: str, prefix for output files
    :param min_length: num default 150, minimum alignment length

    :return: a bed output file name
    """
    input_file = seq_file
    blast_hits_count = 0
    bed_out_file = None
    # define blast outfile
    blast_out_file = os.path.join(output_dir, f"{os.path.basename(input_file)}.b")
    blast_cmd = f"blastn -query {input_file} -db {genome_file} " \
                f"-outfmt \"6 qseqid sseqid pident length mismatch qstart qend sstart send sstrand\" " \
                f"-evalue 1e-40 -qcov_hsp_perc 20 | " \
                f"awk -v 'ml={min_length}' 'BEGIN{{OFS=\"\\t\"}} $4 > ml {{print $0}}' >> {blast_out_file}"
    result = subprocess.run(blast_cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    error_output = result.stderr.decode("utf-8")

    if error_output:
        print(f"Error blast {error_output}")

    # Check the number of blast hits
    with open(blast_out_file) as blast_file:
        for _ in blast_file:
            blast_hits_count += 1

    if blast_hits_count >= 1:
        # Define bed outfile
        # Only convert blast to bed file when hits number is greater than 10
        bed_out_file = os.path.join(output_dir, f"{os.path.basename(input_file)}.b.bed")
        bed_cmd = f"awk 'BEGIN{{OFS=\"\\t\"}} !/^#/ {{if ($10~/plus/){{print $2, $8, $9, $1, $3, \"+\"}} " \
                  f"else {{print $2, $9, $8, $1, $3, \"-\"}}}}' < {blast_out_file} > {bed_out_file}"
        result_awk = subprocess.run(bed_cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        error_output_awk = result_awk.stderr.decode("utf-8")

        if error_output_awk:
            print(f"Error blast {error_output_awk}")
        
        # Update blast hit number to sequence object
        seq_obj.update_blast_hit_n(blast_hits_count)
        
    # if low copy number, check coverage, length and identity
    return bed_out_file, blast_hits_count, blast_out_file


def check_bed_uniqueness(output_dir, bed_file):
    """
    Check if bed file contain repeat lines, if so delete to one

    :param output_dir: str, output directory
    :param bed_file: str, bed_file path
    :return: modified bed file absolute path
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


# Calculate average sequence length in the bed file
def bed_ave_sequence_len(bed_content, start_rank, end_rank):
    """Compute average length for regions within a specified rank range in BED content."""

    # Extracting lengths from the BED content
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
        return "Not enough entries"

    # Extracting lengths of regions within the specified rank range
    selected_lengths = lengths[start_rank - 1:end_rank]

    # Calculating average
    average = sum(selected_lengths) / len(selected_lengths)

    return average


def extract_fasta(input_file, genome_file, output_dir, left_ex, right_ex):
    """
    Extracts fasta sequence from the reference genome using bedtools.

    :param genome_file: str, path to genome file
    :param output_dir: str, prefix for output files
    :param left_ex: int, number of bases to extend the start position of each feature
    :param right_ex: int, number of bases to extend the end position of each feature
    :return: fasta file absolute path derived from bed file
    """
    bed_out_flank_file_dup = os.path.join(output_dir, f"{os.path.basename(input_file)}_{left_ex}_{right_ex}.bed")

    fasta_out_flank_file = os.path.join(output_dir, f"{os.path.basename(input_file)}_{left_ex}_{right_ex}.fa")

    fasta_out_flank_file_nucleotide_clean = os.path.join(output_dir,
                                                            f"{os.path.basename(input_file)}_{left_ex}_{right_ex}_cl.fa")

    bed_cmd = f"bedtools slop -s -i {input_file} -g {genome_file}.length -l {left_ex} -r {right_ex} > {bed_out_flank_file_dup}"
    subprocess.run(bed_cmd, shell=True, check=True)

    bed_out_flank_file = check_bed_uniqueness(output_dir, bed_out_flank_file_dup)

    fasta_cmd = f"bedtools getfasta -s -fi {genome_file} -fo {fasta_out_flank_file} -bed {bed_out_flank_file}"
    subprocess.run(fasta_cmd, shell=True, check=True)

    # Use awk to remove letters that aren't A G C T a g c t
    fasta_nucleotide_clean = f"awk '/^>/ {{print}} !/^>/ {{gsub(/[^AGCTagct]/, \"\"); print}}' {fasta_out_flank_file} > {fasta_out_flank_file_nucleotide_clean}"
    subprocess.run(fasta_nucleotide_clean, shell=True, check=True)

    return fasta_out_flank_file_nucleotide_clean, bed_out_flank_file_dup


def align_sequences(input_file, output_dir):
    """
    Aligns fasta sequences using MAFFT

    :param input_file: str, input file path
    :param output_dir: str, output file directory
    :return: multiple sequence alignment file absolute path
    """
    fasta_out_flank_mafft_file = os.path.join(output_dir, f"{os.path.basename(input_file)}_aln.fa")

    # Construct the command as a list of strings
    mafft_cmd = ["mafft", "--quiet", "--nuc", "--retree", "1", input_file]

    # Execute the command
    result = subprocess.run(mafft_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Check for errors
    if result.returncode != 0:
        raise Exception(f"MAFFT failed with error: {os.path.basename(input_file)}\n{result.stderr.decode('utf-8')}")

    # Write the output to the file
    with open(fasta_out_flank_mafft_file, 'w') as f:
        f.write(result.stdout.decode('utf-8'))

    return fasta_out_flank_mafft_file


def muscle_align(input_file, output_dir):

    output_file = os.path.join(output_dir, f"{os.path.basename(input_file)}_maln.fa")
    muscle_cmd = ["muscle", "-maxiters", "4", "-in", input_file, "-out", output_file]

    # Execute muscle
    result = subprocess.run(muscle_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if result.returncode != 0:
        return False
        # raise Exception(f"Muscle command failed with error: {os.path.basename(input_file)}\n{result.stderr.decode('utf-8')}")

    # Convert the sequences in the output file to lowercase
    sequences = list(SeqIO.parse(output_file, "fasta"))
    for seq in sequences:
        seq.seq = seq.seq.lower()
    SeqIO.write(sequences, output_file, "fasta")

    return output_file


def con_generater(input_file, output_dir, threshold=0.8, ambiguous="N"):
    # Read input file
    alignment = AlignIO.read(input_file, "fasta")
    summary = AlignInfo.SummaryInfo(alignment)
    # Get consensus sequence
    consensus_seq = summary.dumb_consensus(threshold=threshold, ambiguous=ambiguous).upper()

    consensus_record = SeqRecord(consensus_seq, id=f"{os.path.basename(input_file)}", description="")

    # Write to a FASTA file
    output_file = os.path.join(output_dir, f"{os.path.basename(input_file)}_co.fa")
    with open(output_file, "w") as file:
        SeqIO.write(consensus_record, file, "fasta")

    return output_file


def con_generater_no_file(input_file, threshold=0.8, ambiguous="N"):
    # Read input file
    alignment = AlignIO.read(input_file, "fasta")
    summary = AlignInfo.SummaryInfo(alignment)
    # Get consensus sequence
    consensus_seq_str = summary.dumb_consensus(threshold=threshold, ambiguous=ambiguous)

    # Return the consensus sequence string
    return consensus_seq_str


def calc_conservation(col):
    """
    Calculate the conservation of a column as the fraction of the most common nucleotide.
    :param col: one single column content for multiple sequence alignment
    :return: most abundant nucleotide proportion for this column (gap isn't included)
    """
    # Mafft will convert all nucleotide to lowercase
    nucleotide_counts = {'a': 0, 'c': 0, 'g': 0, 't': 0}

    for nucleotide in col:
        if nucleotide in nucleotide_counts:
            nucleotide_counts[nucleotide] += 1

    total_nucleotides = sum(nucleotide_counts.values())
    max_count = max(nucleotide_counts.values())

    return max_count / total_nucleotides


def generate_hmm_from_msa(input_msa_file, output_hmm_file):
    """
    Generate HMM profile using hmmbuild from a multiple sequence alignment file.

    :param input_msa_file: Path to the multiple sequence alignment file (in FASTA or Stockholm format).
    :param output_hmm_file: Path to the output HMM file.
    """

    # Construct the command as a list
    cmd = ["hmmbuild", "--dna", output_hmm_file, input_msa_file]

    try:
        # Execute the command
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Check if the command was successful
        if result.returncode != 0:

            click.echo(f"Error running hmmbuild: {os.path.basename(output_hmm_file)}\n{result.stderr.decode('utf-8')}")

    except FileNotFoundError:
        click.echo("Error: hmmbuild not found. Ensure HMMER is installed and available in the PATH.")


def reverse_complement_seq_file(input_file, output_file):
    """
    Takes an MSA FASTA file and writes the reverse complemented sequences to a new file.
    """
    reverse_complemented_records = []

    for record in SeqIO.parse(input_file, "fasta"):
        # Create a new record with the reverse complemented sequence
        rev_comp_seq = record.seq.reverse_complement()
        rev_comp_record = SeqRecord(rev_comp_seq, id=record.id, description="")
        reverse_complemented_records.append(rev_comp_record)

    # Write the reverse complemented sequences to the output file
    SeqIO.write(reverse_complemented_records, output_file, "fasta")
    return output_file


def remove_gaps(input_file, output_dir, threshold=0.8, min_nucleotide=5):
    """
    Remove gaps when gap percentage is bigger than threshold. Remove columns when nucleotide number is less than 5
    :param input_file: str, input file path
    :param output_dir: str, output file directory
    :param threshold: num (0-1) default 0.8, column with gap percentage that is greater (not include equal)
    than threshold will be removed
    :param min_nucleotide: num (>0) default 5, columns with less than "min_nucleotide" nucleotides will be removed
    :return: absolut path of gap removed multiple sequence alignment file
    """
    keep_list = []
    fasta_out_flank_mafft_gap_filter_file = os.path.join(output_dir,
                                                            f"{os.path.basename(input_file)}_g.fa")
    MSA_mafft = AlignIO.read(input_file, "fasta")

    column_mapping = {}  # Stores the mapping of column indices from filtered MSA to original MSA

    for col_idx in range(MSA_mafft.get_alignment_length()):
        col = MSA_mafft[:, col_idx]
        gap_count = col.count('-')
        gap_fraction = gap_count / len(col)

        # Check if the count of nucleotides in the column is at least 5
        nt_count = len(col) - gap_count
        if nt_count < min_nucleotide:
            continue

        # If the gap fraction is less than the threshold, add the column to the filtered MSA
        if gap_fraction <= threshold:
            keep_list.append(col_idx)
            # Be careful when you want to use len for indexing. Because len start from 1
            column_mapping[
                len(keep_list) - 1] = col_idx  # Store the mapping of original MSA index to filtered MSA index

    # The index in python won't include the second value. That it is to say the returned end_posit from
    # DefineBoundary() will one more index than the final len(keep_list) -1] for this reason, you have to
    # add one more value
    column_mapping[len(keep_list)] = column_mapping[len(keep_list) - 1] + 1
    # Keep the columns when they meet requirements
    MSA_mafft_filtered = MSA_mafft[:, keep_list[0]:keep_list[0] + 1]
    for i in keep_list[1:]:
        MSA_mafft_filtered += MSA_mafft[:, i:i + 1]

    # Write the filtered MSA to the output file
    with open(fasta_out_flank_mafft_gap_filter_file, 'w') as f:
        AlignIO.write(MSA_mafft_filtered, f, 'fasta')
    return fasta_out_flank_mafft_gap_filter_file, column_mapping


def remove_gaps_block(input_file, output_dir, threshold=0.8, conservation_threshold=0.5, min_nucleotide=5):
    """
    Remove gaps that are only from conserved regions. Remove columns when nucleotide number is less than 5
    :param input_file: str, input file path
    :param output_dir: str, output file directory
    :param threshold: num (0-1) default 0.8, column with gap percentage that is greater (not include equal)
    than threshold will be removed
    :param conservation_threshold: num (0-1) default 0.5, connected single gaps will be classified as one block, two columns
    before and after this block will be checked. If the most abundant nucleotides proportion from one beginning and
    one end columns are greater than conservation_threshold this gap block will be removed
    :param min_nucleotide: num (>0) default 5, columns with less than "min_nucleotide" nucleotides will be removed
    :return: gap block removed file absolute path
    """
    # The keep_list is a list of tuples containing start and end position for each gap block
    gap_blocks = []
    fasta_out_flank_mafft_gap_filter_file = os.path.join(output_dir,
                                                            f"{os.path.basename(input_file)}_gb.fa")
    MSA_mafft = AlignIO.read(input_file, "fasta")

    # Initialize the start of the block to None
    block_start = None

    for col_idx in range(MSA_mafft.get_alignment_length()):
        col = MSA_mafft[:, col_idx]
        gap_count = col.count('-')
        gap_fraction = gap_count / len(col)
        # Check the count of nucleotides in the column
        nt_count = len(col) - gap_count

        # If the column meets the gap requirement, and we're not already in a block, start a new one
        if (gap_fraction > threshold or nt_count < min_nucleotide) and block_start is None:
            block_start = col_idx
        # If the column doesn't meet the gap requirement, and we're in a block, end the block
        elif gap_fraction <= threshold and nt_count >= min_nucleotide and block_start is not None:
            gap_blocks.append((block_start, col_idx))
            block_start = None

    # If we're still in a block at the end of the MSA, end it
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

    # When no gap block is detected, return the original MSA
    else:
        # No blocks to delete, return the original MSA. Write file again to convert file name
        with open(fasta_out_flank_mafft_gap_filter_file, 'w') as f:
            AlignIO.write(MSA_mafft, f, 'fasta')
        return fasta_out_flank_mafft_gap_filter_file

    if delete_blocks:

        # Identify the ranges to keep, which are in between the blocks to delete
        keep_ranges = [(0, delete_blocks[0][0])]
        for i in range(len(delete_blocks) - 1):
            keep_ranges.append((delete_blocks[i][1], delete_blocks[i + 1][0]))
        keep_ranges.append((delete_blocks[-1][1], MSA_mafft.get_alignment_length()))

        # Concatenate the columns in each keep_range to create the final filtered MSA
        MSA_mafft_filtered = MSA_mafft[:, keep_ranges[0][0]:keep_ranges[0][1]]
        for keep_range in keep_ranges[1:]:
            MSA_mafft_filtered += MSA_mafft[:, keep_range[0]:keep_range[1]]

        # Write the filtered MSA to the output file
        with open(fasta_out_flank_mafft_gap_filter_file, 'w') as f:
            AlignIO.write(MSA_mafft_filtered, f, 'fasta')

        return fasta_out_flank_mafft_gap_filter_file
    else:
        # No blocks to delete, return the original MSA. Write file again to convert file name
        with open(fasta_out_flank_mafft_gap_filter_file, 'w') as f:
            AlignIO.write(MSA_mafft, f, 'fasta')
        return fasta_out_flank_mafft_gap_filter_file


def remove_gaps_with_similarity_check(input_file, output_dir, gap_threshold=0.8,
                                        simi_check_gap_thre=0.4, similarity_threshold=0.7, min_nucleotide=5):
    """
    Remove gaps when gap percentage is bigger than threshold. Remove columns when nucleotide number is less than 5.
    When gap percentage is equal or bigger than "simi_check_gap_thre". It will calculate most abundant nucleotide
    proportion for this column. If it is less than "similarity_threshold" this column will be removed
    :param input_file: str, input file path
    :param output_dir: str, output file directory
    :param gap_threshold: num (0-1) default 0.8, column with gap percentage that is greater (not include equal)
    than threshold will be removed
    :param simi_check_gap_thre: num (0-1) default 0.4, if gap greater or equal this gap threshold, column similarity
    will be checked
    :param similarity_threshold: num (0-1) default 0.7, when column gap is between "simi_check_gap_thre" and
    "gap_threshold" and most abundant nucleotide proportion is less than "similarity_threshold" remove this column
    :param min_nucleotide: num (>0) default 5, columns with less than "min_nucleotide" nucleotides will be removed
    :return: gap removed file absolute path
    """
    keep_list = []
    fasta_out_flank_mafft_gap_filter_file = os.path.join(output_dir,
                                                            f"{os.path.basename(input_file)}_gs.fa")
    MSA_mafft = AlignIO.read(input_file, "fasta")

    column_mapping = {}  # Stores the mapping of column indices from filtered MSA to original MSA

    if len(MSA_mafft) < 5:
        raise ValueError("Number of sequences is less than 5. Cannot remove gaps.")

    for col_idx in range(MSA_mafft.get_alignment_length()):
        col = MSA_mafft[:, col_idx]
        gap_count = col.count('-')
        gap_fraction = gap_count / len(col)

        # Check if the number of nucleotides in a column is less than 5
        nt_count = len(col) - gap_count
        if nt_count < min_nucleotide:
            continue

        # If the gap fraction is less than the gap_threshold, consider the column for further analysis
        if gap_fraction <= gap_threshold:
            # If the gap fraction is between 0.4 and 0.8, check for the similarity of the remaining nucleotides
            if simi_check_gap_thre > gap_fraction:
                keep_list.append(col_idx)
                column_mapping[len(keep_list) - 1] = col_idx
            elif simi_check_gap_thre <= gap_fraction:
                nt_fraction = calc_conservation(col)

                # If the nucleotides are less than similarity_threshold, skip this column
                if nt_fraction < similarity_threshold:
                    continue
                else:
                    # If the column passes all checks, add it to the keep_list
                    keep_list.append(col_idx)
                    # Store the mapping of original MSA index to filtered MSA index
                    column_mapping[len(keep_list) - 1] = col_idx

    # The index in python won't include the second value. That it is to say the returned end_posit from
    # DefineBoundary() will one more index than the final len(keep_list) -1] for this reason, you have to
    # add one more value
    column_mapping[len(keep_list)] = column_mapping[len(keep_list) - 1] + 1

    # Keep the columns
    MSA_mafft_filtered = MSA_mafft[:, keep_list[0]:keep_list[0] + 1]
    for i in keep_list[1:]:
        MSA_mafft_filtered += MSA_mafft[:, i:i + 1]

    # Write the filtered MSA to the output file
    with open(fasta_out_flank_mafft_gap_filter_file, 'w') as f:
        AlignIO.write(MSA_mafft_filtered, f, 'fasta')
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
    :param gap_threshold: num (0-1) default 0.8, column with gap percentage that is greater (not include equal)
    than threshold will be removed
    :param simi_check_gap_thre: num (0-1) default 0.4, if gap greater or equal this gap threshold, column similarity
    will be checked
    :param similarity_threshold: num (0-1) default 0.7, when column gap is between "simi_check_gap_thre" and
    "gap_threshold" and most abundant nucleotide proportion is less than "similarity_threshold" remove this column
    :param conservation_threshold: num (0-1) default 0.6, connected single gaps will be classified as one block, two columns
    before and after this block will be checked. If the most abundant nucleotides proportion from one beginning and
    one end columns are greater than conservation_threshold this gap block will be removed
    :param min_nucleotide: num (>0) default 5, columns with less than "min_nucleotide" nucleotides will be removed
    :return: gap block removed file absolute path
    """
    # Identify blocks of gaps
    gap_blocks = []

    # Output file name
    fasta_out_flank_mafft_gap_filter_file = os.path.join(output_dir, f"{os.path.basename(input_file)}_gbs.fa")

    # Load the MSA file
    MSA_mafft = AlignIO.read(input_file, "fasta")

    # Raise an error if the number of sequences is less than 5
    if len(MSA_mafft) < 5:
        raise ValueError("Number of sequences is less than 5. Cannot remove gaps.")

    # Define the starting point of a block
    block_start = None

    # Go through each column in the alignment
    for col_idx in range(MSA_mafft.get_alignment_length()):

        # Define some metrics for the column
        col = MSA_mafft[:, col_idx]
        gap_count = col.count('-')
        gap_fraction = gap_count / len(col)
        nt_fraction = calc_conservation(col)
        nt_count = len(col) - gap_count

        # Check if the column meets the criteria to start a new block
        if ((simi_check_gap_thre <= gap_fraction <= gap_threshold and nt_fraction <= similarity_threshold) or
            gap_fraction > gap_threshold or nt_count < min_nucleotide) and block_start is None:
            block_start = col_idx

        # Check if the column meets the criteria to end a block
        elif simi_check_gap_thre <= gap_fraction <= gap_threshold and nt_fraction > \
                similarity_threshold and nt_count >= min_nucleotide and block_start is not None:
            gap_blocks.append((block_start, col_idx))
            block_start = None

        # Check if the block should end due to a column with too few gaps
        elif gap_fraction <= simi_check_gap_thre and nt_count >= min_nucleotide and block_start is not None:
            gap_blocks.append((block_start, col_idx))
            block_start = None

    # If a block was started but not ended, end it at the last column
    if block_start is not None:
        gap_blocks.append((block_start, MSA_mafft.get_alignment_length()))

    # Identify the blocks that should be deleted based on the conservation of the surrounding columns
    delete_blocks = []

    # Check if gap_blocks is empty
    if gap_blocks:
        for block in gap_blocks:
            # Check the columns divergence before and after the block
            col_before = MSA_mafft[:, block[0] - 1] if block[0] > 0 else None
            col_before_2 = MSA_mafft[:, block[0] - 2] if block[0] > 1 else None
            col_after = MSA_mafft[:, block[1]] if block[1] < MSA_mafft.get_alignment_length() else None
            col_after_2 = MSA_mafft[:, block[1] + 1] if block[1] < MSA_mafft.get_alignment_length() - 1 else None
            # If the surrounding columns are conserved, mark the block for deletion
            if col_before is not None and (calc_conservation(col_before) > conservation_threshold or
                                            calc_conservation(
                                                col_before_2) > conservation_threshold) and col_after is not None and \
                    (calc_conservation(col_after) > conservation_threshold or
                        calc_conservation(col_after_2) > conservation_threshold):
                delete_blocks.append(block)
    else:
        # No blocks to delete, return the original MSA. Write file again to convert file name
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
    Return: A list containing column numbers used for cluster
    """
    # Identify blocks of gaps
    gap_blocks = []

    # Load the MSA file
    MSA_mafft = AlignIO.read(input_file, "fasta")

    # Raise an error if the number of sequences is less than 5
    if len(MSA_mafft) < 5:
        raise ValueError("Number of sequences is less than 5. Cannot remove gaps.")

    # Define the starting point of a block
    block_start = None

    # Variables to keep track of the gap count in the previous column and the condition met count
    prev_gap_count = None
    condition_met_count = 0

    # Go through each column in the alignment. Start from 0
    for col_idx in range(MSA_mafft.get_alignment_length()):
        # Define some metrics for the column
        col = MSA_mafft[:, col_idx]
        gap_count = col.count('-')
        gap_fraction = gap_count / len(col)
        nt_fraction = calc_conservation(col)
        nt_count = len(col) - gap_count

        # Check if the column meets the criteria to start a new block
        if simi_check_gap_thre <= gap_fraction and nt_fraction >= similarity_threshold and \
                nt_count >= min_nucleotide and block_start is None:
            block_start = col_idx
            condition_met_count = 0  # Reset the count for a new block

        # Additional check for gap variance
        elif prev_gap_count is not None and (
                gap_count < 0.8 * prev_gap_count or gap_count > 1.2 * prev_gap_count) and block_start is not None:
            gap_blocks.append((block_start, col_idx, condition_met_count))
            block_start = None

        # Stop the gap block when nucleotide similarity is smaller than the threshold
        elif (
                simi_check_gap_thre <= gap_fraction and nt_fraction < similarity_threshold and block_start is not None) or \
                (nt_count < min_nucleotide and block_start is not None):
            # Check the gap count in the next column if it exists
            next_gap_count = None
            if col_idx + 1 < MSA_mafft.get_alignment_length():
                next_col = MSA_mafft[:, col_idx + 1]
                next_gap_count = next_col.count('-')

            # If gap count is not similar to both previous and next columns, stop the block
            if prev_gap_count is not None and next_gap_count is not None:
                if not (0.8 * prev_gap_count <= gap_count <= 1.2 * prev_gap_count and
                        0.8 * next_gap_count <= gap_count <= 1.2 * next_gap_count):
                    condition_met_count += 1

        # Stop the gap block when the gap number is too less
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

            # Calculate the maximum allowable condition_met_count based on the block length
            block_length = end - start
            if block_length < 10:
                max_condition_count = 1
            else:
                max_condition_count = min(0.1 * block_length, 50)

            # Check if the condition_met_count exceeds the maximum allowable value for the block
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


def select_start_end_and_join( input_file, output_dir, start, end, window_size=50):
    """
    Select start and end columns of MSA
    :param input_file: str, absolute input file path
    :param output_dir: str, absolute output directory
    :param start: int: start point, based on start point, right side columns will be selected
    :param end: int: end point, based on end point, left side columns will be selected
    :param window_size: int default 50, columns size will be selected for each side
    :return: Selected MSA object (not file path)
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
    Select columns block for MSA
    :param input_file: str, the absolute input file path
    :param output_dir: str, the absolute output directory
    :param start_point: int, start point to select column block
    :param direction: int "right" or "left", direction to select column block
    :param window_size: int default 50, column block size will be selected
    :return: Select MSA object (not file path)
    """
    # Read the MSA file using Biopython's AlignIO
    alignment = AlignIO.read(input_file, "fasta")

    # Get the total number of columns in the alignment
    total_columns = alignment.get_alignment_length()

    # Determine the start and end column indices based on the direction
    if direction.lower() == 'left':
        start_col = max(0, start_point - window_size)
        end_col = start_point
    elif direction.lower() == 'right':
        start_col = start_point
        end_col = min(start_point + window_size, total_columns)  # get as many columns as possible
    else:
        raise ValueError("Invalid direction. Please choose 'left' or 'right'.")

    # Select the window columns from the alignment
    selected_alignment = alignment[:, start_col:end_col]

    # Convert selected_alignment to a MultipleSeqAlignment object
    selected_alignment = MultipleSeqAlignment(selected_alignment)

    output_file = os.path.join(output_dir, f"{os.path.basename(input_file)}_{direction}_{window_size}.fa")
    AlignIO.write(selected_alignment, output_file, "fasta")
    # Return the new alignment
    return selected_alignment

"""
The asterisk (*) before alignments in the function definition is used to allow the function to accept any number 
of positional arguments. These arguments will be gathered into a tuple called alignments. This is particularly 
useful when you don't know beforehand how many arguments will be passed to the function.
"""


def concatenate_alignments(*alignments, input_file_name, output_dir):
    """
    Concatenate MSA files, which have same sequence number
    :param alignments: str, alignments object (not file path) will be concatenated by order.
    :param input_file_name: str, Used for generating output file name
    :param output_dir: str, the absolute output directory
    :return: the absolute output file path, the start point and end point of the concatenated MSA
    """
    # Check if at least two alignments are provided
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

    # alignments[0].get_alignment_length() will return 50, but python start from 0
    # The last position of alignments[0] is alignments[0].get_alignment_length() - 1
    concat_start = alignments[0].get_alignment_length()
    concat_end = concatenated_alignment.get_alignment_length() - alignments[-1].get_alignment_length() - 1

    # Write to a file
    output_file = os.path.join(output_dir, f"{os.path.basename(input_file_name)}_me.fa")
    AlignIO.write(concatenated_alignment, output_file, "fasta")

    # Return the concatenated alignment
    return output_file, concat_start, concat_end


def change_permissions_recursive(input_dir, mode):
    try:
        for dirpath, dirnames, filenames in os.walk(input_dir):
            os.chmod(dirpath, mode)
            for filename in filenames:
                os.chmod(os.path.join(dirpath, filename), mode)
    except PermissionError:
        click.echo("TE Trimmer don't have right to change permissions. Pleas use sudo to run TE Trimmer")
        return False
    return True


def cd_hit_est(input_file, output_file, identity_thr=0.8, aL=0.9, aS=0.9, s=0.9, thread=10):

    click.echo("cd-hit-est is running\n")

    command = [
        "cd-hit-est",
        "-i", input_file,
        "-o", output_file,
        "-c", str(identity_thr),
        "-aL", str(aL),
        "-aS", str(aS),
        "-M", "0",
        "-T", str(thread),
        "-l", "80",
        "-d", "0",
        "-s", str(s)
    ]

    result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    error_output = result.stderr.decode("utf-8")

    if error_output:

        click.echo(f"Error executing cd-hit-est command {error_output}\n")
        return False

    else:
        click.echo("cd-hit-est is finished")
        return True


def repeatmasker(genome_file, library_file, output_dir, thread=1, classify=False):
    """
    Run RepeatMasker with the provided parameters.
    """

    # Construct the RepeatMasker command
    if classify:
        command = ["RepeatMasker",
                   genome_file,
                   "-lib", library_file,
                   "-s",    # Slow search; 0-5% more sensitive, 2-3 times slower than default
                   "-dir", output_dir,
                   "-pa", str(thread)
                    ]
    else: 
        command = ["RepeatMasker",
                    genome_file,
                    "-lib", library_file,
                    "-pa", str(thread),
                    "-dir", output_dir,
                    "-s",    # Slow search; 0-5% more sensitive, 2-3 times slower than default
                    "-gff",  # Creates an additional Gene Feature Finding format output
                    "-xm",   # Creates an additional output file in cross_match format (for parsing)
                    "-a",    # Writes alignments in .align output file
                    ]

    result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    error_output = result.stderr.decode("utf-8")

    if error_output:

        click.echo(f"Error executing RepeatMasker command{error_output}\n.")
        return False

    else:
        return True


def repeatmasker_output_classify(repeatmasker_out, progress_file, min_iden=70, min_len=80, min_cov=0.8):

    # Read RepeatMasker out file into a DataFrame
    # The regex '\s+' matches one or more whitespace characters
    try:
        df = pd.read_csv(repeatmasker_out, delim_whitespace=True, header=None, skiprows=3)
    except pandas.errors.EmptyDataError:
        return False

    # Rename columns for easier reference
    df.columns = [
        'score', 'perc_div', 'perc_del', 'perc_ins', 'query_name',
        'query_start', 'query_end', 'query_left', 'strand',
        'repeat_name', 'repeat_class', 'repeat_start',
        'repeat_end', 'repeat_left', 'ID'
    ]

    # Filter rows based on query_iden
    df = df[df['perc_div'] <= (100 - min_iden)]

    # Calculate coverage length and add to a column
    df['cov_len'] = abs(df['query_start'] - df['query_end'])

    # Select dataframe columns
    df_filter = df[["query_name", "repeat_name", "repeat_class", "cov_len"]]

    # Group by the columns and sum 'cov_len'
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
                seq_name = line.split('#')[0][1:].strip()  # Remove '>' and split at '#', then take the first part

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
        # Extract the part before #
        consensus_name = filename.split('#')[0]
        for key, value in reclassified_dict.items():
            if consensus_name == key:
                # Extract the part after # and before the first .
                te_type = filename.split('#')[1].split('.')[0].replace('__', '/')

                # If it doesn't match the value in the dictionary, rename the file
                if te_type != value:
                    new_te_type = value.replace('/', '__')
                    new_filename = filename.replace(te_type.replace('/', '__'), new_te_type)
                    old_file_path = os.path.join(directory, filename)
                    new_file_path = os.path.join(directory, new_filename)

                    # If seq_name is True, rename the sequence header in the FASTA file
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


def remove_files_with_start_pattern(input_dir, start_pattern):
    # Remove files and folder start with give pattern
    for filename in os.listdir(input_dir):
        if filename.startswith(start_pattern):
            file_path = os.path.join(input_dir, filename)
            if os.path.isfile(file_path):  # This check ensures you're only removing files, not directories
                os.remove(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)


# Define a function to handle sequence skipping and removal of files
def handle_sequence_skipped(seq_obj, progress_file, keep_intermediate, MSA_dir, classification_dir):

    try:
        seq_name = seq_obj.get_seq_name()
        seq_obj.update_status("skipped", progress_file)
        if not keep_intermediate:
            remove_files_with_start_pattern(MSA_dir, seq_name)
            remove_files_with_start_pattern(classification_dir, seq_name)
    except Exception as e:
        raise Exception(
            f"An error occurred while handling skipped sequence {seq_name}: {e}")


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


# if the seq_obj is low copy, append to consensus_file or final_unknown_con_file file
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

        # Write all consensus sequence to final_cons_file.
    with open(consensus_file, "a") as f:
        f.write(">" + seq_name + "#" + te_type + "\n" + sequence + "\n")

    low_copy_single_fasta_file = os.path.join(proof_dir, f"{seq_name}#{te_type_modified}.fa")
    low_copy_te_aid_pdf_file = os.path.join(proof_dir, f"{seq_name}#{te_type_modified}.pdf")

    shutil.copy(input_fasta, low_copy_single_fasta_file)
    shutil.copy(te_aid_pdf, low_copy_te_aid_pdf_file)



# Classify single fasta
def classify_single(consensus_fasta):
    """
    Run RepeatClassifier with the provided parameters.
    """

    # Store the current working directory
    original_dir = os.getcwd()

    # Change the working directory to the directory of the consensus_fasta
    os.chdir(os.path.dirname(consensus_fasta))

    # Define RepeatClassifier command, the output file will store at the same directory of the consensus_fasta
    command = ["RepeatClassifier", "-consensi", consensus_fasta]

    # Run RepeatClassifier using subprocess
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Change the working directory back to the original directory
    os.chdir(original_dir)

    if result.returncode != 0:
        click.echo(f"RepeatClassifier error for {os.path.basename(consensus_fasta)}\n{result.stderr.decode('utf-8')}")
        return False

    classified_file = f'{consensus_fasta}.classified'

    # Get the first record of classified file
    record = next(SeqIO.parse(classified_file, "fasta"))

    # seq_name = record.id.split("#")[0]
    seq_TE_type = record.id.split("#")[-1]

    return seq_TE_type






