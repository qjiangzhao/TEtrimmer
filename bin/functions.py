import shutil
import subprocess
import os
import click
import shutil
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment, AlignInfo
import pandas as pd
import pandas.errors
from seqclass import SeqObject
import numpy as np
from PyPDF2 import PdfMerger, PdfFileReader, PdfFileWriter


def prcyan(text):
    click.echo(click.style(text, fg='cyan'))


def prgre(text):
    click.echo(click.style(text, fg='green'))


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


def check_database(genome_file, search_type="blast"):
    """
    Checks if the blast database and genome length file exist.
    If they don't exist, create them.

    :param genome_file: str, path to genome file (contain genome name)
    """
    if search_type == "blast":
        blast_database_file = genome_file + ".nin"
        if not os.path.isfile(blast_database_file):
            print(f"\nBlast database doesn't exist. Running makeblastdb")

            try:
                makeblastdb_cmd = f"makeblastdb -in {genome_file} -dbtype nucl -out {genome_file} "
                subprocess.run(makeblastdb_cmd, shell=True, check=True, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, text=True)
            except FileNotFoundError:
                prcyan("'makeblastdb' command not found. Please ensure 'makeblastdb' is correctly installed.")
                raise Exception

            except subprocess.CalledProcessError as e:
                prcyan(f"makeblastdb failed with exit code {e.returncode}")
                prcyan(e.stderr)
                return False
    elif search_type == "mmseqs":
        mmseqs_database_dir = genome_file + "_db"
        if not os.path.isdir(mmseqs_database_dir):
            print(f"\nMMseqs2 database doesn't exist. Running MMseqs2 database creation...\n")

            try:
                mmseqs_createdb_cmd = f"mmseqs createdb {genome_file} {mmseqs_database_dir}"
                subprocess.run(mmseqs_createdb_cmd, shell=True, check=True,stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, text=True)

            except subprocess.CalledProcessError as e:
                prcyan(f"mmseqs failed with exit code {e.returncode}")
                prcyan(e.stderr)

            # Check if MMseqs2 database files have been created
            if os.path.exists(f"{mmseqs_database_dir}.dbtype"):
                print(f"MMseqs2 database created successfully: {mmseqs_database_dir}")

                # Create index for the MMseqs2 database
                print("Creating index for MMseqs2 database...")
                mmseqs_createindex_cmd = f"mmseqs createindex {mmseqs_database_dir} {mmseqs_database_dir}_tmp " \
                                         f"--search-type 3"
                result = subprocess.run(mmseqs_createindex_cmd, shell=True, check=True, stderr=subprocess.PIPE)
                error_output = result.stderr.decode("utf-8")
                if error_output:
                    print(f"Error creating MMseqs2 index: {error_output}\n")
                else:
                    print("MMseqs2 index created successfully.")
            else:
                print(f"Error: MMseqs2 database files not found for {genome_file}, you can build it by yourself to"
                      f"aovid this error")

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

        try:
            subprocess.run(faidx_cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        except FileNotFoundError:
            prcyan("'samtools' command not found. Please ensure 'samtools' is correctly installed.")
            return False

        except subprocess.CalledProcessError as e:
            prcyan(f"\nsamtool faidx failed with error code {e.returncode}")
            prcyan(e.stderr)
            prgre("Please check if your samtools works well. Or you can build the genome index file by yourself\n"
                  "samtools faidx <your_genome>")
            return False
    return True


def separate_sequences(input_file, output_dir, continue_analysis=False):
    """
    separates input file into single fasta file and creates object for each input sequence
    """
    os.makedirs(output_dir, exist_ok=True)
    seq_list = []

    if not continue_analysis:
        print(
            "TE Trimmer is modifying sequence names; All '/', '-', ':', '...', '|' and empty space before '#' will "
            "be converted to '_'\n"
            "You can find the original and modified name relationship from 'Name_mapping.txt' file under "
            "the output directory.\n ")
        # Initialize the name mapping file
        name_mapping_file = os.path.join(os.path.dirname(output_dir), "Name_mapping.txt")

        detected_pound = False
        with open(input_file, 'r') as fasta_file, open(name_mapping_file, 'w') as mapping_file:
            # Write header to the mapping file
            mapping_file.write("original_input_seq_name\tTETrimmer_modified_seq_name\n")
            id_list = []
            # Have to add 'fasta' at the end, this pattern will be used for file deletion
            for record in SeqIO.parse(fasta_file, 'fasta'):
                # Check if # is in the seq.id. If # is present, the string before # is the seq_name, and the string
                # after # is the seq_TE_type
                if len(record.id.split("#")) > 1:
                    detected_pound = True
                    sanitized_id = record.id.split("#")[0].replace('/', '_').replace(' ', '_').\
                        replace('-', '_').replace(':', '_').replace('...', '_').replace('|', '_')
                    te_type = record.id.split("#")[0]

                    # Normally SeqIO.parse only takes content before " " as record.id. Separate with " " to make
                    # the code stronger
                    te_type = te_type.split(" ")[0]

                else:
                    sanitized_id = record.id.replace('/', '_').replace(' ', '_').replace('-', '_')\
                        .replace(':', '_').replace('...', '_')
                    te_type = "Unknown"
                    # modify header to add #Unknown 
                    record.id = f"{record.id}#{te_type}"
                    record.description = record.id

                # double check if sanitized_id is unique. If not, modify sanitized_id
                if sanitized_id not in id_list:
                    id_list.append(sanitized_id)
                else:
                    # print(f"duplicated seq_name {sanitized_id} during separate_sequences.")
                    id_list.append(sanitized_id)
                    count = id_list.count(sanitized_id)
                    sanitized_id = f"{sanitized_id}_n{count}"

                # Write original and modified names to the mapping file
                mapping_file.write(f"{record.id}\t{sanitized_id}\n")

                # Define output file name
                output_filename = os.path.join(output_dir, f"{sanitized_id}.fasta")
                seq_obj = SeqObject(str(sanitized_id), str(output_filename), len(record.seq), te_type)

                # Store all input file information (object manner) to seq_list
                seq_list.append(seq_obj)

                # Convert sequence name to sanitized_id
                record.id = sanitized_id

                # Write single fasta file with sanitized name
                with open(output_filename, 'w') as output_file:
                    SeqIO.write(record, output_file, 'fasta')

            if detected_pound:
                print("TE Trimmer detects # in your input fasta sequence. The string before # is denoted as the "
                      "seq_name, and the string after # is denoted as the TE type\n")

        print("\nFinish to generate single sequence files.\n")

    elif continue_analysis:
        # When continue_analysis is true, generate seq_list based on single fasta files
        for filename in os.listdir(output_dir):
            file = os.path.join(output_dir, filename)
            with open(file, 'r') as fasta_file:
                for record in SeqIO.parse(fasta_file, 'fasta'):
                    # Get sanitized_id from single fasta file name
                    sanitized_id = os.path.splitext(filename)[0]

                    # Note: single fasta file name is different with record.id
                    te_type = record.id.split("#")[-1]
                    seq_obj = SeqObject(str(sanitized_id), str(file), len(record.seq), te_type)
                    seq_list.append(seq_obj)
        print("\nFinish to read in single sequence files generated from previous analysis.\n")
    return seq_list


def blast(seq_file, genome_file, output_dir, min_length=150, search_type="blast", task="blastn", seq_obj=None):
    """
    Runs BLAST with specified task type and saves the results in a bed file.

    :param seq_file: str, path to input fasta file
    :param genome_file: str, path to genome file
    :param output_dir: str, prefix for output files
    :param min_length: num, default 150, minimum alignment length
    :param task: str, default "blastn", BLAST task type ("blastn", "dc-megablast", etc.)
    :param seq_obj: object, optional, sequence object to update with blast hits

    :return: tuple, bed output file name and blast hits count
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
            prcyan("'blastn' command not found. Please ensure 'blastn' is correctly installed.")
            raise Exception

        except subprocess.CalledProcessError as e:
            prcyan(f"\nBLAST is failed for {input_file_n} with error code {e.returncode}")
            prcyan(e.stderr)
            raise Exception

    elif search_type == "mmseqs":
        # MMseqs2 search
        mmseqs_cmd = (f"mmseqs easy-search {seq_file} {genome_file}_db {blast_out_file} {blast_out_file}_tmp --search-type 3 "
                      f"--min-seq-id 0.6 --format-output \"query,target,pident,alnlen,mismatch,qstart,qend,tstart,tend,evalue,qcov\" "
                      f"--cov-mode 4 -c 0.5 --e-profile 1e-40 --threads 1 --min-aln-len {min_length}")
        result = subprocess.run(mmseqs_cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        error_output = result.stderr.decode("utf-8")

        if error_output:
            print(f"Error during MMseqs2: {error_output}")

    # Check the number of BLAST hits
    with open(blast_out_file) as blast_file:
        for _ in blast_file:
            blast_hits_count += 1

    if blast_hits_count > 0:
        # Define bed outfile
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
            prcyan(f"\nBLAST to bed file conversion failed for {input_file_n} with error code {e.returncode}")
            prcyan(e.stderr)
            raise Exception

        # Update BLAST hit number to sequence object when seq_obj is supplied
        if seq_obj is not None:
            seq_obj.update_blast_hit_n(blast_hits_count)

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


def extract_fasta(input_file, genome_file, output_dir, left_ex, right_ex, nameonly=False):
    """
    Extracts fasta sequence from the reference genome using bedtools.

    :param genome_file: str, path to genome file
    :param output_dir: str, prefix for output files
    :param left_ex: int, number of bases to extend the start position of each feature
    :param right_ex: int, number of bases to extend the end position of each feature
    :param nameonly: boolean Default: False: only use bed file name filed for the fasta header
    :return: fasta file absolute path derived from bed file
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
        prcyan("'bedtools slop' command not found. Please ensure 'bedtools' is correctly installed.")
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
        prcyan("'bedtools getfasta' command not found. Please ensure 'bedtools' is correctly installed.")
        raise Exception

    except subprocess.CalledProcessError as e:
        prcyan(f"\nbedtools getfasta failed for {input_file_n} with error code {e.returncode}")
        prcyan(e.stderr)
        raise Exception

    # Use awk to remove letters that aren't A G C T a g c t
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
    Aligns fasta sequences using MAFFT

    :param input_file: str, input file path
    :param output_dir: str, output file directory
    :return: multiple sequence alignment file absolute path
    """
    input_file_n = os.path.basename(input_file)
    fasta_out_flank_mafft_file = os.path.join(output_dir, f"{os.path.basename(input_file)}_aln.fa")

    # Read the top 5 sequences and determine if any are longer than 7000
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

    # If any of the top 5 sequences are longer than 7000, add --memsave to salve memory
    if long_sequences:
        mafft_cmd.insert(1, "--memsave")  # Insert after 'mafft'

    try:
        # Execute the command
        result = subprocess.run(mafft_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    except FileNotFoundError:
        prcyan("'mafft' command not found. Please ensure 'mafft' is correctly installed.")
        raise Exception

    except subprocess.CalledProcessError as e:
        error_message = e.stderr
        prcyan(f"MAFFT failed for {input_file_n} with error code {e.returncode}\n")
        prcyan(e.stderr)

        if "Killed" in error_message and "disttbfast" in error_message and "memopt" in error_message:
            prgre(f"It seems not enough 'RAM' was given for Mafft multiple sequence alignment. Please assign more "
                  f"RAM or lower the thread number to solve this problem.\n")
        raise Exception

    # Write the output to the file
    with open(fasta_out_flank_mafft_file, 'w') as f:
        f.write(result.stdout)

    return fasta_out_flank_mafft_file


def muscle_align(input_file, output_dir, ite_times=4):
    output_file = os.path.join(output_dir, f"{os.path.basename(input_file)}_maln.fa")
    muscle_cmd = ["muscle", "-maxiters", str(ite_times), "-in", input_file, "-out", output_file]

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

        # Calculate the proportion of each nucleotide
        total = sum(counts.values())

        # Check if total is zero
        if total == 0:
            max_props.append(0)
        else:
            proportions = {nucleotide: count / total for nucleotide, count in counts.items()}
            max_props.append(max(proportions.values()))  # Append the maximum proportion

    return np.array(max_props)  # Convert the list to a NumPy array and return it


# Function help to define the threshold of crop end by similarity
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


def generate_hmm_from_msa(input_msa_file, output_hmm_file, error_file):
    """
    Generate HMM profile using hmmbuild from a multiple sequence alignment file.

    :param input_msa_file: Path to the multiple sequence alignment file (in FASTA or Stockholm format).
    :param output_hmm_file: Path to the output HMM file.
    """

    # Construct the command as a list
    cmd = ["hmmbuild", "--dna", output_hmm_file, input_msa_file]

    try:
        # Execute the command
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    except FileNotFoundError:
        prcyan("'hmmbuild' command not found. Please ensure 'hmmbuild' is correctly installed.")
        prgre("This hmm error won't affect the final TE consensus library.")
        pass

    except subprocess.CalledProcessError as e:
        with open(error_file, 'a') as f:
            f.write(f"\nhmm file generation failed for {os.path.basename(output_hmm_file)} with error code {e.returncode}")
            f.write('\n' + e.stderr)
        prcyan(f"\nhmm file generation failed for {os.path.basename(output_hmm_file)} with error code {e.returncode}")
        prgre("\nThis hmm error won't affect the final TE consensus library, you can ignore it."
              "\nFor traceback text, please refer to 'error_file.txt' under 'Multiple_sequence_alignment' folder\n")
        pass


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


def remove_short_seq(input_file, output_file, threshold=0.1):
    alignment = AlignIO.read(input_file, "fasta")

    # Get MSA length
    len = alignment.get_alignment_length()

    # Initialize an empty list to store the filtered sequences
    filtered_alignment = []

    # Loop each sequence in the alignment
    for record in alignment:
        # Get the length of the sequence without gaps
        seq_length = len(record.seq.ungap("-"))

        # Check if the sequence length is greater than or equal the length threshold
        if seq_length >= len * threshold:
            filtered_alignment.append(record)

    # Write the filtered list record to a file
    AlignIO.write(filtered_alignment, output_file)


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


def remove_gaps_with_similarity_check(input_file, output_dir, gap_threshold=0.8, simi_check_gap_thre=0.4,
                                      similarity_threshold=0.7, min_nucleotide=5, return_map=False):
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


# Select MSA columns according to the given star and end position
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
    """
    -l	length of throw_away_sequences, default 10
    -d	length of description in .clstr file, default 20 if set to 0, it takes the fasta defline and
    stops at first space
    -s	length difference cutoff, default 0.0 if set to 0.9, the shorter sequences need to be at least 90% length
    of the representative of the cluster
    -aL	alignment coverage for the longer sequence, default 0.0 if set to 0.9, the alignment must
    covers 90% of the sequence
    -aS	alignment coverage for the shorter sequence, default 0.0 if set to 0.9, the alignment must
    covers 90% of the sequence
    Note: -s only consider length, but -aL and -aS consider alignment
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
        "-l", "50",
        "-d", "0",
        "-s", str(s)
    ]

    try:
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    except FileNotFoundError:
        prcyan("'\ncd-hit-est' command not found. Please ensure 'cd-hit-est' is correctly installed.\n")
        raise Exception

    except subprocess.CalledProcessError as e:
        prcyan(f"\ncd-hit-est failed for {os.path.basename(input_file)} with error code {e.returncode}")
        prcyan(e.stderr)
        raise Exception


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
        prcyan("'RepeatMasker' command not found. Please ensure 'RepeatMasker' is correctly installed.")
        raise Exception

    except subprocess.CalledProcessError as e:
        if classify:
            prcyan(f"\nRepeatMasker failed during final classification step with error code {e.returncode}")
            prcyan(e.stderr)
            prgre("This will not affect the final result but only the classification of TE might not be affect")
            raise Exception
        else:
            prcyan(f"\nRepeatMasker failed during final whole genome TE annotation with error code {e.returncode}")
            prcyan(e.stderr)
            prgre("This won't affect the final TE consensus library at all. You can do the final genome TE annotation"
                  " by yourself with RepeatMasker")
            raise Exception


def repeatmasker_output_classify(repeatmasker_out, progress_file, min_iden=70, min_len=80, min_cov=0.8):
    # Read RepeatMasker out file into a DataFrame
    # The regex '\s+' matches one or more whitespace characters
    # error_bad_lines=False to skip errors
    try:
        df = pd.read_csv(repeatmasker_out, delim_whitespace=True, header=None, skiprows=3, usecols=range(15))
    except pandas.errors.EmptyDataError:
        return False
    except pd.errors.ParserError:
        df = pd.read_csv(repeatmasker_out, delim_whitespace=True, header=None, skiprows=3,
                         error_bad_lines=False, usecols=range(15))

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


def remove_files_with_start_pattern(input_dir, start_pattern, if_seq_name=True):
    # Remove files and folder start with give pattern
    # Add .fasta to the end of start_pattern
    if if_seq_name:
        start_pattern = f"{start_pattern}.fasta"
    for filename in os.listdir(input_dir):
        if filename.startswith(start_pattern):
            file_path = os.path.join(input_dir, filename)
            if os.path.isfile(file_path):  # This check ensures you're only removing files, not directories
                os.remove(file_path)
            elif os.path.isdir(file_path):
                try:
                    shutil.rmtree(file_path)
                except Exception:
                    # File deletion doesn't affect the final consensus sequence. Skip when error happens
                    pass


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
            remove_files_with_start_pattern(MSA_dir, seq_name)
            remove_files_with_start_pattern(classification_dir, seq_name)
    except Exception as e:
        click.echo(f"\nAn error occurred while handling low copy sequence {seq_name}:\n {e}\n")


def handle_sequence_skipped(seq_obj, progress_file, debug, MSA_dir, classification_dir, plot_skip=True,
                            te_aid_plot=None, orf_plot=None, skip_proof_dir=None):
    seq_name = seq_obj.get_seq_name()
    te_type = seq_obj.get_old_TE_type()
    te_type_modified = te_type.replace("/", "__")
    try:
        seq_obj.update_status("skipped", progress_file)
        if plot_skip and (te_aid_plot is not None or orf_plot is not None) and skip_proof_dir is not None:
            merge_pdfs(skip_proof_dir, f"{seq_name}#{te_type_modified}", te_aid_plot, orf_plot)
        if not debug:
            remove_files_with_start_pattern(MSA_dir, seq_name)
            remove_files_with_start_pattern(classification_dir, seq_name)
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


# if the seq_obj is low copy, append to consensus_file and final_unknown_con_file used for classification
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
    #low_copy_te_aid_pdf_file = os.path.join(proof_dir, f"{seq_name}#{te_type_modified}_TE_Aid.pdf")

    shutil.copy(input_fasta, low_copy_single_fasta_file)

    # Sometimes TE Aid can't be plotted properly because the input sequence quality. Only move plot when it is existed.
    #if os.path.exists(te_aid_pdf) and os.path.getsize(te_aid_pdf) > 0:
        #shutil.copy(te_aid_pdf, low_copy_te_aid_pdf_file)


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

    try:
        # Run RepeatClassifier using subprocess
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    except FileNotFoundError:
        prcyan("'RepeatClassifier' command not found. Please ensure 'RepeatModeler' is correctly installed.")
        return False

    except subprocess.CalledProcessError as e:
        prcyan(f"RepeatClassifier error for {os.path.basename(consensus_fasta)} with error code {e.returncode}")
        prcyan(e.stderr)
        prgre("This only affect classification but not consensus sequence. "
              "You can run 'RepeatClassifier -consensi <your_consensus_file>' to test")
        # Change the working directory back to the original directory
        os.chdir(original_dir)
        return False

    # Change the working directory back to the original directory
    os.chdir(original_dir)

    classified_file = f'{consensus_fasta}.classified'

    # Get the first record of classified file
    record = next(SeqIO.parse(classified_file, "fasta"))

    # seq_name = record.id.split("#")[0]
    seq_TE_type = record.id.split("#")[-1]

    return seq_TE_type


def check_terminal_repeat(input_file, output_dir, teaid_blast_out=None, TIR_adj=2000, LTR_adj=3000):
    """
    output_dir: used to store self blast database
    """
    # Read input file and get sequence length
    record = SeqIO.read(input_file, "fasta")
    record_len = len(record.seq)
    found_match = False
    LTR_boundary = None
    TIR_boundary = None

    # teaid output is not given or the given file doesn't exist or is empty, do self blast
    if teaid_blast_out is None or not file_exists_and_not_empty(teaid_blast_out):
        os.makedirs(output_dir, exist_ok=True)
        database_file = os.path.join(output_dir, "Tem_blast_database")
        makeblastdb_cmd = f"makeblastdb -in {input_file} -dbtype nucl -out {database_file}"

        # If error encountered, return directly
        try:
            subprocess.run(makeblastdb_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError:
            return None, None, False

        # Proof check if database is built, otherwise return directly
        if not file_exists_and_not_empty(f"{database_file}.nhr") or not \
                file_exists_and_not_empty(f"{database_file}.nin") or not file_exists_and_not_empty(f"{database_file}.nsq"):
            return None, None, False

        blast_cmd = f"blastn -query {input_file} -db {database_file} " \
                    f"-outfmt \"6 qseqid qstart qend sstart send \" " \
                    f"-evalue 0.05"  # Set higher e-value for self-blast

        try:
            result = subprocess.run(blast_cmd, shell=True, check=True, text=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            blast_out = result.stdout

            # If self blast result is empty, return directly
            if blast_out.strip() == "":
                return None, None, False
        except subprocess.CalledProcessError as e:
            prcyan(f"\nTerminal repeat detection failed for {os.path.basename(input_file)} with error code {e.returncode}")
            prcyan('\n' + e.stderr + '\n')
            raise Exception

        # Define blast out file
        blast_out_file = os.path.join(output_dir, "tem_blast_out.txt")

        with open(blast_out_file, 'w') as f:
            # Give a header to blast result
            f.write("qseqid\tqstart\tqend\tsstart\tsend\n")
            f.write(blast_out)
        df = pd.read_csv(blast_out_file, sep="\t", header=None, skiprows=1)
    else:
        blast_out_file = teaid_blast_out

        # TEaid self blast output default format is separated by white space
        df = pd.read_csv(blast_out_file, sep='\s+', header=None, skiprows=1)

    # Return None directly when the self blast result is empty, this could happen when too many ambiguous letters
    # are there like N.
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

        # Check if the terminal repeat spans the most part of the query sequence. Because the query is after extension,
        # assuming the maximum redundant extension for left and right side are both 2000.
        if abs(LTR_largest[4] - LTR_largest[1]) >= (record_len - LTR_adj):
            # Because blast use index start from 1, modify the start position
            LTR_boundary = [LTR_largest[1] - 1, LTR_largest[4]]
            found_match = "LTR"
        else:
            LTR_boundary = None

        # Return directly when LTR is found
        return LTR_boundary, TIR_boundary, found_match

    df_TIR = df[(df.iloc[:, 2] - df.iloc[:, 1] >= 50) &
            (df.iloc[:, 1] != df.iloc[:, 3]) &
            (df.iloc[:, 2] != df.iloc[:, 4]) &
            (df.iloc[:, 3] > df.iloc[:, 4]) &
            (df.iloc[:, 4] > df.iloc[:, 1])].copy()
    
    if not df_TIR.empty:
        df_TIR["5"] = df.iloc[:, 3] - df.iloc[:, 1]
        df_TIR.reset_index(drop=True, inplace=True)

        # Same like LTR check the terminal repeat spanning region
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
    Merge PDF files to one single file. the file path order in *pdfs is the file order in the final single file
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

    # Define after TETrimmer treatment file name
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
        "-gxtitle", "After TETrimmer treatment",
        "-gytitle", "Without TETrimmer treatment",
        "-graph", "ps",
        "-goutfile", sequence1
    ]

    try:
        subprocess.run(dotmatcher_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    except FileNotFoundError:
        prcyan("'\ndotmatcher' command not found. Please ensure 'emboss' is correctly installed.")
        prgre("\ndotmatcher won't affect the final consensus sequence. You can choose to ignore this error\n")
        raise Exception

    except subprocess.CalledProcessError as e:
        prcyan(f"\ndotmatcher failed for {n_after_tetrimmer} with error code {e.returncode}")
        prgre("\ndotmatcher won't affect the final consensus sequence. You can choose to ignore this error")
        raise Exception

    # Define command to convert ps to pdf
    ps2pdf_command = [
        "ps2pdf",
        f"{sequence1}.ps",
        pdf_out
    ]

    try:
        subprocess.run(ps2pdf_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    except FileNotFoundError:
        prcyan("'\nps2pdf' command not found. Please install it by 'sudo apt-get install ghostscript'")
        prgre("ps2pdf won't affect the final consensus file. You can choose to ignore it.")
        raise Exception

    except subprocess.CalledProcessError as e:
        prcyan(f"\nps2pdf failed for {n_after_tetrimmer} with error code {e.returncode}")
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


