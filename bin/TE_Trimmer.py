# Standard library imports
import os
import traceback
from datetime import timedelta, datetime
import multiprocessing as mp
import click
import concurrent.futures
import shutil
import json
import csv

# Local imports
from Function_separate_fasta import separate_sequences
from Function_blast_extension_mafft import blast, remove_files_with_start_pattern, change_permissions_recursive, \
    repeatmasker, check_database, check_bed_uniqueness, extract_fasta, cd_hit_est, handle_sequence_skipped, \
    repeatmasker_output_classify, update_cons_file, update_low_copy_cons_file
from Class_bed_filter import BEDFile
from Function_def_boundary_and_crop import find_boundary_and_crop
from Class_TE_aid import check_self_alignment
from Function_clean_and_clauster_MSA import clean_and_cluster_MSA
from Class_orf_domain_prediction import prepare_pfam_database


# Define a function to check progress file, which will be used for continue analysis
def check_progress_file(progress_file_path):

    skipped_count = 0
    low_copy_count = 0
    local_completed_sequences = []

    # Gather all 'input_name' values using DictReader
    for row in csv.DictReader(open(progress_file_path, 'r')):
        local_completed_sequences.append(row['input_name'])

        # Calculate skipped and low copy element number
        if row['status'] == 'skipped':
            skipped_count += 1
        if row['low_copy'] == 'True':
            low_copy_count += 1

    # Remove duplicates for processed file IDs
    local_completed_sequences = list(set(local_completed_sequences))

    return local_completed_sequences, skipped_count, low_copy_count

# Print iterations progress
def printProgressBar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█', printEnd="\r"):
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
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}|{iteration}/{total}={percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()

#####################################################################################################
# Code block: Define analyze_sequence function
#####################################################################################################

def analyze_sequence_helper(params):
    return analyze_sequence(*params)


def analyze_sequence(seq_obj, single_file_dir, genome_file, MSA_dir, min_blast_len, min_seq_num, max_msa_lines,
                     top_mas_lines, max_cluster_num, cons_thr, ext_thr, ex_step, classification_dir,
                     max_extension, gap_thr, gap_nul_thr, crop_end_thr, crop_end_win, crop_end_gap_thr,
                     crop_end_gap_win, start_patterns, end_patterns, output_dir, pfam_dir, mini_orf,
                     single_fasta_n, hmm, check_extension_win, keep_intermediate, progress_file,
                     classify_unknown, classify_all, final_con_file, final_unknown_con_file, final_classified_con_file):

    #####################################################################################################
    # Code block: Elongate query sequence when it is too short
    #####################################################################################################
    try:
        # Get query fasta file path
        seq_name = seq_obj.get_seq_name()
        seq_type = seq_obj.get_old_TE_type()
        seq_file = seq_obj.get_input_fasta()

        # Due to DNA element will be much shorter than LTR and LINE elements, set different parameters
        if "DNA" in seq_type:
            ex_step = 500
            max_extension = 7000
            min_blast_len = 150
            crop_end_gap_win = 100
            check_extension_win = 50

        # The average length of SINE element is around 500 bp, give different default parameters
        if "SINE" in seq_type:
            ex_step = 200
            max_extension = 1400
            min_blast_len = 80
            crop_end_gap_win = 50
            check_extension_win = 50

        if "Helitron" in seq_type:
            ex_step = 500
            max_extension = 7000
            min_blast_len = 150
            crop_end_gap_win = 100
            check_extension_win = 50

        if "MITE" in seq_type:
            ex_step = 100
            max_extension = 500
            min_blast_len = 50
            crop_end_gap_win = 40
            check_extension_win = 50

        # run blast for each single fasta file and return a bed file absolute path
        bed_out_file_dup, blast_hits_count, blast_out_file = blast(seq_obj, seq_file, genome_file, MSA_dir,
                                                                             min_length=min_blast_len)

    except Exception as e:
        click.echo(f"Error while running blast for sequence: {seq_name}. Main Error: {str(e)}")
        return

    try:

        # Check if blast hit number is equal 0, then skip this sequence
        if blast_hits_count == 0:

            click.echo(f"\n{seq_name} in main is skipped due to blast hit number is 0\n")
            handle_sequence_skipped(seq_obj, progress_file, keep_intermediate, MSA_dir, classification_dir)

            return

        # Check if blast hit number is smaller than "min_seq_num", not include "min_seq_num"
        elif blast_hits_count != 0 and blast_hits_count < min_seq_num:

            check_low_copy, blast_full_length_n, found_match = check_self_alignment(seq_obj, seq_file, MSA_dir,
                                                                                    genome_file, blast_hits_count,
                                                                                    blast_out_file)
            if check_low_copy is True:

                # Update terminal repeat and blast full length number
                seq_obj.set_old_terminal_repeat(found_match)
                seq_obj.set_old_blast_full_n(blast_full_length_n)
                seq_obj.update_status("processed", progress_file)
                update_low_copy_cons_file(seq_obj, final_con_file, final_unknown_con_file, final_classified_con_file)
            else:
                click.echo(f"\n{seq_name} in main is skipped due to check_low_copy is {check_low_copy}\n")

                # handle_sequence_skipped will update skipped status
                handle_sequence_skipped(seq_obj, progress_file, keep_intermediate, MSA_dir, classification_dir)

            return  # when blast hit number is smaller than 10, code will execute next fasta file

    except Exception as e:
        click.echo(f"Error while checking uniqueness for sequence: {seq_name}. Error: {str(e)}")
        return

    try:
        # remove duplicated lines
        bed_out_file = check_bed_uniqueness(MSA_dir, bed_out_file_dup)

        # test bed_out_file line number and extract top longest lines
        # return bed_out_filter_file absolute path
        bed_out_filter = BEDFile(bed_out_file)
        # for process_lines() function. threshold represent the maximum number to keep for MSA
        # top_longest_lines_count means the number of sequences with top length
        # for example if threshold = 100, top_longest_lines_count = 50, then 50 sequences will be
        # randomly chose from the rest of sequences
        # top_mas_lines has to be equal or smaller than max_mas_lines
        bed_out_filter_file = bed_out_filter.process_lines(MSA_dir, threshold=max_msa_lines,
                                                           top_longest_lines_count=top_mas_lines)

        # extract fast from bed_out_filter_file
        # return fasta_out_flank_file absolute path
        # because have to group MSA the first round extend for left and right side are both 0
        fasta_out_flank_file, bed_out_flank_file = extract_fasta(
            bed_out_filter_file, genome_file, MSA_dir, left_ex=0, right_ex=0)

        # Return False when cluster number is 0. Return True when divergent column number is smaller than 100
        # Otherwise it will return the subset bed and alignment file
        cluster_MSA_result = clean_and_cluster_MSA(fasta_out_flank_file, bed_out_filter_file, MSA_dir,
                                                   clean_column_threshold=0.08,
                                                   min_length_num=min_seq_num, cluster_num=max_cluster_num,
                                                   cluster_col_thr=100)
    except Exception as e:
        click.echo(
            f"Error during processing lines, extracting fasta, or clustering MSA for sequence: {seq_name}. Error: {str(e)}")
        traceback.print_exc()
        return

    try:
        # cluster false means no cluster, TE Trimmer will skip this sequence.
        if cluster_MSA_result is False:

            check_low_copy, blast_full_length_n, found_match = check_self_alignment(seq_obj, seq_file, MSA_dir,
                                                                                    genome_file, blast_hits_count,
                                                                                    blast_out_file)
            if check_low_copy is True:

                # Update terminal repeat and blast full length number
                seq_obj.set_old_terminal_repeat(found_match)
                seq_obj.set_old_blast_full_n(blast_full_length_n)
                seq_obj.update_status("processed", progress_file)
                update_low_copy_cons_file(seq_obj, final_con_file, final_unknown_con_file, final_classified_con_file)

            else:
                click.echo(
                    f"\n{seq_name} is skipped due to cluster_MSA_result sequence number in each cluster is smaller "
                    f"than {min_seq_num} and check_low_copy is {check_low_copy}\n")
                handle_sequence_skipped(seq_obj, progress_file, keep_intermediate, MSA_dir, classification_dir)

            return

        # Cluster True means not necessary to cluster MSA, perform find_boundary_and_crop directly
        elif cluster_MSA_result is True:
            find_boundary_result = find_boundary_and_crop(
                    bed_out_filter_file, genome_file, MSA_dir, pfam_dir, seq_obj, hmm,
                    classify_all, classify_unknown,
                    cons_threshold=cons_thr, ext_threshold=ext_thr,
                    ex_step_size=ex_step, max_extension=max_extension,
                    gap_threshold=gap_thr, gap_nul_thr=gap_nul_thr,
                    crop_end_thr=crop_end_thr, crop_end_win=crop_end_win,
                    crop_end_gap_thr=crop_end_gap_thr, crop_end_gap_win=crop_end_gap_win,
                    start_patterns=start_patterns, end_patterns=end_patterns, mini_orf=mini_orf,
                    define_boundary_win=check_extension_win
            )
            if find_boundary_result == "Short_sequence":
                click.echo(f"\n{seq_name} is skipped due to too short length\n")

                handle_sequence_skipped(seq_obj, progress_file, keep_intermediate, MSA_dir, classification_dir)

                return False

            elif not find_boundary_result:  # This means the errors happen in the function
                return

        # else means need cluster
        else:
            cluster_bed_files_list = cluster_MSA_result

            # cluster_pattern_alignment_list have the same index with cluster_bed_files_list
            all_inner_skipped = True
            for i in range(len(cluster_bed_files_list)):

                # Based on the bed file list, extract fasta file
                inner_fasta_out_flank_file, inner_bed_out_flank_file = extract_fasta(
                    cluster_bed_files_list[i], genome_file, MSA_dir, left_ex=0, right_ex=0)

                inner_cluster_MSA_result = clean_and_cluster_MSA(inner_fasta_out_flank_file,
                                                                 cluster_bed_files_list[i], MSA_dir,
                                                                 clean_column_threshold=0.08,
                                                                 min_length_num=min_seq_num,
                                                                 cluster_num=max_cluster_num,
                                                                 cluster_col_thr=100)
                # inner_cluster_MSA_result is false means this cluster sequence number is too less
                if inner_cluster_MSA_result is False:
                    continue
                elif inner_cluster_MSA_result is True:  # Means don't need to cluster

                    inner_find_boundary_result = find_boundary_and_crop(
                            cluster_bed_files_list[i], genome_file, MSA_dir, pfam_dir, seq_obj,
                            hmm, classify_all, classify_unknown, cons_threshold=cons_thr, ext_threshold=ext_thr,
                            ex_step_size=ex_step, max_extension=max_extension,
                            gap_threshold=gap_thr, gap_nul_thr=gap_nul_thr,
                            crop_end_thr=crop_end_thr, crop_end_win=crop_end_win,
                            crop_end_gap_thr=crop_end_gap_thr, crop_end_gap_win=crop_end_gap_win,
                            start_patterns=start_patterns, end_patterns=end_patterns,
                            mini_orf=mini_orf, define_boundary_win=check_extension_win)

                    if inner_find_boundary_result == "Short_sequence":
                        continue
                    elif inner_find_boundary_result:
                        all_inner_skipped = False
                    elif not inner_find_boundary_result:  # This means the errors happen in the function
                        continue

                else:
                    inner_cluster_bed_files_list = inner_cluster_MSA_result

                    for j in range(len(inner_cluster_bed_files_list)):
                        inner_inner_find_boundary_result = find_boundary_and_crop(
                                inner_cluster_bed_files_list[j], genome_file, MSA_dir,
                                pfam_dir, seq_obj, hmm, classify_all, classify_unknown,
                                cons_threshold=cons_thr, ext_threshold=ext_thr,
                                ex_step_size=ex_step, max_extension=max_extension,
                                gap_threshold=gap_thr, gap_nul_thr=gap_nul_thr,
                                crop_end_thr=crop_end_thr, crop_end_win=crop_end_win,
                                crop_end_gap_thr=crop_end_gap_thr, crop_end_gap_win=crop_end_gap_win,
                                start_patterns=start_patterns, end_patterns=end_patterns,
                                mini_orf=mini_orf, define_boundary_win=check_extension_win)
                        if inner_inner_find_boundary_result == "Short_sequence":
                            continue
                        elif inner_inner_find_boundary_result:
                            all_inner_skipped = False
                        elif not inner_inner_find_boundary_result:  # This means the errors happen in the function
                            continue

            # Check the flag after the loop. If all inner clusters were skipped, write the progress file
            if all_inner_skipped:
                check_low_copy, blast_full_length_n, found_match = check_self_alignment(seq_obj, seq_file, MSA_dir,
                                                                                        genome_file, blast_hits_count,
                                                                                        blast_out_file)

                if check_low_copy is True:

                    # Update terminal repeat and blast full length number
                    seq_obj.set_old_terminal_repeat(found_match)
                    seq_obj.set_old_blast_full_n(blast_full_length_n)
                    seq_obj.update_status("processed", progress_file)
                    update_low_copy_cons_file(seq_obj, final_con_file, final_unknown_con_file, final_classified_con_file)
                else:
                    handle_sequence_skipped(seq_obj, progress_file, keep_intermediate, MSA_dir, classification_dir)
                    click.echo(
                        f"\n{seq_name} is skipped due to sequence number in second round each cluster is "
                        f"smaller than {min_seq_num} or the sequence is too short\n")

                return

    except Exception as e:
        click.echo(f"Error during boundary finding and cropping for sequence: {seq_name}. Error: {str(e)}\n")
        traceback.print_exc()
        return

    # After all processing is done, change status to process and write the name of the file to the progress file
    seq_obj.update_status("processed", progress_file)

    # If all this sequence is finished remove all files contain this name
    if not keep_intermediate:

        remove_files_with_start_pattern(MSA_dir, seq_name)
        remove_files_with_start_pattern(classification_dir, seq_name)

    # Read and count sequences from Finished_sequence_name.txt
    completed_sequence, skipped_count, low_copy_count = check_progress_file(progress_file)

    # Calculate the total count
    processed_count = len(completed_sequence)

    # Calculate sequences number that hasn't been processed by TE Trimmer
    rest_sequence = single_fasta_n - processed_count

    printProgressBar(processed_count, single_fasta_n, prefix='Progress:', suffix='Complete', length=50)

#####################################################################################################
# Code block: Import json species_config file and define the default parameters
#####################################################################################################


species_config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'species_config.json')

# Load the JSON configuration file
with open(species_config_path, "r") as config_file:
    species_config = json.load(config_file)


@click.command(context_settings=dict(max_content_width=120))
@click.option('--input_file', '-i', required=True, type=str,
              help='TE consensus fasta file. Use the output of RepeatModeler, EDTA, or REPET et al.')
@click.option('--genome_file', '-g', required=True, type=str,
              help='Genome file path.')
@click.option('--output_dir', '-o', default=os.getcwd(), type=str,
              help='Output directory. Default: current directory.')
@click.option('--species', '-s', required=True, default='fungi', type=click.Choice(species_config.keys()),
              help='Select the species for which you want to run TE Trimmer.')
@click.option('--continue_analysis', default=False, is_flag=True,
              help='Continue to analysis based on interrupted results.')
@click.option('--merge', default=False, is_flag=True,
              help='Merge input file to remove duplicate sequences.')
@click.option('--genome_anno', default=False, is_flag=True,
              help='Perform genome TE annotation based on TE Trimmer curated database at the end.')
@click.option('--hmm', default=False, is_flag=True,
              help='Generate HMM files for each consensus sequences.')
@click.option('--keep_intermediate', default=False, is_flag=True,
              help='Keep all raw files. WARNING: Many files will be produced.')
@click.option('--pfam_dir', default=None, type=str,
              help="Pfam database directory. Leave this option when you don't have Pfam database, "
                   "TE Trimmer will download automatically")
@click.option('--cons_thr', type=float,
              help='Threshold used for the final consensus sequence generation. Default: 0.8')
@click.option('--max_msa_lines', type=int,
              help='Set the maximum sequences number for multiple sequence alignment. Default: 100')
@click.option('--top_mas_lines', type=int,
              help='When the sequence number of multiple sequence alignment (MSA) is greater than '
                   '"max_mas_lines". It will order sequences by length and choose '
                   '"top_msa_lines" number of sequences. Then randomly choose the rest number of sequences '
                   'Default: 100')
@click.option('--min_seq_num', type=int,
              help='The minimum sequence number for each multiple sequence alignment. Default: 10')
@click.option('--min_blast_len', type=int,
              help='The minimum hit sequence length for blast. Default: 150')
@click.option('--max_cluster_num', type=int,
              help='The maximum cluster number for each multiple sequence alignment. Each multiple sequence alignment '
                   'can be divided into different clusters. TE Trimmer will sort cluster by sequence number and choose'
                   'the top --max_cluster_num of clusters for the further analysis. Default: 2')
@click.option('--ext_thr', type=float,
              help="threshold used for define the extension extent. The smaller number means it become easier "
                   "to have a final longer extension for each side of the sequence. Default: 0.7")
@click.option('--ex_step', type=int,
              help='Number of nucleotides will be added to the left or right side of multiple sequence alignment. '
                   'TE_Trimmer will iteratively add --ex_step number of nucleotide until finding the boundary. '
                   'Default: 1000')
@click.option('--max_extension', type=int,
              help='The maximum extension number for the right and left side. For example, if --ex_step is 1000, '
                   'it can only add seven times to the MSA left or right side. Default: 7000')
@click.option('--gap_thr', type=float,
              help='If columns have a larger gap proportion than --gap_thr and the most common nucleotide proportion '
                   'in this column is less than --gap_nul_thr, this column will be removed. Default: 0.4')
@click.option('--gap_nul_thr', type=float,
              help='Set nucleotide proportion to decide if remove this column. Coupled with --gap_nul_thr option. '
                   'Default: 0.7')
@click.option('--crop_end_win', type=int,
              help='Window size used for crop end process. Coupled with --crop_end_thr option. Default: 20')
@click.option('--crop_end_thr', type=int,
              help='Crop end function will convert each nucleotide in MSA into proportion number. This function will '
                   'check from the beginning and end of each sequence from MSA by iteratively choosing a slide window '
                   'and sum up the proportion numbers. It will stop until the summary of proportion is larger than '
                   '--crop_end_thr. The nucleotide do not match --crop_end_thr will be converted to -. The recommended '
                   'number is 0.8 * --crop_end_win. Default: 16')
@click.option('--crop_end_gap_win', type=int,
              help='Define window size used to crop end by gap, coupled with --crop_end_gap_thr option. Default: 150')
@click.option('--crop_end_gap_thr', type=float,
              help='Crop end by gap function will check from the beginning and end of each sequence from MSA by '
                   'iteratively choosing a slide window. It will stop until the gap proportion in this slide window is '
                   'smaller than --crop_end_gap_thr. The nucleotide do not match the requirement will be converted to -'
                   ' Default: 0.1')
@click.option('--start_patterns', type=str,
              help='LTR elements will always start with fixed pattern. TE Trimmer will check if it starts with those '
                   'patterns. If not, it will seek around the start point, if the pattern is found, the start point '
                   'will be converted to there. Note: if you want to give multiple start patterns, separate them by comma. '
                   'Like: TG,TA,TC (No space between them). The order of the given patterns is matter. Default: TG')
@click.option('--end_patterns', type=str,
              help='LTR elements will always end with fixed pattern. TE Trimmer will check if it end with those '
                   'patterns. If not, it will seek around the end point, if the pattern is found, the start point will '
                   'be converted to there. Note: if you want to give multiple end patterns, separate them by comma. Like: '
                   'CA,TA,GA (No space between them). The order of the given patterns is matter. Default: AC')
@click.option('--mini_orf', type=int,
              help='Set the minimum ORF length that will be predicted by TE Trimmer. Default: 200')
@click.option('--check_extension_win', type=str,
              help='Define check windows size for extension. Deafault: 150')
@click.option('--num_threads', '-t', default=10, type=int,
              help='Threads numbers used for TE Trimmer. Default: 10')
@click.option('--classify_unknown', default=False, is_flag=True,
              help='Use RepeatClassfier to classify the consensus sequence if the input sequence is not classfied or is unknown. Default: False')
@click.option('--classify_all', default=False, is_flag=True,
              help='Use RepeatClassfier to classify every consensus sequence.  WARNING: it will take longer. Default: False')
def main(input_file, genome_file, output_dir, continue_analysis, pfam_dir, min_blast_len, num_threads, max_msa_lines,
         top_mas_lines, min_seq_num, max_cluster_num, cons_thr, ext_thr, ex_step,
         max_extension, gap_thr, gap_nul_thr, crop_end_thr, crop_end_win, crop_end_gap_thr, crop_end_gap_win,
         start_patterns, end_patterns, mini_orf, species, check_extension_win, merge, genome_anno, hmm,
         keep_intermediate, classify_unknown, classify_all):
    """
        ###########################################################
        TE Trimmer v1.1 (22/SEP/2023)
        Email: jqian@bio1.rwth-aachen.de
        https://github.com/qjiangzhao/TE-Trimmer
        ###########################################################

        python ./path_to_TE_Trimmer_bin/main.py -i <TE_consensus_file> -o <genome_file>


        TE Trimmer is designed to replace transposable element (TE) manual curation. Two mandatory arguments
        are required including <genome file> and <TE consensus file> from TE annotation software like RepeatModeler, EDTA,
        and REPET et al. TE Trimmer can do blast, extension, multiple sequence alignment, and defining TE boundaries.

    """

    start_time = datetime.now()
    print(f"\nTE Trimmer started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    #####################################################################################################
    # Code block: Change permissions of Aliview and TE_Aid
    #####################################################################################################

    # Change TE_Aid permission
    TE_aid_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "TE-Aid-master")
    # Change permissions of the directory and all its content to 755
    # 755 in octal corresponds to rwxr-xr-x
    change_permission = change_permissions_recursive(TE_aid_path, 0o755)

    if not change_permission:
        return

    # Change Aliview permission
    aliview_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "aliview")
    # Change permissions of the directory and all its content to 755
    # 755 in octal corresponds to rwxr-xr-x
    change_permissions_recursive(aliview_path, 0o755)

    #####################################################################################################
    # Code block: Define the default options according to the given species
    #####################################################################################################

    # Load species-specific default values from the JSON config
    default_values = species_config.get(species, {})

    # Override default values with command-line options if provided

    if cons_thr is None:
        cons_thr = default_values.get("cons_thr")

    if max_msa_lines is None:
        max_msa_lines = default_values.get("max_msa_lines")

    if top_mas_lines is None:
        top_mas_lines = default_values.get("top_mas_lines")

    if min_seq_num is None:
        min_seq_num = default_values.get("min_seq_num")

    if min_blast_len is None:
        min_blast_len = default_values.get("min_blast_len")

    if max_cluster_num is None:
        max_cluster_num = default_values.get("max_cluster_num")

    if ext_thr is None:
        ext_thr = default_values.get("ext_thr")

    if ex_step is None:
        ex_step = default_values.get("ex_step")

    if max_extension is None:
        max_extension = default_values.get("max_extension")

    if gap_thr is None:
        gap_thr = default_values.get("gap_thr")

    if gap_nul_thr is None:
        gap_nul_thr = default_values.get("gap_nul_thr")

    if crop_end_thr is None:
        crop_end_thr = default_values.get("crop_end_thr")

    if crop_end_win is None:
        crop_end_win = default_values.get("crop_end_win")

    if crop_end_gap_thr is None:
        crop_end_gap_thr = default_values.get("crop_end_gap_thr")

    if crop_end_gap_win is None:
        crop_end_gap_win = default_values.get("crop_end_gap_win")

    if start_patterns is None:
        start_patterns = default_values.get("start_patterns")

    if end_patterns is None:
        end_patterns = default_values.get("end_patterns")

    if mini_orf is None:
        mini_orf = default_values.get("mini_orf")

    if check_extension_win is None:
        check_extension_win = default_values.get("check_extension_win")

    #####################################################################################################
    # Code block: Define input file, output directory, genome, check blast database
    #####################################################################################################

    # bin_py_path contains all classes and bash code
    # so.path.abspath(__file__) will return the current executable python file
    bin_py_path = os.path.dirname(os.path.abspath(__file__))

    # check if path exist otherwise create one
    os.makedirs(output_dir, exist_ok=True)
    output_dir = os.path.abspath(output_dir)  # get absolute path

    # Check if output directory is empty when --continue_analysis is False
    if os.listdir(output_dir) and not continue_analysis:

        click.echo(
            f"WARNING: The output directory {output_dir} is not empty. Please empty your output directory or "
            f"choose another empty directory\n")

        # Stop the whole program when the output directory is not empty
        return

    # Make a new folder for single fasta sequence
    single_file_dir = os.path.join(output_dir, "Single_fasta_files")
    if not os.path.exists(single_file_dir):
        os.mkdir(single_file_dir)

    # Make a new folder for MSA
    MSA_dir = os.path.join(output_dir, "Multiple_sequence_alignment")
    if not os.path.exists(MSA_dir):
        os.mkdir(MSA_dir)

    # Make a new folder for classification
    classification_dir = os.path.join(output_dir, "Classification")
    if not os.path.exists(classification_dir):
        os.mkdir(classification_dir)

    # make a new folder for HMM file
    if hmm:
        hmm_dir = os.path.join(output_dir, "HMM_files")
        if not os.path.exists(hmm_dir):
            os.mkdir(hmm_dir)

    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"The fasta file {input_file} does not exist.")
    input_file = os.path.abspath(input_file)  # get input absolute path

    # check if genome file exist
    if not os.path.isfile(genome_file):
        raise FileNotFoundError(f"The genome fasta file {genome_file} does not exist.")
    genome_file = os.path.abspath(genome_file)  # get genome absolute path

    # Define the progress_file, finished seq IDs will be stored here
    progress_file = os.path.join(output_dir, "Finished_sequence_name.txt")
    # Check and create progress_file if it doesn't exist
    if not os.path.exists(progress_file):
        with open(progress_file, 'a') as f:

            f.write("input_name,consensus_name,blast_hit_n,cons_MSA_seq_n,cons_full_blast_n,input_length,cons_length,"
                    "input_TE_type,reclassified_type,terminal_repeat,low_copy,status\n")

    # If pfam database isn't provided, create pfam database at TE Trimmer software folder,
    # the database will be downloaded into there.
    # If the pfam_dir is given, but the database can't be found there, TE Trimmer will download pfam database there
    # and generate index file
    if pfam_dir is None:
        pfam_dir = os.path.join(os.path.dirname(bin_py_path), "pfam_database")
        if not os.path.exists(pfam_dir):
            os.mkdir(pfam_dir)

    prepare_pfam_database(pfam_dir)

    # consensus file
    final_con_file = os.path.join(output_dir, "TE_Trimmer_consensus.fasta")
    final_unknown_con_file = os.path.join(classification_dir, "temp_TE_Trimmer_unknown_consensus.fasta")
    final_classified_con_file = os.path.join(classification_dir, "temp_TE_Triimer_classifed_consensus.fasta")

    #####################################################################################################
    # Code block: Merge input file and generate single fasta file
    #####################################################################################################

    # Generate single files when continue_analysis is false
    if not continue_analysis:

        # Do cd-hit-est merge when merge is true and continue_analysis is false
        if merge:
            click.echo("\nTE Trimmer is merging input sequences, this might take some time.\n")
            merge_output = os.path.join(output_dir, f"{input_file}_cd_hit.fa")

            # Set lower identity threshold for the query, this can increase sensitive
            # Merge step will remove single LTR but nested TEs can mask other TEs
            cd_hit_est(input_file, merge_output, identity_thr=0.9, aL=0, aS=0.85, s=0, thread=num_threads)

            # Convert input_file to merged input_file
            input_file = merge_output
            click.echo("Merge finished.\n")

        # Separate fasta to single files, if fasta header contain "/" or " " or ":" convert them to "_"
        # Call this function to separate to single fasta files and create objects from input file
        seq_list = separate_sequences(input_file, single_file_dir, continue_analysis=False)

        # Calculate the total sequence number 
        single_fasta_n = len(seq_list)
        click.echo(f"{single_fasta_n} sequences are detected from the input file")

        # Create new object to check blast database availability
        # Check if blast database and genome length files are available, otherwise create them at the
        # same directory of genome file
        check_database(genome_file)
        # Initial call to print 0% progress
        printProgressBar(0, single_fasta_n, prefix='Progress:', suffix='Complete', length=50)

    else:
        # Check if the can perform continue analysis
        if not os.listdir(single_file_dir):
            
            click.echo("\nWARNING: TE Trimmer can't do continue analysis, please make sure the output directory is same"
                       " with your previous analysis.\n")
            return

        else:

            click.echo("\nTE Trimmer will continue to analyze based on previous results.\n")

            # Create seq_list, which contain sequence objects using the single fasta files.
            seq_list = separate_sequences(input_file, single_file_dir, continue_analysis=True)
            single_fasta_n = len(seq_list)

            # Check which sequences have already been processed
            complete_sequences, skipped_count, low_copy_count = check_progress_file(progress_file)

            # Filter out already complete sequences from the total sequences
            seq_list = [seq for seq in seq_list if seq.name not in complete_sequences]
            click.echo(f"\n{single_fasta_n - len(seq_list)} sequences has been processed previously. "
                       f"TE Trimmer will continue to analyze based on previous results\n")

    #####################################################################################################
    # Code block: Enable multiple threads
    #####################################################################################################
    analyze_sequence_params = [
        (seq, single_file_dir, genome_file, MSA_dir, min_blast_len, min_seq_num, max_msa_lines,
         top_mas_lines, max_cluster_num, cons_thr, ext_thr, ex_step, classification_dir,
         max_extension, gap_thr, gap_nul_thr, crop_end_thr, crop_end_win, crop_end_gap_thr, crop_end_gap_win,
         start_patterns, end_patterns, output_dir, pfam_dir, mini_orf, single_fasta_n, hmm,
         check_extension_win, keep_intermediate, progress_file, classify_unknown, classify_all,
         final_con_file, final_unknown_con_file, final_classified_con_file
         ) for seq in seq_list]

    # Using a ProcessPoolExecutor to run the function in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
         executor.map(analyze_sequence_helper, analyze_sequence_params)
    # multiprocessing
    #with mp.Pool(processes=10) as p:
        #p.starmap(analyze_sequence, analyze_sequence_params)

    #####################################################################################################
    # Code block: Check if all sequences are finished
    #####################################################################################################

    cd_hit_merge_output_final = None
    if os.path.exists(final_con_file):

        # Do cd-hit-est merge when finish all sequence
        cd_hit_merge_output_final = os.path.join(output_dir, "TE_Trimmer_consensus_merged.fasta")
        # According to 80-80-80 rule to filter final consensus sequences

        cd_hit_est(final_con_file, cd_hit_merge_output_final, identity_thr=0.8, aL=0.8, aS=0.8, s=1, thread=num_threads)

    # At the end of the program, check if all sequences have been processed
    completed_sequence, skipped_count, low_copy_count = check_progress_file(progress_file)

    # Calculate the total count
    processed_count = len(completed_sequence)

    if processed_count == single_fasta_n:
        click.echo(f"All sequences have been analysed!\n"
                   f"In the analysed sequences {skipped_count} are skipped.\n"
                   f"In the analysed sequences {low_copy_count} are identified as low copy TE")

    else:
        remaining = single_fasta_n - processed_count
        click.echo(f"{remaining} sequences have not been analysed.")

    t_end_time = datetime.now()
    print("end time:", t_end_time)

    #####################################################################################################
    # Code block: Finish classifying unknown consensus file and writing sequences back to consensus file
    #####################################################################################################

    updated_type = {}
    temp_repeatmasker_dir = os.path.join(classification_dir, "temp_repeatmakser_classification")

    if os.path.exists(final_unknown_con_file) and os.path.exists(final_classified_con_file):

        os.makedirs(temp_repeatmasker_dir, exist_ok=True)
        classification_out = repeatmasker(final_unknown_con_file, final_classified_con_file, temp_repeatmasker_dir,
                                          thread=num_threads, classify=True)

        if classification_out:

            repeatmasker_out = os.path.join(temp_repeatmasker_dir, "temp_TE_Trimmer_unknown_consensus.fasta.out")
            reclassified_dict = repeatmasker_output_classify(repeatmasker_out, progress_file)
            print(reclassified_dict)
    else:
        click.echo("one of the consensus file does not exist, RepeatMasker reclassify does not run")

    #update_cons_file(updated_type, final_unknown_con_file, final_con_file)

    #####################################################################################################
    # Code block: merge consensus_file
    #####################################################################################################

    cd_hit_merge_output_final = None
    if os.path.exists(final_con_file):

        # Do cd-hit-est merge when finish all sequence
        cd_hit_merge_output_final = os.path.join(output_dir, "TE_Trimmer_consensus_merged.fasta")
        # According to 80-80-80 rule to filter final consensus sequences
        cd_hit_est(final_con_file, cd_hit_merge_output_final, identity_thr=0.8, aL=0.8, aS=0.8, s=1, thread=num_threads)

    # Delete MSA_dir and Classification if they are empty
    if not os.listdir(MSA_dir):
        os.rmdir(MSA_dir)

    if not os.listdir(classification_dir):
        os.rmdir(classification_dir)

    # If 95% of the query sequences are processed, RepeatMasker is allowed to be performed
    if processed_count >= single_fasta_n*0.9:
        # Run RepeatMasker
        if genome_anno and cd_hit_merge_output_final:
            click.echo("TE Trimmer is performing whole genome TE annotation by RepeatMasker")
            # make a new folder for RepeatMasker output
            repeatmasker_dir = os.path.join(output_dir, "RepeatMasker_result")
            if not os.path.exists(repeatmasker_dir):
                os.mkdir(repeatmasker_dir)
            genome_anno_result = \
                repeatmasker(genome_file, cd_hit_merge_output_final, repeatmasker_dir, thread=num_threads)
            if genome_anno_result:
                click.echo("Finished whole genome TE annotation by RepeatMasker")

    else:
        click.echo("Less than 90% of the query sequences processed, TE Trimmer can't perform whole genome TE annotation")

    end_time = datetime.now()
    duration = end_time - start_time

    # Remove microseconds from the duration
    duration_without_microseconds = timedelta(days=duration.days, seconds=duration.seconds)
    if not keep_intermediate:
        # Remove all single files when all the sequences are processed
        shutil.rmtree(single_file_dir)

    print(f"\nTE Trimmer finished at {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    print(f"TE Trimmer runtime was {duration_without_microseconds}")
    # printProgressBar(single_fasta_n, single_fasta_n, prefix = 'Progress:', suffix = 'Sequences Processed', length = 50)
    
        
# The following is necessary to make the script executable, i.e., python myscript.py.
if __name__ == '__main__':
    main()
