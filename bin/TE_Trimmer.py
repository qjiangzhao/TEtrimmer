# Standard library imports
import os
import traceback
from datetime import timedelta, datetime
import multiprocessing as mp
import click
import concurrent.futures
import json
import pandas as pd
from Bio import SeqIO

# Local imports
from Function_separate_fasta import separate_sequences
from Function_blast_extension_mafft import blast, remove_files_with_start_pattern, change_permissions_recursive, \
    repeatmasker, check_database, check_bed_uniqueness, extract_fasta, cd_hit_est, handle_sequence_low_copy, \
    handle_sequence_skipped, repeatmasker_output_classify, rename_cons_file, update_low_copy_cons_file, \
    rename_files_based_on_dict
from Class_bed_filter import BEDFile
from Function_def_boundary_and_crop import find_boundary_and_crop
from Class_TE_aid import check_self_alignment
from Function_clean_and_clauster_MSA import clean_and_cluster_MSA
from Function_orf_domain_prediction import prepare_pfam_database


# Define a function to check progress file, which will be used for continue analysis
def check_progress_file(progress_file_path):
    # Read the progress file into a pandas DataFrame
    df = pd.read_csv(progress_file_path)

    # Calculate skipped and low copy element number
    skipped_count = df[df['status'].str.strip().str.lower() == 'skipped'].shape[0]
    low_copy_count = df[df['low_copy'].astype(str).str.strip().str.lower() == 'true'].shape[0]
    unknown_n = df[df['reclassified_type'].str.contains('Unknown', na=False)].shape[0]
    classifid_n = df[df['reclassified_type'].str.contains('/', na=False)].shape[0]

    # Calculate classified_pro
    if classifid_n != 0 and unknown_n != 0:
        classified_pro = classifid_n / (unknown_n + classifid_n)
    else:
        classified_pro = 0  # Set a default value (or any other value you deem appropriate)

    # Get unique 'input_name' values
    local_completed_sequences = df['input_name'].unique().tolist()

    return local_completed_sequences, skipped_count, low_copy_count, classified_pro


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
    click.echo(f'\r{prefix} |{bar}| {iteration}/{total} = {percent}% {suffix}', nl=False)

    # Print New Line on Complete
    if iteration == total:
        click.echo()

#####################################################################################################
# Code block: Define analyze_sequence function
#####################################################################################################


def analyze_sequence_helper(params):
    return analyze_sequence(*params)


def analyze_sequence(seq_obj, genome_file, MSA_dir, min_blast_len, min_seq_num, max_msa_lines,
                     top_mas_lines, max_cluster_num, cons_thr, ext_thr, ext_step, classification_dir,
                     max_ext, gap_thr, gap_nul_thr, crop_end_div_thr, crop_end_div_win, crop_end_gap_thr,
                     crop_end_gap_win, start_patterns, end_patterns, output_dir, pfam_dir, mini_orf,
                     single_fasta_n, hmm, ext_check_win, debug, progress_file,
                     classify_unknown, classify_all, final_con_file, final_unknown_con_file,
                     final_classified_con_file, low_copy_dir, fast_mode, error_files, plot_skip, skipped_dir):

    #####################################################################################################
    # Code block: Set different extension parameters for DNA, SIAN, Helitron, and MITE elements
    #####################################################################################################
    try:
        # Get query fasta file path
        seq_name = seq_obj.get_seq_name()
        seq_type = seq_obj.get_old_TE_type()
        seq_file = seq_obj.get_input_fasta()

        # Due to DNA element will be much shorter than LTR and LINE elements, set different parameters
        if "DNA" in seq_type:
            ext_step = 500
            max_ext = 3500
            min_blast_len = 150
            crop_end_gap_win = 100
            ext_check_win = 50

        # The average length of SINE element is around 500 bp, give different default parameters
        if "SINE" in seq_type:
            ext_step = 200
            max_ext = 1400
            min_blast_len = 80
            crop_end_gap_win = 50
            ext_check_win = 50

        if "Helitron" in seq_type:
            ext_step = 500
            max_ext = 3500
            min_blast_len = 150
            crop_end_gap_win = 100
            ext_check_win = 50

        if "MITE" in seq_type:
            ext_step = 100
            max_ext = 500
            min_blast_len = 50
            crop_end_gap_win = 40
            ext_check_win = 50

        # run blast for each single fasta file and return a bed file absolute path
        bed_out_file_dup, blast_hits_count, blast_out_file = blast(seq_file, genome_file, MSA_dir,
                                                                   min_length=min_blast_len, seq_obj=seq_obj)

    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"Error while running blast for sequence: {seq_name}\n")
            f.write(tb_content + '\n\n')
        click.echo(f"Error while running blast for sequence: {seq_name}. Main Error: {str(e)}")
        return

    try:

        # Check if blast hit number is equal 0, then skip this sequence
        if blast_hits_count == 0:
            click.echo(f"\n{seq_name} is skipped due to blast hit number is 0\n")
            handle_sequence_skipped(seq_obj, progress_file, debug, MSA_dir, classification_dir)
            return

        # Check if blast hit number is smaller than "min_seq_num", not include "min_seq_num"
        elif blast_hits_count != 0 and blast_hits_count < min_seq_num:
            check_low_copy, blast_full_length_n, found_match, TE_aid_plot = check_self_alignment(
                seq_obj, seq_file, MSA_dir, genome_file, blast_hits_count, blast_out_file, plot_skip=plot_skip)

            if check_low_copy is True:
                # Update terminal repeat and blast full length number, remove low copy intermediate files
                handle_sequence_low_copy(seq_obj, progress_file, debug, MSA_dir,
                                         classification_dir, found_match=found_match,
                                         blast_full_length_n=blast_full_length_n)

                # Integrate low copy element sequence into consensus file
                update_low_copy_cons_file(seq_obj, final_con_file, final_unknown_con_file,
                                          final_classified_con_file, low_copy_dir, TE_aid_plot)
            else:
                click.echo(f"\n{seq_name} is skipped due to blast hit number is smaller than {min_seq_num} "
                           f"and check_low_copy is {check_low_copy}\n")

                # handle_sequence_skipped will update skipped status
                handle_sequence_skipped(seq_obj, progress_file, debug, MSA_dir, classification_dir,
                                        plot_skip=plot_skip, te_aid_plot=TE_aid_plot, skip_proof_dir=skipped_dir)
            return  # when blast hit number is smaller than 10, code will execute next fasta file

    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"\nError while checking low copy for sequence: {seq_name}\n")
            f.write(tb_content + '\n\n')
        click.echo(f"\nError while checking low copy for sequence: {seq_name}. Error: {str(e)}\n")
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
                                                   cluster_col_thr=100, fast_mode=fast_mode)
    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"Error when group MSA: {seq_name}\n")
            f.write(tb_content + '\n\n')
        click.echo(f"\nError when group MSA for sequence: {seq_name}. Error: {str(e)}\n")
        return

    try:
        # cluster false means no cluster, TE Trimmer will skip this sequence.
        if cluster_MSA_result is False:
            check_low_copy, blast_full_length_n, found_match, TE_aid_plot = check_self_alignment(
                seq_obj, seq_file, MSA_dir, genome_file, blast_hits_count, blast_out_file, plot_skip=plot_skip)

            if check_low_copy is True:
                # Update terminal repeat and blast full length number, remove low copy intermediate files
                handle_sequence_low_copy(seq_obj, progress_file, debug, MSA_dir,
                                         classification_dir, found_match=found_match,
                                         blast_full_length_n=blast_full_length_n)
                update_low_copy_cons_file(seq_obj, final_con_file, final_unknown_con_file,
                                          final_classified_con_file, low_copy_dir, TE_aid_plot)

            else:
                click.echo(
                    f"\n{seq_name} is skipped due to sequence number in each cluster is smaller "
                    f"than {min_seq_num} and check_low_copy is {check_low_copy}\n")
                handle_sequence_skipped(seq_obj, progress_file, debug, MSA_dir, classification_dir,
                                        plot_skip=plot_skip, te_aid_plot=TE_aid_plot, skip_proof_dir=skipped_dir)
            return

        # Cluster True means not necessary to cluster MSA, perform find_boundary_and_crop directly
        elif cluster_MSA_result is True:
            find_boundary_result = find_boundary_and_crop(
                bed_out_filter_file, genome_file, MSA_dir, pfam_dir, seq_obj, hmm,
                classify_all, classify_unknown, error_files,
                cons_threshold=cons_thr, ext_threshold=ext_thr,
                ex_step_size=ext_step, max_extension=max_ext,
                gap_threshold=gap_thr, gap_nul_thr=gap_nul_thr,
                crop_end_thr=crop_end_div_thr, crop_end_win=crop_end_div_win,
                crop_end_gap_thr=crop_end_gap_thr, crop_end_gap_win=crop_end_gap_win,
                start_patterns=start_patterns, end_patterns=end_patterns, mini_orf=mini_orf,
                define_boundary_win=ext_check_win, fast_mode=fast_mode
            )
            if find_boundary_result == "Short_sequence":
                click.echo(f"\n{seq_name} is skipped due to too short length\n")
                handle_sequence_skipped(seq_obj, progress_file, debug, MSA_dir, classification_dir)
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
                                                                 cluster_col_thr=100, fast_mode=fast_mode)
                # inner_cluster_MSA_result is false means this cluster sequence number is too less
                if inner_cluster_MSA_result is False:
                    continue
                elif inner_cluster_MSA_result is True:  # Means don't need to cluster

                    inner_find_boundary_result = find_boundary_and_crop(
                        cluster_bed_files_list[i], genome_file, MSA_dir, pfam_dir, seq_obj,
                        hmm, classify_all, classify_unknown, error_files, cons_threshold=cons_thr,
                        ext_threshold=ext_thr, ex_step_size=ext_step, max_extension=max_ext,
                        gap_threshold=gap_thr, gap_nul_thr=gap_nul_thr,
                        crop_end_thr=crop_end_div_thr, crop_end_win=crop_end_div_win,
                        crop_end_gap_thr=crop_end_gap_thr, crop_end_gap_win=crop_end_gap_win,
                        start_patterns=start_patterns, end_patterns=end_patterns,
                        mini_orf=mini_orf, define_boundary_win=ext_check_win, fast_mode=fast_mode)

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
                            pfam_dir, seq_obj, hmm, classify_all, classify_unknown, error_files,
                            cons_threshold=cons_thr, ext_threshold=ext_thr,
                            ex_step_size=ext_step, max_extension=max_ext,
                            gap_threshold=gap_thr, gap_nul_thr=gap_nul_thr,
                            crop_end_thr=crop_end_div_thr, crop_end_win=crop_end_div_win,
                            crop_end_gap_thr=crop_end_gap_thr, crop_end_gap_win=crop_end_gap_win,
                            start_patterns=start_patterns, end_patterns=end_patterns,
                            mini_orf=mini_orf, define_boundary_win=ext_check_win, fast_mode=fast_mode)
                        if inner_inner_find_boundary_result == "Short_sequence":
                            continue
                        elif inner_inner_find_boundary_result:
                            all_inner_skipped = False
                        elif not inner_inner_find_boundary_result:  # This means the errors happen in the function
                            continue

            # Check the flag after the loop. If all inner clusters were skipped, write the progress file
            if all_inner_skipped:
                check_low_copy, blast_full_length_n, found_match, TE_aid_plot = check_self_alignment(
                    seq_obj, seq_file, MSA_dir, genome_file, blast_hits_count, blast_out_file, plot_skip=plot_skip)

                if check_low_copy is True:
                    # Update terminal repeat and blast full length number, remove low copy intermediate files
                    handle_sequence_low_copy(seq_obj, progress_file, debug, MSA_dir,
                                             classification_dir, found_match=found_match,
                                             blast_full_length_n=blast_full_length_n)
                    update_low_copy_cons_file(seq_obj, final_con_file, final_unknown_con_file,
                                              final_classified_con_file, low_copy_dir, TE_aid_plot)
                else:
                    handle_sequence_skipped(seq_obj, progress_file, debug, MSA_dir, classification_dir,
                                            plot_skip=plot_skip, te_aid_plot=TE_aid_plot, skip_proof_dir=skipped_dir)
                    click.echo(
                        f"\n{seq_name} is skipped due to sequence number in second round each cluster is "
                        f"smaller than {min_seq_num} or the sequence is too short and check_low_copy is {check_low_copy}\n")
                return

    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"Error during boundary finding and cropping for sequence: {seq_name}\n")
            f.write(tb_content + '\n\n')
        click.echo(f"Error during boundary finding and cropping for sequence: {seq_name}. Error: {str(e)}\n")
        return

    # After all processing is done, change status to process and write the name of the file to the progress file
    seq_obj.update_status("processed", progress_file)

    # If all this sequence is finished remove all files contain this name
    if not debug:
        remove_files_with_start_pattern(MSA_dir, seq_name)
        remove_files_with_start_pattern(classification_dir, seq_name)

    # Read and count sequences from Finished_sequence_name.txt
    completed_sequence, skipped_count, low_copy_count, classified_pro = check_progress_file(progress_file)

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


@click.command(context_settings=dict(max_content_width=120),
               help="""\b
               ##########################################################################################
               \b
               ████████╗███████╗    ████████╗██████╗ ██╗███╗   ███╗███╗   ███╗███████╗██████╗ 
               ╚══██╔══╝██╔════╝    ╚══██╔══╝██╔══██╗██║████╗ ████║████╗ ████║██╔════╝██╔══██╗
                  ██║   █████╗         ██║   ██████╔╝██║██╔████╔██║██╔████╔██║█████╗  ██████╔╝
                  ██║   ██╔══╝         ██║   ██╔══██╗██║██║╚██╔╝██║██║╚██╔╝██║██╔══╝  ██╔══██╗
                  ██║   ███████╗       ██║   ██║  ██║██║██║ ╚═╝ ██║██║ ╚═╝ ██║███████╗██║  ██║
                  ╚═╝   ╚══════╝       ╚═╝   ╚═╝  ╚═╝╚═╝╚═╝     ╚═╝╚═╝     ╚═╝╚══════╝╚═╝  ╚═╝
                  
                Version: v1.1 (27/OCT/2023) 

                Github: https://github.com/qjiangzhao/TE-Trimmer

                Developers:                                                                                                       
                Jiangzhao Qian;  RWTH Aachen University;                Email: jqian@bio1.rwth-aachen.de                          
                Hang Xue;        University of California, Berkeley;    Email: hang_xue@berkeley.edu

                Funding source:                                                                                         
                Panstruga's Lab. Website: https://www.bio1.rwth-aachen.de/PlantMolCellBiology/index.html                 

                Many thanks to Dr. Stefan Kusch                                                           

                ##########################################################################################              

                python ./path_to_TE_Trimmer_bin/TE_Trimmer.py -i <TE_consensus_file> -o <genome_file>

                TE Trimmer is designed to replace transposable element (TE) manual curation. 

                Two mandatory arguments are required including <genome file> and <TE consensus file> from TE 
                annotation software like RepeatModeler or EDTA et al. TE Trimmer can do blast, extension, multiple sequence alignment, and defining TE boundaries.

""")
@click.option('--input_file', '-i', required=True, type=str,
              help='Path to TE consensus file (FASTA format). Use the output from RepeatModeler or EDTA et al.')
@click.option('--genome_file', '-g', required=True, type=str,
              help='Path to genome FASTA file.')
@click.option('--output_dir', '-o', default=os.getcwd(), type=str,
              help='Output directory. Default: current directory.')
@click.option('--species', '-s', default='fungi', type=click.Choice(species_config.keys()),
              help='Select the species for which you want to run TE Trimmer.')
@click.option('--continue_analysis', '-ca', default=False, is_flag=True,
              help='Continue to analysis after interruption.')
@click.option('--dedup', default=False, is_flag=True,
              help='Remove duplicate sequences in input file.')
@click.option('--genome_anno', default=False, is_flag=True,
              help='Perform genome TE annotation using the TE Trimmer curated database. Requires RepeatMasker.')
@click.option('--hmm', default=False, is_flag=True,
              help='Generate HMM files for each consensus sequences.')
@click.option('--debug', default=False, is_flag=True,
              help='Open debug mode. This will keep all raw files. WARNING: Many files will be produced.')
@click.option('--fast_mode', default=False, is_flag=True,
              help='Reduce running time but at the cost of lower accuracy and specificity.')
@click.option('--plot_skip', default=False, is_flag=True,
              help='Perform TE_Aid plot for skipped elements')
@click.option('--pfam_dir', default=None, type=str,
              help='Pfam database directory. Omit this option if you do not have a local PFAM database. '
                   'TE Trimmer will download the database automatically in this case.')
@click.option('--cons_thr', type=float,
              help='Threshold used for the final consensus sequence generation. Default: 0.8')
@click.option('--max_msa_lines', type=int,
              help='Set the maximum sequences number for multiple sequence alignment. Default: 100')
@click.option('--top_mas_lines', type=int,
              help='When the sequence number of multiple sequence alignment (MSA) is greater than <--max_msa_lines>, '
                   'TE Trimmer will sort sequences by length and choose <--top_msa_lines> number '
                   'of sequences. Then, TE Trimmer will randomly select sequences from all remaining BLAST hits until '
                   '<--max_msa_lines> sequences are found for the multiple sequence alignment. Default: 70')
@click.option('--min_seq_num', type=int,
              help='The minimum sequence number for each multiple sequence alignment. Default: 10')
@click.option('--min_blast_len', type=int,
              help='The minimum sequence length for blast hits. Default: 150')
@click.option('--max_cluster_num', type=int,
              help='The maximum cluster number for each multiple sequence alignment. Each multiple '
                   'sequence alignment can be divided into different clusters. TE Trimmer will sort '
                   'clusters by sequence number and choose the top <--max_cluster_num> of clusters for '
                   'further analysis. WARNING: Big number will dramatically increase running time. Default: 2')
@click.option('--ext_thr', type=float,
              help="threshold used for define the extension extent. The lower the value of <--ext_thr>, the easier the "
                   "extensions on both ends be longer. Reduce <--ext_thr> if TE Trimmer fails to determine the correct "
                   "ends of repeat elements. Default: 0.7")
@click.option('--ext_check_win', type=str,
              help='Define check windows size for extension. Default: 150')
@click.option('--ext_step', type=int,
              help='Number of nucleotides to be added to the left and right ends of the multiple sequence alignment. '
                   'TE_Trimmer will iteratively add <--ext_step> number of nucleotides until finding the boundary. '
                   'Default: 1000')
@click.option('--max_ext', type=int,
              help='The maximum extension in nucleotides at both ends of the multiple sequence alignment. Default: 7000')
@click.option('--gap_thr', type=float,
              help='If multiple sequence alignment positions (columns) have a gap proportion larger than <--gap_thr> '
                   'and the proportion of the most common nucleotide in this column is less than <--gap_nul_thr>, '
                   'this column will be removed from the consensus. Default: 0.4')
@click.option('--gap_nul_thr', type=float,
              help='Set nucleotide proportion threshold for keeping the column of the multiple sequence alignment. '
                   'Used with the <--gap_thr> option. Default: 0.7')
@click.option('--crop_end_div_thr', type=float,
              help='The crop end by divergence function will convert each nucleotide in the multiple sequence '
                   'alignment into a proportion value. This function will iteratively choose a sliding window from '
                   'each end of each sequence of the MSA and sum up the proportion numbers in this window. '
                   'The cropping will continue until the sum of proportions is larger than <--crop_end_div_thr>. '
                   'Cropped nucleotides will be converted to -. Default: 0.8')
@click.option('--crop_end_div_win', type=int,
              help='Window size used for the end-cropping process. Used with --crop_end_div_thr option. Default: 20')
@click.option('--crop_end_gap_thr', type=float,
              help='The crop end by gap function will iteratively choose a sliding window from '
                   'each end of each sequence of the MSA and calculate the gap proportion in this window. '
                   'The cropping will continue until the sum of gap proportions is smaller than <--crop_end_gap_thr>. '
                   'Cropped nucleotides will be converted to -. Default: 0.1')
@click.option('--crop_end_gap_win', type=int,
              help='Define window size used to crop end by gap, used with <--crop_end_gap_thr> option. Default: 250')
@click.option('--start_patterns', type=str,
              help='LTR elements always start with a conserved sequence pattern. TE Trimmer searches the '
                   'beginning of the consensus sequence for these patterns. If the pattern is not found, '
                   'it will extend the search of <--start_patterns> to up to 15 nucleotides from the beginning '
                   'of the consensus sequence and redefine the start of the consensus sequence '
                   'if the pattern is found. Note: The user can provide multiple LTR start patterns in a '
                   'comma-separated list, like: TG,TA,TC (no spaces; the order of patterns determines '
                   'the priority for the search). Default: TG')
@click.option('--end_patterns', type=str,
              help='LTR elements always end with a conserved sequence pattern. TE Trimmer searches the '
                   'end of the consensus sequence for these patterns. If the pattern is not found, '
                   'it will extend the search of <--end_patterns> to up to 15 nucleotides from the end '
                   'of the consensus sequence and redefine the end of the consensus sequence '
                   'if the pattern is found. Note: The user can provide multiple LTR end patterns in a '
                   'comma-separated list, like: CA,TA,GA (no spaces; the order of patterns determines '
                   'the priority for the search). Default: CA')
@click.option('--mini_orf', type=int,
              help='Define the minimum ORF length that will be predicted by TE Trimmer. Default: 200')
@click.option('--num_threads', '-t', default=10, type=int,
              help='Threads numbers used for TE Trimmer. Default: 10')
@click.option('--classify_unknown', default=False, is_flag=True,
              help='Use RepeatClassifier to classify the consensus sequence if the input sequence is not classified or '
                   'is unknown. Default: False')
@click.option('--classify_all', default=False, is_flag=True,
              help='Use RepeatClassifier to classify every consensus sequence.  WARNING: it will take longer time. '
                   'Default: False')
def main(input_file, genome_file, output_dir, continue_analysis, pfam_dir, min_blast_len, num_threads, max_msa_lines,
         top_mas_lines, min_seq_num, max_cluster_num, cons_thr, ext_thr, ext_step, plot_skip,
         max_ext, gap_thr, gap_nul_thr, crop_end_div_thr, crop_end_div_win, crop_end_gap_thr, crop_end_gap_win,
         start_patterns, end_patterns, mini_orf, species, ext_check_win, dedup, genome_anno, hmm,
         debug, fast_mode, classify_unknown, classify_all):

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

        # Convert string "False" to boolean
        if max_cluster_num == "False":
            max_cluster_num = False

        if fast_mode:
            max_cluster_num = 2

    if ext_thr is None:
        ext_thr = default_values.get("ext_thr")

    if ext_step is None:
        ext_step = default_values.get("ext_step")

    if max_ext is None:
        max_ext = default_values.get("max_ext")

    if gap_thr is None:
        gap_thr = default_values.get("gap_thr")

    if gap_nul_thr is None:
        gap_nul_thr = default_values.get("gap_nul_thr")

    if crop_end_div_thr is None:
        crop_end_div_thr = default_values.get("crop_end_div_thr")

    if crop_end_div_win is None:
        crop_end_div_win = default_values.get("crop_end_div_win")

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

    if ext_check_win is None:
        ext_check_win = default_values.get("ext_check_win")

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

    # Define proof_annotation folder path
    proof_annotation_dir = os.path.join(output_dir, "TE_Trimmer_for_proof_annotation")
    os.makedirs(proof_annotation_dir, exist_ok=True)

    # Define low copy folder
    low_copy_dir = os.path.join(proof_annotation_dir, "Low_copy_TE")
    os.makedirs(low_copy_dir, exist_ok=True)

    # Define skipped folder if it is required
    if plot_skip:
        skipped_dir = os.path.join(proof_annotation_dir, "Skipped_TE")
        os.makedirs(skipped_dir, exist_ok=True)
    else:
        skipped_dir = None

    # Define proof annotation evaluation folders
    perfect_proof = os.path.join(proof_annotation_dir, "Perfect_annotation")
    good_proof = os.path.join(proof_annotation_dir, "Good_annotation")
    intermediate_proof = os.path.join(proof_annotation_dir, "Recommend_check_annotation")
    need_check_proof = os.path.join(proof_annotation_dir, "Need_check_annotation")
    os.makedirs(perfect_proof, exist_ok=True)
    os.makedirs(good_proof, exist_ok=True)
    os.makedirs(intermediate_proof, exist_ok=True)
    os.makedirs(need_check_proof, exist_ok=True)

    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"The fasta file {input_file} does not exist.")
    input_file = os.path.abspath(input_file)  # get input absolute path

    # check if genome file exist
    if not os.path.isfile(genome_file):
        raise FileNotFoundError(f"The genome fasta file {genome_file} does not exist.")
    genome_file = os.path.abspath(genome_file)  # get genome absolute path

    # Define the progress_file, finished seq IDs will be stored here
    progress_file = os.path.join(output_dir, "Finished_sequence_recording.txt")
    # Check and create progress_file if it doesn't exist
    if not os.path.exists(progress_file):
        with open(progress_file, 'a') as f:
            f.write("input_name,consensus_name,blast_hit_n,cons_MSA_seq_n,cons_full_blast_n,input_length,cons_length,"
                    "input_TE_type,reclassified_type,terminal_repeat,low_copy,evaluation,status\n")

    # Define error files to store not mandatory function errors including RepeatClassified classification,
    # RepeatMasker classification, PFAM scanning, muscle alignment
    error_files = os.path.join(output_dir, "error_file.txt")

    # If pfam database isn't provided, create pfam database at TE Trimmer software folder,
    # the database will be downloaded into there.
    # If the pfam_dir is given, but the database can't be found there, TE Trimmer will download pfam database there
    # and generate index file
    if pfam_dir is None:
        pfam_dir = os.path.join(os.path.dirname(bin_py_path), "pfam_database")
        if not os.path.exists(pfam_dir):
            os.mkdir(pfam_dir)
    try:
        if_pfam = prepare_pfam_database(pfam_dir)

        if not if_pfam:  # Check if if_pfam is False
            raise ValueError("PFAM database preparation failed.")  # Raise an exception to be caught below

    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"PFAM database building error\n")
            f.write(tb_content + '\n\n')

        click.echo(f"Note: Can't download PFAM database from internet, please use your own PFAM database\n"
                   f"For example: --pfam_dir <your_PFAM_directory>\n"
                   f"Your PFAM directory should contain: \n"
                   f"Pfam-A.hmm\n"
                   f"Pfam-A.hmm.h3f\n"
                   f"Pfam-A.hmm.h3m\n"
                   f"Pfam-A.hmm.dat\n"
                   f"Pfam-A.hmm.h3i\n"
                   f"Pfam-A.hmm.h3p\n\n"
                   f"PFAM prediction will be used to determine the direction of TEs, it is necessary to make it run\n")
        return

    # Define consensus files. temp files will be used for the final RepeatMasker classification
    final_con_file = os.path.join(output_dir, "TE_Trimmer_consensus.fasta")
    final_unknown_con_file = os.path.join(classification_dir, "temp_TE_Trimmer_unknown_consensus.fasta")
    final_classified_con_file = os.path.join(classification_dir, "temp_TE_Trimmer_classified_consensus.fasta")

    #####################################################################################################
    # Code block: Merge input file and generate single fasta file
    #####################################################################################################

    # Generate single files when continue_analysis is false
    if not continue_analysis:

        # Do cd-hit-est merge when merge is true and continue_analysis is false
        if dedup:
            click.echo("\nTE Trimmer is merging input sequences, this might take some time.\n")
            merge_output = os.path.join(output_dir, f"{input_file}_cd_hit.fa")

            # Set lower identity threshold for the query, this can increase sensitive
            # Merge step will remove single LTR but nested TEs can mask other TEs
            cd_iden_thr = 0.9
            cd_alin_s = 0.85

            if fast_mode:
                cd_iden_thr = 0.8
                cd_alin_s = 0.8

            cd_hit_est(input_file, merge_output, identity_thr=cd_iden_thr, aL=0, aS=cd_alin_s, s=0, thread=num_threads)

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
            complete_sequences, skipped_count, low_copy_count, classified_pro = check_progress_file(progress_file)

            # Filter out already complete sequences from the total sequences
            seq_list = [seq for seq in seq_list if seq.name not in complete_sequences]
            click.echo(f"\n{single_fasta_n - len(seq_list)} sequences has been processed previously.\n")

    #####################################################################################################
    # Code block: Enable multiple threads
    #####################################################################################################
    analyze_sequence_params = [
        (seq, genome_file, MSA_dir, min_blast_len, min_seq_num, max_msa_lines,
         top_mas_lines, max_cluster_num, cons_thr, ext_thr, ext_step, classification_dir,
         max_ext, gap_thr, gap_nul_thr, crop_end_div_thr, crop_end_div_win, crop_end_gap_thr, crop_end_gap_win,
         start_patterns, end_patterns, output_dir, pfam_dir, mini_orf, single_fasta_n, hmm,
         ext_check_win, debug, progress_file, classify_unknown, classify_all,
         final_con_file, final_unknown_con_file, final_classified_con_file, low_copy_dir, fast_mode, error_files,
         plot_skip, skipped_dir
         ) for seq in seq_list]

    # Using a ProcessPoolExecutor to run the function in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        executor.map(analyze_sequence_helper, analyze_sequence_params)
    # multiprocessing
    # with mp.Pool(processes=10) as p:
    # p.starmap(analyze_sequence, analyze_sequence_params)

    #####################################################################################################
    # Code block: Check if all sequences are finished
    #####################################################################################################

    # At the end of the program, check if all sequences have been processed
    completed_sequence, skipped_count, low_copy_count, classified_pro = check_progress_file(progress_file)

    # Calculate the total count
    processed_count = len(completed_sequence)

    if processed_count == single_fasta_n:
        click.echo(f"All sequences have been analysed!\n"
                   f"In the analysed sequences {skipped_count} are skipped.\n"
                   f"In the analysed sequences {low_copy_count} are identified as low copy TE\n"
                   f"TE Trimmer is doing the final classification, merging, or whole genome annotation step.")

    else:
        remaining = single_fasta_n - processed_count
        click.echo(f"\n\n{remaining} sequences have not been analysed.\n")

    #####################################################################################################
    # Code block: Finish classifying unknown consensus file and writing sequences back to consensus file
    #####################################################################################################

    # Suppress RepeatMasker final classification under fast_mode
    if fast_mode:
        classified_pro = 0.01

    # Final RepeatMasker classification isn't necessary, skip it when errors are there
    try:
        if 0.3 <= classified_pro < 0.99:
            temp_repeatmasker_dir = os.path.join(classification_dir, "temp_repeatmakser_classification")

            if os.path.exists(final_unknown_con_file) and os.path.exists(final_classified_con_file):
                os.makedirs(temp_repeatmasker_dir, exist_ok=True)
                classification_out = repeatmasker(final_unknown_con_file, final_classified_con_file,
                                                  temp_repeatmasker_dir,
                                                  thread=num_threads, classify=True)

                if classification_out:
                    repeatmasker_out = os.path.join(temp_repeatmasker_dir,
                                                    "temp_TE_Trimmer_unknown_consensus.fasta.out")
                    reclassified_dict = repeatmasker_output_classify(repeatmasker_out, progress_file,
                                                                     min_iden=60, min_len=80, min_cov=0.5)
                    if reclassified_dict:
                        click.echo(
                            f"{len(reclassified_dict)} TE elements were re-classified by final classification module")

                        # Update final consensus file
                        rename_cons_file(final_con_file, reclassified_dict)
                        rename_files_based_on_dict(proof_annotation_dir, reclassified_dict)
                        rename_files_based_on_dict(perfect_proof, reclassified_dict)
                        rename_files_based_on_dict(good_proof, reclassified_dict)
                        rename_files_based_on_dict(intermediate_proof, reclassified_dict)
                        rename_files_based_on_dict(need_check_proof, reclassified_dict)
                        rename_files_based_on_dict(low_copy_dir, reclassified_dict, seq_name=True)
                        if hmm:
                            rename_files_based_on_dict(hmm_dir, reclassified_dict)

            else:
                click.echo("One of the consensus file does not exist, RepeatMasker reclassify does not run, "
                           "This won't affect the final consensus sequences.")
    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"\nFinal RepeatMasker classification is wrong.\n")
            f.write(tb_content + '\n\n')

        click.echo("Note: The final TE classification is skipped, this won't affect the TE consensus result.")

    #####################################################################################################
    # Code block: merge consensus_file to remove duplications
    #####################################################################################################

    try:
        click.echo("TE Trimmer is removing sequence duplications")

        # Do first round cd-hit-est
        cd_hit_merge_output_round1 = os.path.join(classification_dir, "TE_Trimmer_consensus_merged_round1.fasta")
        cd_hit_merge_output_round1_clstr = f"{cd_hit_merge_output_round1}.clstr"

        # Round 1 merge only require that the alignment coverage for the shorter sequence is greater than 0.8
        # and the similarity is greater than 0.85
        cd_hit_est(final_con_file, cd_hit_merge_output_round1, identity_thr=0.85, aS=0.8, s=0, thread=num_threads)

        # Read progress file
        progress_df = pd.read_csv(progress_file)

        # Create a dictionary with sequence names as keys
        sequence_info = {}
        for index, row in progress_df.iterrows():
            sequence_name = row["consensus_name"]
            evaluation = row["evaluation"] if pd.notna(row["evaluation"]) else "Unknown"  # Default value for NaN
            te_type = row["reclassified_type"] if pd.notna(row["reclassified_type"]) else "Unknown"
            length = row["cons_length"] if pd.notna(row["cons_length"]) else 0  # Default value for NaN
            sequence_info[sequence_name] = {"evaluation": evaluation, "type": te_type, "length": length}

        # Parse cd-hit-est result
        clusters = {}
        current_cluster = []
        with open(cd_hit_merge_output_round1_clstr, "r") as f:
            for line in f:
                if line.startswith(">Cluster"):
                    cluster_name = line.strip().replace(" ", "")  # Remove the empty space in the cluster name
                    if current_cluster:  # If the current_cluster isn't empty
                        clusters[cluster_name] = current_cluster
                        current_cluster = []
                else:
                    seq_info = line.split(">")[1].split("...")[0].split("#")[0]
                    current_cluster.append(seq_info)
            if current_cluster:
                clusters[cluster_name] = current_cluster

        # Check if "Perfect" and "Good" level sequences are included inside cluster and choose the longest one
        best_sequences = []  # Define list to store Perfect or Good sequence names

        # Define list to store sequence in clusters that don't contain Perfect or Good sequences
        sequence_for_round2 = []
        for cluster_name, sequences in clusters.items():
            best_seq = None
            best_length = 0
            for seq in sequences:
                if seq in sequence_info:
                    evaluation = sequence_info[seq]["evaluation"]
                    length = sequence_info[seq]["length"]
                    if evaluation in ["Perfect", "Good"] and length > best_length:
                        best_length = length
                        best_seq = seq
            if best_seq:
                best_sequences.append(best_seq)
            else:
                sequence_for_round2.extend(clusters[cluster_name])  # extend() will make a flat a list

        # Read the original consensus file
        consensus_sequences = SeqIO.parse(final_con_file, "fasta")

        # Define temporary file to store Perfect or Good sequences
        temp_consensus_round1 = os.path.join(classification_dir, "temp_consensus_round1.fasta")

        # Define temporary file to store rest sequence for round2 cd-hit-est
        temp_consensus_round2_input = os.path.join(classification_dir, "temp_consensus_round2_input.fasta")

        # Write sequences to files
        with open(temp_consensus_round1, "w") as high_quality_file, \
                open(temp_consensus_round2_input, 'w') as round2_file:

            for seq_record in consensus_sequences:

                # Sequence names in best_sequences and sequence_for_round2 dont' contain classification
                seq_id = seq_record.id.split("#")[0]

                if seq_id in best_sequences:
                    SeqIO.write(seq_record, high_quality_file, "fasta")
                elif seq_id in sequence_for_round2:
                    SeqIO.write(seq_record, round2_file, 'fasta')

        # Do second round cd-hit-est based on temp_consensus_round2_input
        cd_hit_merge_output_round2 = os.path.join(classification_dir, "TE_Trimmer_consensus_merged_round2.fasta")

        # Round 2 merge require that the alignment coverage for the long and short sequence are both greater than 0.8
        # and the similarity is greater than 0.85
        cd_hit_est(temp_consensus_round2_input, cd_hit_merge_output_round2, identity_thr=0.85, aL=0.8, aS=0.8, s=0.8,
                   thread=num_threads)

        # Define merged file
        cd_hit_est_final_merged = os.path.join(output_dir, "TE_Trimmer_consensus_merged.fasta")

        # Combine the two files into merged file
        with open(temp_consensus_round1, 'r') as file1, \
                open(cd_hit_merge_output_round2, 'r') as file2, \
                open(cd_hit_est_final_merged, 'w') as combined_file:
            # Write contents of the first file
            for line in file1:
                combined_file.write(line)

            # Write contents of the second file
            for line in file2:
                combined_file.write(line)

        # Find sequence names that aren't included inside in cd_hit_est_final_merged file
        # Parse the sequences in the original and merged files
        original_sequences = SeqIO.parse(final_con_file, "fasta")
        merged_sequences = SeqIO.parse(cd_hit_est_final_merged, "fasta")

        # Extract sequence IDs from both files and store to set
        original_ids = {seq_record.id.split("#")[0] for seq_record in original_sequences}
        merged_ids = {seq_record.id.split("#")[0] for seq_record in merged_sequences}

        # Find the difference between the two sets to get sequence names not included in the merged file
        missing_ids = original_ids - merged_ids
        click.echo(missing_ids)

        # Based on missing_ids delete files in proof annotation folder and HMM folder
        for missing_id in missing_ids:
            evaluation_leve = sequence_info[missing_id]["evaluation"]
            if evaluation_leve == "Perfect":
                remove_files_with_start_pattern(perfect_proof, missing_id)
            elif evaluation_leve == "Good":
                remove_files_with_start_pattern(good_proof, missing_id)
            elif evaluation_leve == "Reco_check":
                remove_files_with_start_pattern(intermediate_proof, missing_id)
            elif evaluation_leve == "Need_check":
                remove_files_with_start_pattern(need_check_proof, missing_id)
            else:
                remove_files_with_start_pattern(low_copy_dir, missing_id)

        if hmm:
            remove_files_with_start_pattern(hmm_dir, missing_ids)
        click.echo("\nFinished to remove sequence duplications.\n")

    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"Final cd-hit-est deduplication error\n")
            f.write(tb_content + '\n\n')
            click.echo("\nThe final cd-hit-est merge step can't be performed. Please remove duplicated sequence"
                       " by yourself\n")

    #####################################################################################################
    # Code block: Whole genome TE annotation
    #####################################################################################################

    # Delete MSA_dir and Classification if they are empty
    if not os.listdir(MSA_dir):
        os.rmdir(MSA_dir)

    if not os.listdir(classification_dir):
        os.rmdir(classification_dir)

    try:
        # If 90% of the query sequences are processed, RepeatMasker is allowed to be performed whole genome annotation
        #if processed_count >= single_fasta_n * 0.9:

        # Run RepeatMasker
        if genome_anno:
            click.echo("\nTE Trimmer is performing whole genome TE annotation by RepeatMasker\n")

            if os.path.exists(cd_hit_est_final_merged):
                repeatmakser_lib = cd_hit_est_final_merged
            else:
                repeatmakser_lib = final_con_file

            # make a new folder for RepeatMasker output
            repeatmasker_dir = os.path.join(output_dir, "RepeatMasker_result")
            if not os.path.exists(repeatmasker_dir):
                os.mkdir(repeatmasker_dir)
            genome_anno_result = \
                repeatmasker(genome_file, repeatmakser_lib, repeatmasker_dir, thread=num_threads)
            if genome_anno_result:
                click.echo("\nFinished whole genome TE annotation by RepeatMasker\n")

        #else:
            #click.echo(
                #"Less than 90% of the query sequences processed, TE Trimmer can't perform whole genome TE annotation")
    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"Genome TE annotation error.\n")
            f.write(tb_content + '\n\n')
            click.echo("Final genome annotation can't be performed. This won't affect the TE consensus library.")

    end_time = datetime.now()
    duration = end_time - start_time

    # Remove microseconds from the duration
    duration_without_microseconds = timedelta(days=duration.days, seconds=duration.seconds)

    # if not debug:
    # Remove all single files when all the sequences are processed
    # shutil.rmtree(single_file_dir)

    print(f"\nTE Trimmer finished at {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    print(f"TE Trimmer runtime was {duration_without_microseconds}")


# The following is necessary to make the script executable, i.e., python myscript.py.
if __name__ == '__main__':
    main()
