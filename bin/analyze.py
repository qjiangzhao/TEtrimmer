# Standard library imports
import os
import time
import traceback
import click
import pandas as pd

# Local imports
from functions import blast, remove_files_with_start_pattern, check_bed_uniqueness, \
    extract_fasta, handle_sequence_low_copy, handle_sequence_skipped, update_low_copy_cons_file, prcyan, prgre
# from Class_bed_filter import BEDFile
import bedfilter
from boundarycrop import find_boundary_and_crop
from TEaid import check_self_alignment
from MSAcluster import clean_and_cluster_MSA
from orfdomain import prepare_pfam_database, PlotPfam


# Define a function to check the progress file, which will be used to continue analysis if program exited prematurely
def check_progress_file(progress_file_path):
    # Read the progress file into a pandas DataFrame
    df = pd.read_csv(progress_file_path)

    # Calculate skipped and low copy element number
    skipped_count = df[df['status'].str.strip().str.lower() == 'skipped'].shape[0]
    low_copy_count = df[df['low_copy'].astype(str).str.strip().str.lower() == 'true'].shape[0]
    unknown_n = df[df['reclassified_type'].str.contains('Unknown', na=False)].shape[0]
    classifid_n = df[df['reclassified_type'].str.contains('/', na=False)].shape[0]

    # Calculate classified proportion
    if classifid_n != 0 and unknown_n != 0:
        classified_pro = classifid_n / (unknown_n + classifid_n)
    else:
        classified_pro = 0  # Set a default value (or any other value deemed appropriate)

    # Get unique 'input_name' values
    local_completed_sequences = df['input_name'].unique().tolist()

    return local_completed_sequences, skipped_count, low_copy_count, classified_pro


# Print iterations progress
def printProgressBar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█', printEnd="\r", final = False):
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


#####################################################################################################
# Code block: Define analyze_sequence function
#####################################################################################################

def analyze_sequence_helper(params):
    return analyze_sequence(*params)


def analyze_sequence(seq_obj, genome_file, MSA_dir, min_blast_len, min_seq_num, max_msa_lines,
                     top_mas_lines, max_cluster_num, cons_thr, ext_thr, ex_step, classification_dir,
                     max_extension, gap_thr, gap_nul_thr, crop_end_thr, crop_end_win, crop_end_gap_thr,
                     crop_end_gap_win, start_patterns, end_patterns, output_dir, pfam_dir, mini_orf,
                     single_fasta_n, hmm, check_extension_win, debug, progress_file,
                     classify_unknown, classify_all, final_con_file, final_unknown_con_file,
                     final_classified_con_file, low_copy_dir, fast_mode, error_files, plot_skip, skipped_dir,
                     plot_query, engine):
    #####################################################################################################
    # Code block: Set different elongation number for different elements and do BLAST search
    #####################################################################################################
    try:
        # Get query fasta file path
        seq_name = seq_obj.get_seq_name()
        seq_type = seq_obj.get_old_TE_type()
        seq_file = seq_obj.get_input_fasta()  # Return complete file path

        # Since DNA element are significantly shorter than LTR and LINE elements, adjust default parameters
        if "DNA" in seq_type:
            ex_step = 500
            max_extension = 7000
            min_blast_len = 150
            crop_end_gap_win = 100
            check_extension_win = 50

        # The average length of SINE elements is around 500 bp, adjust default parameters
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

        # run BLAST search for each FASTA file and return a BED file absolute path
        bed_out_file_dup, blast_hits_count, blast_out_file = blast(seq_file, genome_file, MSA_dir,
                                                                   min_length=min_blast_len, task="blastn",
                                                                   seq_obj=seq_obj, search_type=engine)

    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"Error while running blast for sequence: {seq_name}\n")
            f.write(tb_content + '\n\n')
        prcyan(f"Error while running blast for sequence: {seq_name}. Main Error: {str(e)}. \n"
               f"Trace back content: {tb_content}\n")
        return

    #####################################################################################################
    # Code block: Perform ORF and PFAM prediction for input sequences
    #####################################################################################################

    try:
        input_orf_pfam_obj = PlotPfam(seq_file, MSA_dir, pfam_database_dir=pfam_dir, mini_orf=mini_orf,
                                      after_tetrimmer=False)
        input_orf_domain_plot = None

        # "run_getorf()" function will return 'True' if any ORF was detected. Otherwise, it will return 'False'.
        if input_orf_pfam_obj.run_getorf():
            # "run_pfam_scan()" will return 'True' if any PFAM domains were found. Otherwise, it will return 'False'.
            pfam_scan_result = input_orf_pfam_obj.run_pfam_scan()
            input_orf_domain_plot = input_orf_pfam_obj.orf_domain_plot()

    except Exception as e:
        with open(error_files, "a") as f:
            tb_content = traceback.format_exc()
            f.write(f"Error when doing ORF and PFAM prediction for input sequence {seq_name}\n")
            f.write(tb_content + '\n\n')
        #prcyan(f"Error while performing input sequence ORF and PFAM predictions: {seq_name}. Main Error: {str(e)}. \n"
               #f"Trace back content: {tb_content}\n")
        pass

    #####################################################################################################
    # Code block: Check BLAST hit number
    #####################################################################################################

    try:
        # Check if BLAST hit number is exactly 0. If so, skip this sequence
        if blast_hits_count == 0:
            click.echo(f"\n{seq_name} is skipped due to blast hit number is 0\n")
            handle_sequence_skipped(seq_obj, progress_file, debug, MSA_dir, classification_dir)
            return

        # Check if BLAST hit number is smaller than "min_seq_num"; do not include "min_seq_num"
        elif blast_hits_count != 0 and blast_hits_count < min_seq_num:
            check_low_copy, blast_full_length_n, found_match, TE_aid_plot = check_self_alignment(
                seq_obj, seq_file, MSA_dir, genome_file, blast_hits_count, blast_out_file, plot_skip=plot_skip)

            if check_low_copy is True:
                # Update terminal repeat and BLAST full length number, remove low-copy intermediate files
                handle_sequence_low_copy(seq_obj, progress_file, debug, MSA_dir,
                                         classification_dir, found_match=found_match,
                                         blast_full_length_n=blast_full_length_n, te_aid_plot=TE_aid_plot,
                                         orf_plot=input_orf_domain_plot, low_copy_dir=low_copy_dir)

                # Integrate low-copy element sequences into consensus file
                update_low_copy_cons_file(seq_obj, final_con_file, final_unknown_con_file,
                                          final_classified_con_file, low_copy_dir, TE_aid_plot)
            else:
                click.echo(f"\n{seq_name} was skipped because the BLAST hit number is smaller than {min_seq_num} "
                           f"and check_low_copy is {check_low_copy}.\n")

                # handle_sequence_skipped will update skipped status
                handle_sequence_skipped(seq_obj, progress_file, debug, MSA_dir, classification_dir,
                                        plot_skip=plot_skip, te_aid_plot=TE_aid_plot, skip_proof_dir=skipped_dir,
                                        orf_plot=input_orf_domain_plot)

            return  # if BLAST hit number is smaller than 10, code will execute next FASTA file

    except Exception as e:
        # Add sequence to skip if check low-copy module returns errors
        handle_sequence_skipped(seq_obj, progress_file, debug, MSA_dir, classification_dir)
        with open(error_files, "a") as f:
            # Return the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"\nError while checking low-copy status for sequence: {seq_name}\n")
            f.write(tb_content + '\n\n')
        prcyan(f"\nError while checking low-copy status for sequence: {seq_name}. Error: {str(e)}\n")
        prgre("\nThe low-copy TE check module is an optional additional analysis. You may ignore this error, as it "
              "will not affect the final result significantly.\n")
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
        # top_mas_lines has to be equal to or smaller than max_mas_lines.
        bed_out_filter_file = bedfilter.process_lines(bed_out_file, MSA_dir, threshold=max_msa_lines,
                                                      top_longest_lines_count=top_mas_lines)

        # Extract FASTA from bed_out_filter_file
        # Return fasta_out_flank_file absolute path of FASTA file
        # It is imperative to group the sequences in the first round of MSA; extend for both ends is 0.
        fasta_out_flank_file, bed_out_flank_file = extract_fasta(
            bed_out_filter_file, genome_file, MSA_dir, left_ex=0, right_ex=0, nameonly=True)

        # Return 'False' if cluster number is 0. Return 'True' if divergent column number is smaller than 100.
        # Otherwise, return the subset BED and alignment files.
        cluster_MSA_result = clean_and_cluster_MSA(fasta_out_flank_file, bed_out_filter_file, MSA_dir,
                                                   clean_column_threshold=0.02,
                                                   min_length_num=min_seq_num, cluster_num=max_cluster_num,
                                                   cluster_col_thr=100, fast_mode=fast_mode)
    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"Error while grouping MSA: {seq_name}\n")
            f.write(tb_content + '\n\n')
        prcyan(f"\nError while grouping MSA for sequence: {seq_name}. Error: {str(e)}\n")
        prcyan('\n' + tb_content + '\n')
        return
    
    #####################################################################################################
    # Code block: Perform find_boundary_and_crop on clustered MSA if necessary
    #####################################################################################################

    try:
        # cluster_false means too few sequences were found in clusters from MSA (all_cluster_size < 10); TE Trimmer will skip this sequence.
        if cluster_MSA_result is False:
            check_low_copy, blast_full_length_n, found_match, TE_aid_plot = check_self_alignment(
                seq_obj, seq_file, MSA_dir, genome_file, blast_hits_count, blast_out_file, plot_skip=plot_skip)

            if check_low_copy is True:
                # Update terminal repeat and BLAST full-length number, remove low-copy intermediate files
                handle_sequence_low_copy(seq_obj, progress_file, debug, MSA_dir,
                                         classification_dir, found_match=found_match,
                                         blast_full_length_n=blast_full_length_n, te_aid_plot=TE_aid_plot,
                                         orf_plot=input_orf_domain_plot, low_copy_dir=low_copy_dir)
                update_low_copy_cons_file(seq_obj, final_con_file, final_unknown_con_file,
                                          final_classified_con_file, low_copy_dir, TE_aid_plot)

            else:
                click.echo(
                    f"\n{seq_name} was skipped because the sequence number in each cluster was smaller "
                    f"than {min_seq_num} and check_low_copy is {check_low_copy}.\n")
                handle_sequence_skipped(seq_obj, progress_file, debug, MSA_dir, classification_dir,
                                        plot_skip=plot_skip, te_aid_plot=TE_aid_plot, skip_proof_dir=skipped_dir,
                                        orf_plot=input_orf_domain_plot)
            return
        else:
            cluster_bed_files_list = cluster_MSA_result

            # cluster_pattern_alignment_list has the same index as cluster_bed_files_list
            all_inner_skipped = True
            for i in range(len(cluster_bed_files_list)):
                try:
                    find_boundary_result = find_boundary_and_crop(
                        cluster_bed_files_list[i], genome_file, MSA_dir, pfam_dir, seq_obj,
                        hmm, classify_all, classify_unknown, error_files, plot_query, cons_threshold=cons_thr,
                        ext_threshold=ext_thr, ex_step_size=ex_step, max_extension=max_extension,
                        gap_threshold=gap_thr, gap_nul_thr=gap_nul_thr,
                        crop_end_thr=crop_end_thr, crop_end_win=crop_end_win,
                        crop_end_gap_thr=crop_end_gap_thr, crop_end_gap_win=crop_end_gap_win,
                        start_patterns=start_patterns, end_patterns=end_patterns,
                        mini_orf=mini_orf, define_boundary_win=check_extension_win,
                        fast_mode=fast_mode, engine=engine, input_orf_pfam=input_orf_domain_plot, debug=debug)
                except Exception as e:
                    return

                if not find_boundary_result:
                    continue
                elif find_boundary_result:
                    all_inner_skipped = False

            # Check the flag after the loop. If all inner clusters were skipped, write the progress file.
            if all_inner_skipped:
                check_low_copy, blast_full_length_n, found_match, TE_aid_plot = check_self_alignment(
                    seq_obj, seq_file, MSA_dir, genome_file, blast_hits_count, blast_out_file, plot_skip=plot_skip)

                if check_low_copy is True:
                    # Update terminal repeat and BLAST full-length number, remove low-copy intermediate files.
                    handle_sequence_low_copy(seq_obj, progress_file, debug, MSA_dir,
                                             classification_dir, found_match=found_match,
                                             blast_full_length_n=blast_full_length_n, te_aid_plot=TE_aid_plot,
                                             orf_plot=input_orf_domain_plot, low_copy_dir=low_copy_dir)
                    update_low_copy_cons_file(seq_obj, final_con_file, final_unknown_con_file,
                                              final_classified_con_file, low_copy_dir, TE_aid_plot)
                else:
                    handle_sequence_skipped(seq_obj, progress_file, debug, MSA_dir, classification_dir,
                                            plot_skip=plot_skip, te_aid_plot=TE_aid_plot, skip_proof_dir=skipped_dir,
                                            orf_plot=input_orf_domain_plot)
                    click.echo(
                        f"\n{seq_name} was skipped because sequence is too short and check_low_copy is {check_low_copy}.\n")
                return

    except Exception as e:
        with open(error_files, "a") as f:
            # Return the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"\nError during boundary finding and cropping for sequence: {seq_name}\n")
            f.write(tb_content + '\n\n')
        prcyan(f"\nError during boundary finding and cropping for sequence: {seq_name}. Error: {str(e)}\n")
        prcyan(tb_content + '\n')
        return

    # After all processing is done, change status to 'process' and write the file name to the progress file
    seq_obj.update_status("processed", progress_file)

    # If analysis of this sequence has been completed, remove all files contain sequence name
    if not debug:
        remove_files_with_start_pattern(MSA_dir, seq_name)

    # Read and count sequences from Finished_sequence_name.txt
    completed_sequence, skipped_count, low_copy_count, classified_pro = check_progress_file(progress_file)

    # Calculate the total count
    processed_count = len(completed_sequence)

    # Calculate number of sequences not processed by TE Trimmer
    rest_sequence = single_fasta_n - processed_count

    printProgressBar(processed_count, single_fasta_n, prefix='Progress:', suffix='Complete', length=50)


def create_dir(continue_analysis, hmm, pfam_dir, output_dir, input_file, genome_file, plot_skip):

    if not os.path.isfile(input_file):
        prcyan(f"The FASTA file {input_file} does not exist. Please check the input file path!")
        raise FileNotFoundError
    input_file = os.path.abspath(input_file)  # get absolute path for input

    # check if genome file exist
    if not os.path.isfile(genome_file):
        prcyan(f"The genome FASTA file {genome_file} does not exist. Please check the genome file path!")
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
        current_time = time.strftime("%Y%m%d_%H%M%S")
        new_output_dir = os.path.join(output_dir, f"TETrimmer_output_{current_time}")
        os.makedirs(new_output_dir, exist_ok=True)
        output_dir = new_output_dir
        prgre(f"\nThe given output directory is not empty. Results will be stored into folder: \n"
              f"{output_dir}\n")

        # Stop the whole program if the output directory is not empty and
        # raise Exception

    # Create a new folder for single FASTA sequences
    single_file_dir = os.path.join(output_dir, "Single_fasta_files")
    os.makedirs(single_file_dir, exist_ok=True)

    # Create a new folder for MSA
    MSA_dir = os.path.join(output_dir, "Multiple_sequence_alignment")
    os.makedirs(MSA_dir, exist_ok=True)

    # Create a new folder for classifications
    classification_dir = os.path.join(output_dir, "Classification")
    os.makedirs(classification_dir, exist_ok=True)

    # Create a new folder for HMM files
    if hmm:
        hmm_dir = os.path.join(output_dir, "HMM_files")
        os.makedirs(hmm_dir, exist_ok=True)
    else:
        hmm_dir = ''

    # Define proof_annotation folder path
    proof_annotation_dir = os.path.join(output_dir, "TE_Trimmer_for_proof_annotation")
    os.makedirs(proof_annotation_dir, exist_ok=True)

    # Define skipped folder if required
    if plot_skip:
        skipped_dir = os.path.join(proof_annotation_dir, "Skipped_TE")
        os.makedirs(skipped_dir, exist_ok=True)
    else:
        skipped_dir = None

    # Define low-copy folder
    low_copy_dir = os.path.join(proof_annotation_dir, "Low_copy_TE")
    os.makedirs(low_copy_dir, exist_ok=True)

    # Define annotation evaluation folders
    perfect_proof = os.path.join(proof_annotation_dir, "Annotations_perfect")
    good_proof = os.path.join(proof_annotation_dir, "Annotations_good")
    intermediate_proof = os.path.join(proof_annotation_dir, "Annotations_check_recommended")
    need_check_proof = os.path.join(proof_annotation_dir, "Annotations_check_required")
    os.makedirs(perfect_proof, exist_ok=True)
    os.makedirs(good_proof, exist_ok=True)
    os.makedirs(intermediate_proof, exist_ok=True)
    os.makedirs(need_check_proof, exist_ok=True)

    # Define the progress_file, finished sequence IDs will be stored here
    progress_file = os.path.join(output_dir, "summary.txt")

    # Check and create progress_file if it does not exist yet
    if not os.path.exists(progress_file):
        with open(progress_file, 'a') as f:
            f.write("input_name,consensus_name,blast_hit_n,cons_MSA_seq_n,cons_full_blast_n,input_length,cons_length,"
                    "input_TE_type,reclassified_type,terminal_repeat,low_copy,evaluation,status\n")

    # Define error files to store non-mandatory function errors including RepeatClassified classification,
    # RepeatMasker classification, PFAM scanning, MUSCLE alignment
    error_files = os.path.join(MSA_dir, "error_file.txt")

    # If a PFAM database was not provided, create PFAM database in the TE Trimmer software folder, so the
    # database can be downloaded here.
    # If the pfam_dir was provided but the database cannot be found there, TE Trimmer will download the PFAM database
    # into the provided directory and generate the index file.
    if pfam_dir is None:
        pfam_dir = os.path.join(os.path.dirname(bin_py_path), "pfam_database")
        os.makedirs(pfam_dir, exist_ok=True)
    try:
        os.makedirs(pfam_dir, exist_ok=True)
        if_pfam = prepare_pfam_database(pfam_dir)

        if not if_pfam:  # Check if if_pfam is 'False'
            raise Exception

    except Exception as e:
        with open(error_files, "a") as f:
            # Return the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"PFAM database building error\n")
            f.write(tb_content + '\n\n')

        prgre(f"Note: Cannot download PFAM database from internet, please use a local PFAM database.\n"
              f"For example: --pfam_dir <your_PFAM_directory>\n"
              f"Your PFAM directory should contain: \n"
              f"Pfam-A.hmm\n"
              f"Pfam-A.hmm.h3f\n"
              f"Pfam-A.hmm.h3m\n"
              f"Pfam-A.hmm.dat\n"
              f"Pfam-A.hmm.h3i\n"
              f"Pfam-A.hmm.h3p\n\n"
              f"PFAM predictions will be used to determine the direction of TEs. The database is therefore mandatory.\n\n"
              f"You can download <Pfam-A.hmm.gz> and <Pfam-A.hmm.dat.gz> from "
              f"https://www.ebi.ac.uk/interpro/download/pfam/\n"
              f"Afterwards, do: \n"
              f"gzip -d Pfam-A.hmm.gz\n"
              f"gzip -d Pfam-A.hmm.dat.gz\n"
              f"hmmpress Pfam-A.hmm\n\n")
        return

    # Define consensus files. Temporary (temp) files will be used for the final RepeatMasker classification.
    final_con_file = os.path.join(output_dir, "TE_Trimmer_consensus.fasta")
    final_unknown_con_file = os.path.join(classification_dir, "temp_TE_Trimmer_unknown_consensus.fasta")
    final_classified_con_file = os.path.join(classification_dir, "temp_TE_Trimmer_classified_consensus.fasta")

    return bin_py_path, output_dir, single_file_dir, MSA_dir, classification_dir, hmm_dir, proof_annotation_dir, low_copy_dir, perfect_proof, \
    good_proof, intermediate_proof, need_check_proof, progress_file, pfam_dir, final_con_file, final_unknown_con_file, final_classified_con_file, \
    error_files, input_file, genome_file, skipped_dir

