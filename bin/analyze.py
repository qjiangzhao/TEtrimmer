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


# Define a function to check progress file, which will be used for continue analysis
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
        classified_pro = 0  # Set a default value (or any other value you deem appropriate)

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
                     single_fasta_n, hmm, hmm_dir, check_extension_win, debug, progress_file,
                     classify_unknown, classify_all, final_con_file, final_con_file_no_low_copy, final_unknown_con_file,
                     final_classified_con_file, low_copy_dir, fast_mode, error_files, plot_skip, skipped_dir,
                     plot_query, engine, proof_annotation_dir):
    #####################################################################################################
    # Code block: Set different elongation number for different elements and do blast
    #####################################################################################################

    try:
        # Get query fasta file path
        seq_name = seq_obj.get_seq_name()
        seq_type = seq_obj.get_old_TE_type()
        seq_file = seq_obj.get_input_fasta()  # Return full file path

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
    # Code block: Do ORF and PFAM prediction for input sequence
    #####################################################################################################

    try:
        input_orf_pfam_obj = PlotPfam(seq_file, MSA_dir, pfam_database_dir=pfam_dir, mini_orf=mini_orf,
                                      after_tetrimmer=False)
        input_orf_domain_plot = None

        # "run_getorf()" function will return True when any ORF are detected. Otherwise, it will be False
        if input_orf_pfam_obj.run_getorf():
            # "run_pfam_scan()" will return True when pfam domains are found, otherwise it will return False
            pfam_scan_result = input_orf_pfam_obj.run_pfam_scan()
            input_orf_domain_plot = input_orf_pfam_obj.orf_domain_plot()

    except Exception as e:
        with open(error_files, "a") as f:
            tb_content = traceback.format_exc()
            f.write(f"Error when doing ORF and PFAM prediction for input sequence {seq_name}\n")
            f.write(tb_content + '\n\n')
        #prcyan(f"Error while doing input sequence ORF and PFAM prediction: {seq_name}. Main Error: {str(e)}. \n"
               #f"Trace back content: {tb_content}\n")

    #####################################################################################################
    # Code block: Check blast hits number
    #####################################################################################################

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
                                         blast_full_length_n=blast_full_length_n, te_aid_plot=TE_aid_plot,
                                         orf_plot=input_orf_domain_plot, low_copy_dir=low_copy_dir)

                # Integrate low copy element sequence into consensus file
                update_low_copy_cons_file(seq_obj, final_con_file, final_unknown_con_file,
                                          final_classified_con_file, low_copy_dir, TE_aid_plot)
            else:
                click.echo(f"\n{seq_name} is skipped due to blast hit number is smaller than {min_seq_num} "
                           f"and check_low_copy is {check_low_copy}\n")

                # handle_sequence_skipped will update skipped status
                handle_sequence_skipped(seq_obj, progress_file, debug, MSA_dir, classification_dir,
                                        plot_skip=plot_skip, te_aid_plot=TE_aid_plot, skip_proof_dir=skipped_dir,
                                        orf_plot=input_orf_domain_plot)

            return  # when blast hit number is smaller than 10, code will execute next fasta file

    except Exception as e:
        # Add sequence to skip when check low copy module gets errors
        handle_sequence_skipped(seq_obj, progress_file, debug, MSA_dir, classification_dir)
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"\nError while checking low copy for sequence: {seq_name}\n")
            f.write(tb_content + '\n\n')
        prcyan(f"\nError while checking low copy for sequence: {seq_name}. Error: {str(e)}\n")
        prgre("\nLow copy TE check module is not very important, you can ignore this error message, which won't affect "
              "your final result too much.\n")
        return

    #####################################################################################################
    # Code block: Separate MSA based on sequence relatedness
    #####################################################################################################

    try:
        # remove duplicated lines
        bed_out_file = check_bed_uniqueness(MSA_dir, bed_out_file_dup)

        # test bed_out_file line number and extract top longest lines
        # for process_lines() function. threshold represent the maximum number to keep for MSA
        # top_longest_lines_count means the number of sequences with top length
        # for example if threshold = 100, top_longest_lines_count = 50, then 50 sequences will be
        # randomly chose from the rest of sequences
        # top_mas_lines has to be equal or smaller than max_mas_lines
        bed_out_filter_file = bedfilter.process_lines(bed_out_file, MSA_dir, threshold=max_msa_lines,
                                                      top_longest_lines_count=top_mas_lines)

        # extract fast from bed_out_filter_file
        # return fasta_out_flank_file absolute path
        # because have to group MSA the first round extend for left and right side are both 0
        fasta_out_flank_file, bed_out_flank_file = extract_fasta(
            bed_out_filter_file, genome_file, MSA_dir, left_ex=0, right_ex=0, nameonly=True)

        # Return False when cluster number is 0. Return True when divergent column number is smaller than 100
        # Otherwise it will return the subset bed and alignment file
        cluster_MSA_result = clean_and_cluster_MSA(fasta_out_flank_file, bed_out_filter_file, MSA_dir,
                                                   clean_column_threshold=0.02,
                                                   min_length_num=min_seq_num, cluster_num=max_cluster_num,
                                                   cluster_col_thr=100, fast_mode=fast_mode)
    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"Error when group MSA: {seq_name}\n")
            f.write(tb_content + '\n\n')
        prcyan(f"\nError when group MSA for sequence: {seq_name}. Error: {str(e)}\n")
        prcyan('\n' + tb_content + '\n')
        return
    
    #####################################################################################################
    # Code block: perform find_boundary_and_crop on clustered MSA when necessary
    #####################################################################################################

    try:
        # cluster false means no sufficient cluster (all cluster size < 10), TETrimmer will skip this sequence.
        if cluster_MSA_result is False:
            check_low_copy, blast_full_length_n, found_match, TE_aid_plot = check_self_alignment(
                seq_obj, seq_file, MSA_dir, genome_file, blast_hits_count, blast_out_file, plot_skip=plot_skip)

            if check_low_copy is True:
                # Update terminal repeat and blast full length number, remove low copy intermediate files
                handle_sequence_low_copy(seq_obj, progress_file, debug, MSA_dir,
                                         classification_dir, found_match=found_match,
                                         blast_full_length_n=blast_full_length_n, te_aid_plot=TE_aid_plot,
                                         orf_plot=input_orf_domain_plot, low_copy_dir=low_copy_dir)
                update_low_copy_cons_file(seq_obj, final_con_file, final_unknown_con_file,
                                          final_classified_con_file, low_copy_dir, TE_aid_plot)

            else:
                click.echo(
                    f"\n{seq_name} is skipped due to sequence number in each cluster is smaller "
                    f"than {min_seq_num} and check_low_copy is {check_low_copy}\n")
                handle_sequence_skipped(seq_obj, progress_file, debug, MSA_dir, classification_dir,
                                        plot_skip=plot_skip, te_aid_plot=TE_aid_plot, skip_proof_dir=skipped_dir,
                                        orf_plot=input_orf_domain_plot)
            return
        else:
            cluster_bed_files_list = cluster_MSA_result

            # cluster_pattern_alignment_list have the same index with cluster_bed_files_list
            all_inner_skipped = True
            for i in range(len(cluster_bed_files_list)):
                try:
                    find_boundary_result = find_boundary_and_crop(
                        cluster_bed_files_list[i], genome_file, MSA_dir, pfam_dir, seq_obj,
                        hmm, classify_all, classify_unknown, error_files, plot_query, classification_dir,
                        final_con_file, final_con_file_no_low_copy, proof_annotation_dir, hmm_dir, cons_threshold=cons_thr,
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

            # Check the flag after the loop. If all inner clusters were skipped, write the progress file
            if all_inner_skipped:
                check_low_copy, blast_full_length_n, found_match, TE_aid_plot = check_self_alignment(
                    seq_obj, seq_file, MSA_dir, genome_file, blast_hits_count, blast_out_file, plot_skip=plot_skip)

                if check_low_copy is True:
                    # Update terminal repeat and blast full length number, remove low copy intermediate files
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
                        f"\n{seq_name} is skipped due to sequence is too short and check_low_copy is {check_low_copy}\n")
                return

    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"\nError during boundary finding and cropping for sequence: {seq_name}\n")
            f.write(tb_content + '\n\n')
        prcyan(f"\nError during boundary finding and cropping for sequence: {seq_name}. Error: {str(e)}\n")
        prcyan(tb_content + '\n')
        return

    # After all processing is done, change status to process and write the name of the file to the progress file
    seq_obj.update_status("processed", progress_file)

    # If all this sequence is finished remove all files contain this name
    if not debug:
        remove_files_with_start_pattern(MSA_dir, seq_name)

    # Read and count sequences from Finished_sequence_name.txt
    completed_sequence, skipped_count, low_copy_count, classified_pro = check_progress_file(progress_file)

    # Calculate the total count
    processed_count = len(completed_sequence)

    # Calculate sequences number that hasn't been processed by TETrimmer
    rest_sequence = single_fasta_n - processed_count

    printProgressBar(processed_count, single_fasta_n, prefix='Progress:', suffix='Complete', length=50)


def create_dir(continue_analysis, hmm, pfam_dir, output_dir, input_file, genome_file, plot_skip):

    if not os.path.isfile(input_file):
        prcyan(f"The fasta file {input_file} does not exist. Please check your input file path!")
        raise FileNotFoundError
    input_file = os.path.abspath(input_file)  # get input absolute path

    # check if genome file exist
    if not os.path.isfile(genome_file):
        prcyan(f"The genome fasta file {genome_file} does not exist. Please check your genome file path!")
        raise FileNotFoundError
    genome_file = os.path.abspath(genome_file)  # get genome absolute path

    # bin_py_path contains all classes and bash code
    # so.path.abspath(__file__) will return the current executable python file
    bin_py_path = os.path.dirname(os.path.abspath(__file__))

    # check if output path exist otherwise create it
    os.makedirs(output_dir, exist_ok=True)
    output_dir = os.path.abspath(output_dir)  # get absolute path

    # Check if output directory is empty when --continue_analysis is False
    if os.listdir(output_dir) and not continue_analysis:

        """
        prcyan(f"\nWARNING: The output directory {output_dir} is not empty. Please empty your output directory or "
               f"choose another empty directory.")
        prgre("\nNOTE: TETrimmer can create output directory when it is not exist.")
        """
        # When the current folder isn't empty, create a new folder with current time
        current_time = time.strftime("%Y%m%d_%H%M%S")
        new_output_dir = os.path.join(output_dir, f"TETrimmer_output_{current_time}")
        os.makedirs(new_output_dir, exist_ok=True)
        output_dir = new_output_dir
        prgre(f"\nThe given output directory is not empty. Results are stored into folder: \n"
              f"{output_dir}\n")

        # Stop the whole program when the output directory is not empty
        #  raise Exception

    # Make a new folder for single fasta sequence
    single_file_dir = os.path.join(output_dir, "Single_fasta_files")
    os.makedirs(single_file_dir, exist_ok=True)

    # Make a new folder for MSA
    MSA_dir = os.path.join(output_dir, "Multiple_sequence_alignment")
    os.makedirs(MSA_dir, exist_ok=True)

    # Make a new folder for classification
    classification_dir = os.path.join(output_dir, "Classification_and_deduplication")
    os.makedirs(classification_dir, exist_ok=True)

    # make a new folder for HMM file
    if hmm:
        hmm_dir = os.path.join(output_dir, "HMM_files")
        os.makedirs(hmm_dir, exist_ok=True)
    else:
        hmm_dir = ''

    # Define proof_annotation folder path
    proof_annotation_dir = os.path.join(output_dir, "TETrimmer_for_proof_annotation")
    os.makedirs(proof_annotation_dir, exist_ok=True)

    # Define clustered proof annotation folder inside proof_annotation_dir
    cluster_proof_anno_dir = os.path.join(proof_annotation_dir, "Clustered_proof_annotation")
    os.makedirs(cluster_proof_anno_dir, exist_ok=True)

    # Define skipped folder if it is required
    if plot_skip:
        skipped_dir = os.path.join(proof_annotation_dir, "Skipped_TE")
        os.makedirs(skipped_dir, exist_ok=True)
    else:
        skipped_dir = None

    # Define low copy folder
    low_copy_dir = os.path.join(proof_annotation_dir, "Low_copy_TE")
    os.makedirs(low_copy_dir, exist_ok=True)

    # Define proof annotation evaluation folders
    perfect_proof = os.path.join(proof_annotation_dir, "Perfect_annotation")
    good_proof = os.path.join(proof_annotation_dir, "Good_annotation")
    intermediate_proof = os.path.join(proof_annotation_dir, "Recommend_check_annotation")
    need_check_proof = os.path.join(proof_annotation_dir, "Need_check_annotation")
    os.makedirs(perfect_proof, exist_ok=True)
    os.makedirs(good_proof, exist_ok=True)
    os.makedirs(intermediate_proof, exist_ok=True)
    os.makedirs(need_check_proof, exist_ok=True)

    # Define the progress_file, finished seq IDs will be stored here
    progress_file = os.path.join(output_dir, "summary.txt")

    # Check and create progress_file if it doesn't exist
    if not os.path.exists(progress_file):
        with open(progress_file, 'a') as f:
            f.write("input_name,consensus_name,blast_hit_n,cons_MSA_seq_n,cons_full_blast_n,input_length,cons_length,"
                    "input_TE_type,reclassified_type,terminal_repeat,low_copy,evaluation,status\n")

    # Define error files to store not mandatory function errors including RepeatClassified classification,
    # RepeatMasker classification, PFAM scanning, muscle alignment
    error_files = os.path.join(MSA_dir, "error_file.txt")

    # If pfam database isn't provided, create pfam database at TETrimmer software folder,
    # the database will be downloaded into there.
    # If the pfam_dir is given, but the database can't be found there, TETrimmer will download pfam database there
    # and generate index file
    if pfam_dir is None:
        pfam_dir = os.path.join(os.path.dirname(bin_py_path), "pfam_database")
    try:
        os.makedirs(pfam_dir, exist_ok=True)
        if_pfam = prepare_pfam_database(pfam_dir)

        if not if_pfam:  # Check if if_pfam is False
            raise Exception

    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"PFAM database building error\n")
            f.write(tb_content + '\n\n')

        prgre(f"Note: Can't download PFAM database from internet, please use your own PFAM database\n"
              f"For example: --pfam_dir <your_PFAM_directory>\n"
              f"Your PFAM directory should contain: \n"
              f"Pfam-A.hmm\n"
              f"Pfam-A.hmm.h3f\n"
              f"Pfam-A.hmm.h3m\n"
              f"Pfam-A.hmm.dat\n"
              f"Pfam-A.hmm.h3i\n"
              f"Pfam-A.hmm.h3p\n\n"
              f"PFAM prediction will be used to determine the direction of TEs, it is mandatory to have it\n\n"
              f"You can download <Pfam-A.hmm.gz> and <Pfam-A.hmm.dat.gz> from "
              f"https://www.ebi.ac.uk/interpro/download/pfam/\n"
              f"After: \n"
              f"gzip -d Pfam-A.hmm.gz\n"
              f"gzip -d Pfam-A.hmm.dat.gz\n"
              f"hmmpress Pfam-A.hmm\n\n")
        return

    # Define consensus files. temp files will be used for the final RepeatMasker classification
    final_con_file = os.path.join(output_dir, "TETrimmer_consensus.fasta")
    final_unknown_con_file = os.path.join(classification_dir, "temp_TETrimmer_unknown_consensus.fasta")
    final_classified_con_file = os.path.join(classification_dir, "temp_TETrimmer_classified_consensus.fasta")

    # Define consensus file without low copy
    final_con_file_no_low_copy = os.path.join(classification_dir, "TETrimmer_consensus_no_low_copy.fasta")

    return bin_py_path, output_dir, single_file_dir, MSA_dir, classification_dir, hmm_dir, proof_annotation_dir, \
        low_copy_dir, perfect_proof, good_proof, intermediate_proof, need_check_proof, progress_file, pfam_dir, \
        final_con_file, final_con_file_no_low_copy, final_unknown_con_file, final_classified_con_file, \
        error_files, input_file, genome_file, skipped_dir, cluster_proof_anno_dir

