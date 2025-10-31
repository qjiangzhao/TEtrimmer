import logging
import os
import shutil
import traceback

import pandas as pd
from Bio import AlignIO, SeqIO
from matplotlib.pyplot import title

import cialign
from boundaryclass import CropEnd, CropEndByGap, DefineBoundary, GenomeBlastCoverage

# Local imports
from functions import (
    align_sequences,
    check_start_and_end_patterns,
    check_terminal_repeat,
    classify_single,
    con_generater,
    con_generater_no_file,
    concatenate_alignments,
    define_crop_end_simi_thr,
    dotplot,
    extract_fasta,
    find_poly_a_end_position,
    generate_hmm_from_msa,
    is_LTR,
    merge_pdfs,
    modify_fasta_headers,
    remove_files_with_start_pattern,
    remove_gaps_with_similarity_check,
    reverse_complement_seq_file,
    scale_single_page_pdf,
    select_start_to_end,
    select_start_end_and_join,
    select_window_columns,
    pairwise_seqs_align
)
from MSAcluster import CleanAndSelectColumn, clean_and_cluster_MSA, process_msa

from orfdomain import PlotPfam, determine_sequence_direction
from TEaid import TEAid


def reset_bed(input_file, output_dir):
    # This function will prepare bed file ready for the left and right side extension

    # calculate alignment length in column 6
    df = pd.read_csv(input_file, sep='\t', header=None)

    """
    # Define a threshold for outlier removal
    threshold = 0.5

    # Identify the top 10 sequence lengths
    top_10_lengths = df.nlargest(10, 6)[6]

    # Calculate the mean and standard deviation of the top 10 lengths
    mean_top_10_lengths = top_10_lengths.mean()
    std_top_10_lengths = top_10_lengths.std()

    # Identify values outside the threshold
    df['difference_from_mean'] = df[6] - mean_top_10_lengths

    # ~ is used to negate the condition, thus keeping rows where the condition is not satisfied.
    filtered_df = df[~((abs(df['difference_from_mean']) > threshold * std_top_10_lengths)
                       & (df['difference_from_mean'] < 0))]
    """
    filtered_df = df

    # Conditionally update values in column 1 and column 2 for right extension
    # 1 adjusted to avoid error message "Error: malformed BED entry at line 91. Start was greater than end"
    reset_right_filtered_df = filtered_df.copy()
    reset_right_filtered_df.loc[reset_right_filtered_df[5] == '+', 1] = (
        reset_right_filtered_df.loc[reset_right_filtered_df[5] == '+', 2] - 1
    )
    reset_right_filtered_df.loc[reset_right_filtered_df[5] == '-', 2] = (
        reset_right_filtered_df.loc[reset_right_filtered_df[5] == '-', 1] + 1
    )

    # Conditionally update values in column 1 and column 2 for left extension
    reset_left_filtered_df = filtered_df.copy()
    reset_left_filtered_df.loc[reset_left_filtered_df[5] == '+', 2] = (
        reset_left_filtered_df.loc[reset_left_filtered_df[5] == '+', 1] + 1
    )
    reset_left_filtered_df.loc[reset_left_filtered_df[5] == '-', 1] = (
        reset_left_filtered_df.loc[reset_left_filtered_df[5] == '-', 2] - 1
    )

    # write new BED file
    reset_right_long_bed = os.path.join(
        output_dir, f'{os.path.basename(input_file)}_reset_right_long.bed'
    )
    reset_left_long_bed = os.path.join(
        output_dir, f'{os.path.basename(input_file)}_reset_left_long.bed'
    )
    reset_right_filtered_df.to_csv(
        reset_right_long_bed, sep='\t', index=False, header=None
    )
    reset_left_filtered_df.to_csv(
        reset_left_long_bed, sep='\t', index=False, header=None
    )

    return reset_left_long_bed, reset_right_long_bed


def extend_end(
    max_extension,
    ex_step_size,
    end,
    input_file,
    genome_file,
    output_dir,
    crop_end_gap_thr,
    crop_end_gap_win,
    ext_threshold,
    ext_buffer,
    define_boundary_win,
    bed_dic,
):
    # end: left or right extension
    if_ex = True
    ex_total = 0

    # The following code calculates the majority of the alignment length and finds the longest alignment for extension
    # while ignoring the shorter ones.
    # long_bed function won't do length selection
    # It generates the proper bed file for left and right extension.
    reset_left_long_bed, reset_right_long_bed = reset_bed(input_file, output_dir)

    # adjust is the overlap region length between adjacent extended MSAs.
    adjust = ex_step_size / 10
    ite = 1
    reset_left_long_bed_copy = reset_left_long_bed
    reset_right_long_bed_copy = reset_right_long_bed
    while if_ex and (ex_total < max_extension):
        # BEDtools will make sure that the extension will not excess the maximum length of the chromosome
        # track extend total
        ex_total += (ex_step_size - adjust)

        if end == 'left':
            left_ex = ex_total
            right_ex = ex_step_size - ex_total + adjust
            bed_fasta, bed_out_flank_file = extract_fasta(
                reset_left_long_bed,
                genome_file,
                output_dir,
                left_ex,
                right_ex,
                nameonly=True,
            )
        if end == 'right':
            left_ex = ex_step_size - ex_total + adjust
            right_ex = ex_total
            bed_fasta, bed_out_flank_file = extract_fasta(
                reset_right_long_bed,
                genome_file,
                output_dir,
                left_ex,
                right_ex,
                nameonly=True,
            )

        # align_sequences() will return extended MSA in absolute file
        bed_fasta_mafft_with_gap = align_sequences(bed_fasta, output_dir)

        if not os.path.isfile(bed_fasta_mafft_with_gap):
            logging.error(
                f'{input_file} encountered a problem during the MAFFT extension step.'
            )
            return False

        # Remove nucleotide whose proportion is smaller than threshold
        bed_fasta_mafft_with_gap_column_clean_object = CleanAndSelectColumn(
            bed_fasta_mafft_with_gap, threshold=0.08
        )
        bed_fasta_mafft_with_gap_column_clean = (
            bed_fasta_mafft_with_gap_column_clean_object.clean_column(output_dir)
        )

        # Remove gaps with similarity check
        bed_fasta_mafft = remove_gaps_with_similarity_check(
            bed_fasta_mafft_with_gap_column_clean,
            output_dir,
            gap_threshold=0.6,
            simi_check_gap_thre=0.4,
            similarity_threshold=0.7,
            min_nucleotide=5,
            return_map=False,
        )

        # bed_fasta_mafft will be false if the MSA column number is smaller than 50 after removing the gaps
        if not bed_fasta_mafft:
            break
        bed_fasta_mafft_object = CropEndByGap(
            bed_fasta_mafft,
            gap_threshold=crop_end_gap_thr,
            window_size=crop_end_gap_win,
        )
        # if the crop_end_by_gap removed more than 90% of the aligned sequences, this sequence ID will be stored in
        # large_crop_ids
        large_crop_ids, remain_n, remain_ids = bed_fasta_mafft_object.find_large_crops()
        cropped_alignment = bed_fasta_mafft_object.crop_alignment()
        bed_fasta_mafft_cop_end_gap = bed_fasta_mafft_object.write_to_file(
            output_dir, cropped_alignment
        )

        # ext_threshold refer to option --ext_thr
        # max_X is float number, means the proportion of X
        bed_boundary = DefineBoundary(
            bed_fasta_mafft_cop_end_gap,
            threshold=ext_threshold,
            check_window=define_boundary_win,
            max_X=0.3,
            if_con_generater=False,
            extension_buffer=ext_buffer,
        )
        # Read bed_out_flank_file
        bed_out_flank_file_df = pd.read_csv(bed_out_flank_file, sep='\t', header=None)

        # Update bed_dic
        for i in bed_dic:
            id_list = large_crop_ids + remain_ids
            if i in id_list:
                matching_rows = bed_out_flank_file_df[bed_out_flank_file_df[3] == i]
                if end == 'left':
                    if bed_dic[i][2] == '+':
                        bed_dic[i][0] = matching_rows.iloc[0, 1]
                    elif bed_dic[i][2] == '-':
                        bed_dic[i][1] = matching_rows.iloc[0, 2]
                if end == 'right':
                    if bed_dic[i][2] == '+':
                        bed_dic[i][1] = matching_rows.iloc[0, 2]
                    elif bed_dic[i][2] == '-':
                        bed_dic[i][0] = matching_rows.iloc[0, 1]

        # If the remaining sequences that need further extension are fewer than 5, stop the extension
        # bed_boundary.if_continue will be false if the extension is enough
        if not bed_boundary.if_continue or remain_n < 5:
            break
        elif bed_boundary.if_continue:
            # Subset BED file to only keep sequences that need further extension
            if end == 'left':
                if_ex = bed_boundary.left_ext
                df = pd.read_csv(reset_left_long_bed, sep='\t', header=None)

                # Filter out rows where the fourth column is in large_crop_ids
                df_filtered = df[~df[3].isin(large_crop_ids)]
                reset_left_long_bed = reset_left_long_bed_copy + f'_{ite}'

                # Write the filtered DataFrame back to reset_left_long_bed
                df_filtered.to_csv(
                    reset_left_long_bed, sep='\t', index=False, header=None
                )

            if end == 'right':
                if_ex = bed_boundary.right_ext
                df = pd.read_csv(reset_right_long_bed, sep='\t', header=None)

                # Filter out rows where the fourth column is in large_crop_ids
                df_filtered = df[~df[3].isin(large_crop_ids)]
                reset_right_long_bed = reset_right_long_bed_copy + f'_{ite}'

                # Write the filtered DataFrame into reset_left_long_bed
                df_filtered.to_csv(
                    reset_right_long_bed, sep='\t', index=False, header=None
                )
            ite += 1
    # if_ex will be true if the extension stopped due to reach to the maximum extension limitation
    return bed_dic, ex_total, if_ex


# fianl_MSA aims to provide the boundaries
def final_MSA(
    bed_final_MSA,
    MSA_seq_n,
    genome_file,
    output_dir,
    gap_nul_thr,
    gap_threshold,
    ext_threshold,
    define_boundary_win,
    crop_end_gap_thr,
    crop_end_gap_win,
    crop_end_thr,
    crop_end_win,
    seq_obj,
    poly_patterns, 
    poly_len,
    blast_database_path
):

    bed_fasta, bed_out_flank_file = extract_fasta(
        bed_final_MSA, genome_file, output_dir, 0, 0, nameonly=True
    )

    # Align extended TE sequences
    # align_sequences() will return extended MSA in absolute file
    bed_fasta_mafft_with_gap = align_sequences(bed_fasta, output_dir)

    if not os.path.isfile(bed_fasta_mafft_with_gap):
        logging.error(
            f'{bed_final_MSA} encountered a problem during the final MAFFT extension step.'
        )
        return False

    # Remove nucleotides whose proportion is smaller than threshold
    bed_fasta_mafft_with_gap_column_clean_object = CleanAndSelectColumn(
        bed_fasta_mafft_with_gap, threshold=0.01
    )
    bed_fasta_mafft_with_gap_column_clean = bed_fasta_mafft_with_gap_column_clean_object.clean_column(output_dir)

    # Remove gaps by similarity check
    bed_fasta_mafft_gap_sim, column_mapping_initial = remove_gaps_with_similarity_check(
        bed_fasta_mafft_with_gap_column_clean, output_dir, return_map=True,
    )

    # bed_fasta_mafft_boundary_crop_for_select is used to see the regions beyond the defined boundaries
    bed_fasta_mafft_boundary_crop_for_select = bed_fasta_mafft_gap_sim

    #####################################################################################################
    # Code block: Generate consensus sequences
    #####################################################################################################

    # Generate consensus sequence before end cropping to perform genome blast coverage calculation
    # this function will write the consensus sequence to a file
    bed_fasta_mafft_gap_sim_con_02 = con_generater(bed_fasta_mafft_gap_sim, output_dir, threshold=0.2, ambiguous='N')

    # Crop end before terminal check, use loose threshold to avoid overhead cropping
    bed_fasta_mafft_gap_sim_cp_object = CropEnd(bed_fasta_mafft_gap_sim, threshold=20, window_size=40)
    bed_fasta_mafft_gap_sim_cp_object.crop_alignment()
    bed_fasta_mafft_gap_sim_cp = bed_fasta_mafft_gap_sim_cp_object.write_to_file(output_dir)

    # Generate consensus sequence, use low threshold to reduce number of N's in the consensus sequence
    # Use more stringent threshold for poly A identification when terminal is not found.
    bed_fasta_mafft_gap_sim_cp_con_08 = con_generater(
        bed_fasta_mafft_gap_sim_cp, output_dir, threshold=0.8, ambiguous='N'
    )

    # read sequence and convert to plain upper-case string
    bed_fasta_mafft_gap_sim_cp_con_08_record = SeqIO.read(bed_fasta_mafft_gap_sim_cp_con_08, "fasta")
    bed_fasta_mafft_gap_sim_cp_con_08_seq = str(bed_fasta_mafft_gap_sim_cp_con_08_record.seq).upper()
    bed_fasta_mafft_gap_sim_cp_con_08_len = len(bed_fasta_mafft_gap_sim_cp_con_08_seq)

    # helper: clip coordinates so they stay inside the sequence
    def clip_con(start, stop, seq, seq_len):
        """
        Return seq[a:b] but keep 0 ≤ a ≤ b ≤ len(seq)
        """
        start = max(0, start)
        stop = min(seq_len, stop)
        return seq[start:stop]

    # Generate consensus sequence, use low threshold to reduce number of N's in the consensus sequence
    # for LTR/TIR checking
    bed_fasta_mafft_gap_sim_cp_con_045 = con_generater(
        bed_fasta_mafft_gap_sim_cp, output_dir, threshold=0.45, ambiguous='N'
    )

    #####################################################################################################
    # Code block: Define boundary by the genome blast coverage
    #####################################################################################################

    # Calculate the genome blast coverage list
    bed_fasta_mafft_gap_sim_con_coverage_obj = GenomeBlastCoverage(
        bed_fasta_mafft_gap_sim_con_02,
        blast_database_path,
        output_dir
    )

    bed_fasta_mafft_gap_sim_con_coverage_obj.calculate_blast_coverage()
    start_posit_cov, end_posit_cov = bed_fasta_mafft_gap_sim_con_coverage_obj.find_boundary_blast_coverage()

    def print_horizontal(lst, cols=8, max_width=20):
        """
        Print list horizontally in 'cols' columns.
        Each line is prefixed with the index of its first element.
        Elements are truncated to max_width chars for readability.
        """
        n = len(lst)
        if n == 0:
            print("[empty]")
            return

        # width for the left label like "[123]:"
        label_w = len(str(n - 1)) + 3

        def short(x):
            s = str(x)
            return s if len(s) <= max_width else s[:max_width - 1] + "…"

        for start in range(0, n, cols):
            end = min(start + cols, n)
            row = [short(lst[i]) for i in range(start, end)]
            print(f"[{start}]".ljust(label_w), "  ".join(row))

    # usage
    #print_horizontal(bed_fasta_mafft_gap_sim_con_coverage_obj.coverage_list, cols=20, max_width=25)

    #####################################################################################################
    # Code block: TE boundary definition by MSA conservation
    #####################################################################################################

    # For LINE elements bed_fasta_mafft_boundary_crop will be changed. Copy alignment to another variable
    if seq_obj.old_TE_type.startswith("LINE"):
        # For highly divergent regions, more gaps can be found. According to this feature, remove
        # high-divergence regions. This function is very useful for dealing with LINE elements.
        # don't crop the right ends of the MSA when poly A was found for the LINE element

        cropped_MSA_by_gap_object = CropEndByGap(
            bed_fasta_mafft_gap_sim_cp,
            gap_threshold=crop_end_gap_thr,
            window_size=crop_end_gap_win
        )
        cropped_alignment_line = cropped_MSA_by_gap_object.crop_alignment()
        bed_fasta_mafft_gap_sim_cp = cropped_MSA_by_gap_object.write_to_file(
            output_dir, cropped_alignment_line
        )

    # Gaps are removed again after processing with CropEnd
    # column_mapping build relation between bed_fasta_mafft_gap_sim_cp and cropped_alignment_output_file_g
    # bed_fasta_mafft_gap_sim_cp shares the same column coordinate like bed_fasta_mafft_boundary_crop_for_select
    cropped_MSA_output_file, cropped_MSA_output_file_gs, column_mapping_MSA = crop_end_and_remove_gap(
        bed_fasta_mafft_gap_sim_cp,
        output_dir,
        crop_end_threshold=crop_end_thr,
        window_size=crop_end_win,
        gap_threshold=0.8
    )

    # This means the sequence length after gap removal is shorter than 50
    if not cropped_MSA_output_file_gs:
        return False

    # Use DefineBoundary to define start position.
    cropped_boundary_obj_MSA = DefineBoundary(
        cropped_MSA_output_file_gs,
        threshold=0.8,
        check_window=4,
        max_X=0
    )

    start_posit_MSA = column_mapping_MSA[cropped_boundary_obj_MSA.start_post]
    end_posit_MSA = column_mapping_MSA[cropped_boundary_obj_MSA.end_post]

    #####################################################################################################
    # Code block: TE right boundary definition by poly-A
    #####################################################################################################

    # Check if poly A can be found from the sequence and return the last A position
    # poly_a will be None if not found
    if seq_obj.old_TE_type.startswith("LINE") or seq_obj.old_TE_type.startswith("SINE"):

        # Poly_a will be None if poly A is not found
        poly_a = find_poly_a_end_position(bed_fasta_mafft_gap_sim_cp_con_08,
                                          poly_patterns=poly_patterns,
                                          min_length=poly_len)

        # right 100-nt window of the poly-A tail
        poly_a_beyond_right_window = clip_con(poly_a, poly_a + 100,
                                              bed_fasta_mafft_gap_sim_cp_con_08_seq,
                                              bed_fasta_mafft_gap_sim_cp_con_08_len
                                              )
        poly_a_beyond_right_window_Ns = poly_a_beyond_right_window.count('N')

        # When the N number is too less, mean the beyond region is still conserved.
        if poly_a_beyond_right_window < 50:
            poly_a = None
    else:
        poly_a = None

    #####################################################################################################
    # Code block: TE boundary definition by LTR TIR
    #####################################################################################################

    # Check terminal repeats
    # Define the output folder to store the temporary BLAST database
    check_terminal_repeat_output = f'{bed_fasta_mafft_gap_sim_cp_con_045}_tem'

    LTR_boundary, TIR_boundary, found_match_crop = check_terminal_repeat(
        bed_fasta_mafft_gap_sim_cp_con_045, check_terminal_repeat_output
    )

    left_posit_repeat = None
    right_posit_repeat = None

    # Check terminal repeats
    if LTR_boundary is not None and len(LTR_boundary) == 2:
        left_posit_repeat, right_posit_repeat = LTR_boundary
    elif TIR_boundary is not None and len(TIR_boundary) == 2:
        left_posit_repeat, right_posit_repeat = TIR_boundary

    # Check the MSA conservation beyond the defined boundaries
    if left_posit_repeat is not None and right_posit_repeat is not None:

        # left 100-nt window
        beyond_left_window = clip_con(left_posit_repeat - 100, left_posit_repeat,
                               bed_fasta_mafft_gap_sim_cp_con_08_seq,
                               bed_fasta_mafft_gap_sim_cp_con_08_len)
        beyond_left_Ns = beyond_left_window.count('N')

        # right 100-nt window
        beyond_right_window = clip_con(right_posit_repeat, right_posit_repeat + 100,
                                       bed_fasta_mafft_gap_sim_cp_con_08_seq,
                                       bed_fasta_mafft_gap_sim_cp_con_08_len
                                       )
        beyond_right_Ns = beyond_right_window.count('N')

        # When the N number is too less, mean the beyond region is still conserved.
        if beyond_left_Ns < 50 and beyond_right_Ns < 50:
            left_posit_repeat = None
            right_posit_repeat = None

    #####################################################################################################
    # Code block: Evaluate boundary result and choose the final_start and final_end
    #####################################################################################################

    if left_posit_repeat is not None and right_posit_repeat is not None:
        final_start = left_posit_repeat
        final_end = right_posit_repeat

    elif poly_a is None:
        if MSA_seq_n >= 25:
            if seq_obj.old_TE_type.startswith("LINE"):
                final_start = start_posit_cov
            else:
                final_start = start_posit_MSA
            final_end = end_posit_MSA
        else:
            final_start = start_posit_cov
            final_end = end_posit_cov

    elif poly_a is not None:
        final_end = poly_a

        if MSA_seq_n >= 25:
            if seq_obj.old_TE_type.startswith("LINE"):
                final_start = start_posit_cov
            else:
                final_start = start_posit_MSA
        else:
            final_start = start_posit_cov
    else:
        final_start = start_posit_MSA
        final_end = end_posit_MSA

    #####################################################################################################
    # Code block: Slice the MSA
    #####################################################################################################
    # Extract MSA based on the defined start and end position
    bed_fasta_mafft_gap_sim_final_start_to_end, not_use_a, not_use_b = select_start_to_end(
        cropped_MSA_output_file,  # This file share the same column index with bed_fasta_mafft_gap_sim
        output_dir,
        final_start,
        final_end
    )



    '''
    Y: mean the difference isn't larger than 20. 
    N: mean the difference is larger than 20.
    Only compare with LTR/TIR boundaries when it is found. Otherwise, compare coverage_blast and MSA_conservation 
    boundaries.
    
    When comparing to LTR/TIR both the left and right boundaries should reach the requirement (difference < 20)
    When comparing the coverage_blast and MSA_conservation, compare the left and right boundaries separately.
    Decide the boundary based on the MSA sequence number.
    
    For the poly-tail information, compare when it is not None. Compare the right boundary from coverage_blast and 
    MSA_conservation, if the distance difference with the poly-tail is smaller than 30 use the poly-tail as the right
    boundary, for this don't consider the MSA sequence number.
    
    LTR/TIR     Coverage_blast      MSA_conservation    Poly-tail       others       Final choise       full_len_b      Evaluation
    Found       Y                   Y                   Not consider                 LTR/TIR             >=2            Good
    Found       Y                   N                   Not consider                 LTR/TIR             >=2            Intermediate
    Found       N                   Y                   Not consider                 LTR/TIR             >=2            Intermediate
    Found       N                   N                   Not consider                 LTR/TIR                            Need_check   
    Not_found   Y                   Y                                   MSA>=25      MSA_conservation    >=2            Intermediate   
    Not_found   Y                   Y                                   MSA<25       Coverage_blast      >=2            Intermediate
    Not_found   N                   N                                   MSA>=25      MSA_conservation                   Need_check   
    Not_found   N                   N                                   MSA<25       Coverage_blast                     Need_check 
    '''

    return (
        bed_out_flank_file,  # the bed file after MSA extension
        bed_fasta_mafft_gap_sim_final_start_to_end,  # the final MSA file after boundary definition
        bed_fasta_mafft_boundary_crop_for_select, # the raw MSA before boundary definition
        found_match_crop,  # can be False, LTR, or TIR
        bed_fasta_mafft_with_gap,
        final_start,
        final_end,
        bed_fasta_mafft_gap_sim_cp_con_08_len
    )


def crop_end_and_remove_gap(
    input_file, output_dir, crop_end_threshold=0.8, window_size=20, gap_threshold=0.8, crop_r=True
):
    # window_size means the checked nucleotide number each time
    # threshold means the nucleotide proportion in window_size must greater than INT
    CropEnd_thres = crop_end_threshold * window_size
    cropped_MSA = CropEnd(input_file, threshold=CropEnd_thres, window_size=window_size, crop_r=crop_r)
    cropped_MSA.crop_alignment()

    # write_to_file() function will return absolute cropped alignment file
    cropped_MSA_output_file = cropped_MSA.write_to_file(output_dir)

    # Remove gaps again after the CropEnd step
    # "column_mapping" is a dictionary; the key is gap removed MSA index, the value is the corresponding original
    # MSA index
    cropped_MSA_output_file_gs, column_mapping = remove_gaps_with_similarity_check(
        cropped_MSA_output_file, output_dir, return_map=True
    )

    return cropped_MSA_output_file, cropped_MSA_output_file_gs, column_mapping


def find_boundary_and_crop(
    bed_file,
    genome_file,
    output_dir,
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
    more_extension_dir,
    hmm_dir,
    poly_patterns, 
    poly_len,    
    cons_threshold=0.8,
    ext_threshold=0.7,
    ex_step_size=1000,
    max_extension=7000,
    extension_buffer=300,
    gap_threshold=0.4,
    gap_nul_thr=0.7,
    crop_end_thr=0.8,
    crop_end_win=40,
    crop_end_gap_thr=0.1,
    crop_end_gap_win=150,
    ltr_start_patterns=None,
    ltr_end_patterns=None,
    helitron_start_patterns=None,
    helitron_end_patterns=None,
    mini_orf=200,
    define_boundary_win=150,
    fast_mode=False,
    engine='blast',
    input_orf_pfam=False,
    debug=False,
    cluster_msa=None,
    perfect_seq_num=30,
    blast_database_path=None,
    mmseqs_database_dir=None,
):

    """
    :param bed_file: str, BED file directory
    :param genome_file: str, directory containing the genome FASTA file
    :param output_dir: str, output directory
    :param pfam_dir: str, PFAM database directory
    :param cons_threshold: float (0 to 1), threshold used for final consensus sequence generation. Default: 0.8
    :param ext_threshold: float (0 to 1), threshold used to define the extent of sequence extension. A higher ratio
    represents a more stringent extension threshold (shorter extensions), while a smaller ratio represents a loosened
    threshold, yielding longer extensions on both sides of the sequence. Default: 0.7
    :param ex_step_size: int, number of nucleotides to bee added to both sides in each extension cycle. Default: 1000
    :param max_extension: int, the maximum number of nucleotides to be added on both ends for extension. Default: 7000
    :param crop_end_thr: int, if ambiguous sides exceed <crop_end_thr> in the extension window, delete the window. Default: 16
    :param crop_end_win: int, window size used for the crop-end process. Default: 20
    :param gap_threshold: float, columns with a greater gap proportion than <gap_threshold>, and if the highest
    nucleotide proportion in this column is less than <gap_nul_thr>, will be removed. Default: 0.4
    :param gap_nul_thr: float, set nucleotide proportion threshold to decide if the column should be removed. Default: 0.07
    :param crop_end_gap_thr: float, set gap threshold to crop end by gap. Default: 0.05
    :param crop_end_gap_win: int, set window size used to crop end by gap. Default: 300.
    :param ltr_start_patterns: str, patterns to check for start points. Default: None
    :param ltr_end_patterns: str, patterns to check for end points. Default: None
    :param min_orf: int, set minimum ORF length for ORF prediction. Default: 200

    """
    #####################################################################################################
    # Code block: Define consensus sequence object
    #####################################################################################################
    try:
        seq_name = seq_obj.name
        seq_file = seq_obj.get_input_fasta()  # Return full file path

        # Define unique sequence names
        # Because seq_obj.create_consi_obj() is generated after the unique name definition, len(seq_obj.consi_obj_list)
        # is the appended number for the unique name.
        consi_n = len(seq_obj.consi_obj_list)

        if consi_n > 0:
            uniq_seq_name = f'{seq_name}_{consi_n:02}'
        else:
            uniq_seq_name = seq_name

        # Create consensus object
        consi_obj = seq_obj.create_consi_obj(uniq_seq_name)

    except Exception as e:
        with open(error_files, 'a') as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f'\nCreate consensus sequence object failed for {seq_name} with error:\n{e}')
            f.write('\n' + tb_content + '\n\n')
        logging.error(f'\nCreate consensus sequence object failed for {seq_name} with error:\n{e}\n{tb_content}\n')
        raise Exception from e


    #####################################################################################################
    # Code block: Get extended MSA
    #####################################################################################################

    # Check if this is a LINE element, if so decrease the ext_threshold number,
    # because LINE have high divergence at the 5' end
    if seq_obj.old_TE_type.startswith("LINE"):
        ext_threshold = ext_threshold - 0.1

    msa_loop_n = 1

    while msa_loop_n <= 2:
        # TEtrimmer will do the second round clustering and extension. The MSA has been extensively extended during
        # the first time of loop, reduce the max_extension to 2000 for the second round, which can prevent to generate
        # very long sequence especially for tandem repeats.
        if msa_loop_n == 2:
            max_extension = 2 * ex_step_size

        try:
            # Read BED file and build dictionary to store sequence position
            df = pd.read_csv(bed_file, sep='\t', header=None)

            # df_row shows the MSA sequence number, which will be used in the final_MSA function
            df_row_n = df.shape[0]

            # Initialize an empty dictionary to store the desired key-value pairs
            bed_dict = {}

            # Iterate through the DataFrame to populate the dictionary
            # Here is an example of bed file content to help understanding the code
            """
            Scaffold_116	341861	341862	1	94.573	+	6468	rnd_1_family_188
            Scaffold_116	412105	412106	3	93.762	+	6492	rnd_1_family_188
            Scaffold_119	245388	245389	4	94.512	-	6469	rnd_1_family_188
            Scaffold_119	220368	220369	5	94.465	-	6468	rnd_1_family_188
            Scaffold_1	238299	238300	8	94.373	-	6469	rnd_1_family_188
            Scaffold_1	219863	219864	9	93.803	-	6471	rnd_1_family_188
            """
            for index, row in df.iterrows():
                key = row.iloc[3]

                # Use sequence position and strand as dictionary values
                value_list = [row.iloc[1], row.iloc[2], row.iloc[5]]

                # Add the key-value pair to the dictionary
                bed_dict[key] = value_list

        except Exception as e:
            with open(error_files, 'a') as f:
                # Get the traceback content as a string
                tb_content = traceback.format_exc()
                f.write(
                    f'\nFind boundary bed_file reading failed for {seq_name} with error:\n{e}'
                )
                f.write('\n' + tb_content + '\n\n')
            logging.error(
                f"\nFind boundary bed_file reading failed for {seq_name} with error:\n{e}"
            )

            raise Exception from e

        #####################################################################################################
        # Code block: Define length of extensions on left and right sides of the MSA and define the boundary
        #####################################################################################################
        try:
            # Extend MSA left end
            # left_ext_to_max will be true if the extension stopped due to reach to the maximum extension limitation
            left_bed_dic, left_ex_total, left_ext_to_max = extend_end(
                max_extension,
                ex_step_size,
                'left',
                bed_file,
                genome_file,
                output_dir,
                crop_end_gap_thr,
                crop_end_gap_win,
                ext_threshold,
                extension_buffer,
                define_boundary_win,
                bed_dict,
            )

            # Extend MSA right end
            final_bed_dic, right_ex_total, right_ext_to_max = extend_end(
                max_extension,
                ex_step_size,
                'right',
                bed_file,
                genome_file,
                output_dir,
                crop_end_gap_thr,
                crop_end_gap_win,
                ext_threshold,
                extension_buffer,
                define_boundary_win,
                left_bed_dic,
            )
        except Exception as e:
            with open(error_files, 'a') as f:
                # Get the traceback content as a string
                tb_content = traceback.format_exc()
                f.write(f'\nMSA extension failed for {seq_name} with error:\n{e}')
                f.write('\n' + tb_content + '\n\n')
            logging.error(f'\nMSA extension failed for {seq_name} with error:\n{e}\n{tb_content}\n')
            raise Exception from e

        try:
            # Update BED file for final MSA
            for key, value in final_bed_dic.items():
                # Find rows where the fourth column matches the key and update position by dictionary value
                df.loc[df[3] == key, [1, 2]] = value[:2]

            # Define BED file name for final MSA
            bed_final_MSA = f'{bed_file}_fm_{msa_loop_n}.bed'
            df.to_csv(bed_final_MSA, sep='\t', index=False, header=False)

            # final_MSA return 'False' if when the start crop point is greater than the end crop point, which means the
            # sequence is too short. If the divergence is too high, CropEnd will return 'False'.

            final_msa_result = final_MSA(
                bed_final_MSA,
                df_row_n,
                genome_file,
                output_dir,
                gap_nul_thr,
                gap_threshold,
                ext_threshold,
                define_boundary_win,
                crop_end_gap_thr,
                crop_end_gap_win,
                crop_end_thr,
                crop_end_win,
                seq_obj,
                poly_patterns,
                poly_len,
                blast_database_path
            )

            # Check if final_msa_result returned 'False', indicating an error or specific condition
            if not final_msa_result:
                # Handle the error or specific condition
                return False
            else:
                # Unpack the returned values if the function executed normally
                (
                    bed_out_flank_file,
                    cropped_boundary_MSA,
                    bed_fasta_mafft_boundary_crop_for_select,
                    found_match_crop,
                    bed_fasta_mafft_with_gap,
                    final_start,
                    final_end,
                    MSA_length
                ) = final_msa_result

        except Exception as e:
            with open(error_files, 'a') as f:
                # Return the traceback content as a string
                tb_content = traceback.format_exc()
                f.write(f'Boundary definition error {seq_name}\n')
                f.write('\n' + tb_content + '\n\n')
            logging.error(f'\nBoundary definition failed for {seq_name} with error:\n{e}\n{tb_content}\n')
            raise Exception from e

        #####################################################################################################
        # Code block: Check the consistency of the final MSA
        #####################################################################################################
        # Plotting is done after PFAM predictions in case consensus/MSA are in the wrong direction.
        # Don't do second round MSA clustering after boundary definition when LTR or TIR is found.
        try:
            if msa_loop_n <= 1:
                final_msa_consistency = clean_and_cluster_MSA(
                    cropped_boundary_MSA,
                    bed_out_flank_file,
                    output_dir,
                    input_msa=cropped_boundary_MSA,
                    clean_column_threshold=0.01,
                    min_length_num=10,
                    cluster_num=2,
                    cluster_col_thr=500,
                )

                # final_msa_consistency will be false if no cluster available. In this case, use the original MSA
                # after boundary definition
                if not final_msa_consistency:
                    break
                else:
                    cluster_bed_files_list, fasta_out_flank_mafft_gap_rm = final_msa_consistency

                # If the clustered MSA contains only 1 element
                if len(cluster_bed_files_list) == 1:
                    df_new = pd.read_csv(cluster_bed_files_list[0], sep='\t', header=None)

                    # Check if only less than 2 sequence was eliminated
                    if len(df_new) >= len(df) - 2:
                        break
                    else:
                        bed_file = cluster_bed_files_list[0]
                else:
                    # Only use the top 1 cluster for further analysis
                    bed_file = cluster_bed_files_list[0]
        except Exception as e:
            with open(error_files, 'a') as f:
                # Get the traceback content as a string
                tb_content = traceback.format_exc()
                f.write(f'Boundary definition clustering error {seq_name}\n')
                f.write('\n' + tb_content + '\n\n')
            logging.error(
                f'\nBoundary definition clustering failed for {seq_name} with error:\n{e}\n{tb_content}\n')
            raise Exception from e

        msa_loop_n += 1

        # Use smaller extension step for the next round extension
        # ex_step_size = 700

    #####################################################################################################
    # Code block: Check if the final MSA contains too many instances of the ambiguous letter "N"
    #####################################################################################################
    try:
        initial_cons = con_generater_no_file(
            cropped_boundary_MSA, threshold=0.55, ambiguous='N'
        )

        # Calculate the proportion of 'N' in the sequence
        n_proportion = initial_cons.count('N') / len(initial_cons)

        # Check if the proportion is greater than 40%
        # If 40% of this consensus sequence is "N", stop analysis for this MSA
        if n_proportion > 0.4:
            return False

        #####################################################################################################
        # Code block: For LTR elements, check if the cropped MSA starts with given patterns like TGA ACA
        #####################################################################################################
        # If both start and end patterns are 'None', and found_match_crop is 'False' (no terminal repeat found),
        # skip this block, because terminal repeats can precisely define the start and end points of TEs.
        if seq_obj.old_TE_type.startswith("LTR"):  # Check if file name contains "LTR"

            if found_match_crop in ("LTR", "TIR"):
                sliding_win_check_for_start_end_pattern = False
            else:
                sliding_win_check_for_start_end_pattern = True

            # Generate consensus sequences
            # use bed_fasta_mafft_boundary_crop_for_select rather cropped_boundary_MSA
            ltr_consensus_seq = con_generater_no_file(
                bed_fasta_mafft_boundary_crop_for_select, threshold=0.7, ambiguous='N'
            )

            # Convert consensus_seq to list
            ltr_consensus_seq_list = list(ltr_consensus_seq)

            # start_matched and end_matched are boolen values
            # start_matched_pattern and end_matched_pattern are string
            (
                start_matched,
                end_matched,
                start_matched_pattern,
                end_matched_pattern,
                check_start,
                check_end
             ) = check_start_and_end_patterns(
                ltr_consensus_seq_list,
                start=final_start,
                end=final_end,
                start_patterns=ltr_start_patterns,
                end_patterns=ltr_end_patterns,
                check_helitron=False,
                sliding_win_check=sliding_win_check_for_start_end_pattern
            )

            # Update  start and end pattern for the consensus sequence object
            consi_obj.set_start_pattern_content(start_matched_pattern)
            consi_obj.set_end_pattern_content(end_matched_pattern)

            # If the new start or end positions are different from previous MSA, replace with the new
            # start or end position in the cropped_boundary object and generate the new MSA file
            if (
                final_start != check_start
                or final_end != check_end
            ):
                final_start = check_start
                final_end = check_end

    except Exception as e:
        with open(error_files, 'a') as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(
                f'\nStart and end pattern definition failed for {seq_name} with error \n{e}\n'
            )
            f.write('\n' + tb_content + '\n\n')
        logging.error(
            f"\nStart and end patterns for LTR elements like 'TGT' or 'ACA' definition failed for {seq_name} with error \n{e}"
        )
        logging.warning(
            '\nThis error will not affect the final result significantly, you can choose to ignore it. '
            "For traceback text, please refer to 'error_file.txt' in the 'Multiple_sequence_alignment' folder.\n"
        )

    #####################################################################################################
    # Code block: For Helitrons, check if the boundary start and end patterns
    #####################################################################################################
    try:

        if "helitron" in seq_obj.old_TE_type.lower():  # Check if file name contains "Helitron or helitron"

            # Generate consensus sequences
            # use bed_fasta_mafft_boundary_crop_for_select rather cropped_boundary_MSA
            helitron_consensus_seq = con_generater_no_file(
                bed_fasta_mafft_boundary_crop_for_select, threshold=0.7, ambiguous='N'
            )

            # Convert consensus_seq to list
            helitron_consensus_seq_list = list(helitron_consensus_seq)

            # start_matched and end_matched are boolen values
            # start_matched_pattern and end_matched_pattern are string
            (
                start_matched,
                end_matched,
                start_matched_pattern,
                end_matched_pattern,
                check_start,
                check_end
            ) = check_start_and_end_patterns(
                helitron_consensus_seq_list,
                start=final_start,
                end=final_end,
                start_patterns=helitron_start_patterns,
                end_patterns=helitron_end_patterns,
                check_helitron=True
            )

            # Update  start and end pattern for the consensus sequence object
            consi_obj.set_start_pattern_content(start_matched_pattern)
            consi_obj.set_end_pattern_content(end_matched_pattern)

            # If the new start or end positions are different from previous MSA, replace with the new
            # start or end position in the cropped_boundary object and generate the new MSA file
            if (
                    final_start != check_start
                    or final_end != check_end
            ):
                final_start = check_start
                final_end = check_end

    except Exception as e:
        with open(error_files, 'a') as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(
                f'\nStart and end pattern definition failed for {seq_name} with error \n{e}\n'
            )
            f.write('\n' + tb_content + '\n\n')
        logging.error(
            f"\nStart and end patterns for Helitron definition failed for {seq_name} with error \n{e}"
        )
        logging.warning(
            '\nThis error will not affect the final result significantly, you can choose to ignore it. '
            "For traceback text, please refer to 'error_file.txt' in the 'Multiple_sequence_alignment' folder.\n"
        )

    #####################################################################################################
    # Code block: Clean the final MSA
    #####################################################################################################
    try:
        extension_enough = True

        # If the defined final_start or final_end position too close to the edge of the
        # bed_fasta_mafft_boundary_crop_for_select MSA, label this as extension not enough
        if final_start <= 5 or final_end >= MSA_length - 5:
            extension_enough = False

        # Extract MSA based on the defined start and end position
        bed_fasta_mafft_gap_sim_final_start_to_end, not_use_a, not_use_b = select_start_to_end(
            bed_fasta_mafft_boundary_crop_for_select, output_dir, final_start, final_end
        )

        # Define the threshold use for cropping end by similarity to make sure the ends of the MSA are not totally removed
        # because the boundary of the MSA has been defined.
        start_pro_mean, end_pro_mean = define_crop_end_simi_thr(
            bed_fasta_mafft_gap_sim_final_start_to_end
        )

        if start_pro_mean < crop_end_win * 0.8:
            # Set crop_end_sim_thr to be smaller than pro_mean to avoid overcropping.
            start_pro_mean = int(start_pro_mean) - 3
        else:
            start_pro_mean = crop_end_win * 0.7

        # Because the divergence for the start and end parts of the MSA can be different, crop separately.
        if end_pro_mean < crop_end_win * 0.8:
            end_pro_mean = int(end_pro_mean) - 3
        else:
            end_pro_mean = crop_end_win * 0.7

        # Crop end by similarity to clean MSA for left side
        bed_fasta_mafft_gap_sim_final_start_to_end_left_obj = CropEnd(
            bed_fasta_mafft_gap_sim_final_start_to_end,
            threshold=start_pro_mean,
            window_size=40,
            crop_l=True,
            crop_r=False,
        )
        bed_fasta_mafft_gap_sim_final_start_to_end_left_obj.crop_alignment()
        bed_fasta_mafft_gap_sim_final_start_to_end_left = (
            bed_fasta_mafft_gap_sim_final_start_to_end_left_obj.write_to_file(output_dir)
        )

        # Crop end by similarity to clean MSA for right side
        bed_fasta_mafft_gap_sim_final_start_to_end_right_obj = CropEnd(
            bed_fasta_mafft_gap_sim_final_start_to_end_left,  # Crop right based on the result of left cropping
            threshold=end_pro_mean,
            window_size=40,
            crop_l=False,
            crop_r=True,
        )
        bed_fasta_mafft_gap_sim_final_start_to_end_right_obj.crop_alignment()
        bed_fasta_mafft_gap_sim_final_start_to_end_cp = bed_fasta_mafft_gap_sim_final_start_to_end_right_obj.write_to_file(
            output_dir)

        # Remove gaps by similarity again
        bed_fasta_mafft_gap_sim_final_start_to_end_cp_g, column_mapping = remove_gaps_with_similarity_check(
            bed_fasta_mafft_gap_sim_final_start_to_end_cp, output_dir, return_map=True
        )

        cropped_boundary_MSA = bed_fasta_mafft_gap_sim_final_start_to_end_cp_g

    except Exception as e:
        with open(error_files, 'a') as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(
                f'\nFinal MSA cleaning failed for {seq_name} with error:\n{e}'
            )
            f.write('\n' + tb_content + '\n\n')
        logging.error(
            f'\nFinal MSA cleaning failed for {seq_name} with error:\n{e}\n{tb_content}\n')
        raise Exception from e

    #####################################################################################################
    # Code block: Select MSA for outer TE boundaries TE-Aid plot
    #####################################################################################################
    try:
        # select_start_to_end will convert out_boundary_msa_start to 0 when it is negative
        # and convert out_boundary_msa_end to the MSA lengh when it is longer than the MSA length
        out_boundary_msa_start = final_start - 150

        out_boundary_msa_end = final_end + 150

        out_boundary_msa_for_teaid, out_boundary_msa_start_new, out_boundary_msa_end_new = select_start_to_end(
            bed_fasta_mafft_boundary_crop_for_select,
            output_dir,
            out_boundary_msa_start,
            out_boundary_msa_end
        )

        start_relate_to_out_boundary_msa_for_teaid = final_start - out_boundary_msa_start_new

        end_relate_to_out_boundary_msa_for_teaid = final_end - out_boundary_msa_start_new

    except Exception as e:
        with open(error_files, 'a') as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(
                f'\nGenerate out boundary MSA failed for {seq_name} with error:\n{e}'
            )
            f.write('\n' + tb_content + '\n\n')
        logging.error(
            f'\nGenerate out boundary MSA failed for {seq_name} with error:\n{e}\n{tb_content}\n')
        raise Exception from e

    #####################################################################################################
    # Code block: Generate MSA for CIAlign plot (entire MSA plotting)
    #####################################################################################################
    try:
        # Add 300 columns to the left of start point
        cropped_boundary_manual_MSA_left = select_window_columns(
            bed_fasta_mafft_boundary_crop_for_select,
            output_dir,
            final_start,
            'left',
            window_size=300,
        )

        # Add 300 columns to the right of the end point
        cropped_boundary_manual_MSA_right = select_window_columns(
            bed_fasta_mafft_boundary_crop_for_select,
            output_dir,
            final_end,
            'right',
            window_size=300,
        )

        # Concatenate MSA for manual curation
        # The concat_start and concat_end points correspond with cropped_boundary.start_post and cropped_boundary.end_post
        # Because the input file of "concatenate_alignment()" is an alignment object and not a file path, read "cropped_boundary_MSA"
        # as an alignment object.
        input_file_name = (
            str(os.path.basename(bed_fasta_mafft_boundary_crop_for_select))
            + '_proof_anno'
        )
        cropped_boundary_MSA_alignment = AlignIO.read(cropped_boundary_MSA, 'fasta')
        cropped_boundary_manual_MSA_concatenate, concat_start_man, concat_end_man = (
            concatenate_alignments(
                cropped_boundary_manual_MSA_left,
                cropped_boundary_MSA_alignment,
                cropped_boundary_manual_MSA_right,
                input_file_name=input_file_name,
                output_dir=output_dir,
            )
        )

        #####################################################################################################
        # Code block: Concatenate beginning and end part of MSA and plot
        #####################################################################################################

        # "sequence_len" represent the length of the final cropped MSA
        # Extract the beginning and end columns of the cropped MSA, then join them by "----------" (pseudo-gaps).
        # Return MSA length, which will be used for plotting
        cropped_boundary_plot_select_start_end_and_joint, sequence_len = select_start_end_and_join(
                cropped_boundary_MSA,
                output_dir,
            )

        # "column_mapping" is a dictionary, the key represents the cropped nucleotide position and the value represents the
        # original MSA nucleotide position.
        # Add 50 columns to the left of the start point
        cropped_boundary_plot_left = select_window_columns(
            bed_fasta_mafft_boundary_crop_for_select,
            output_dir,
            final_start,
            'left',
            window_size=50,
        )

        # Add 50 columns to the right of the end point
        cropped_boundary_plot_right = select_window_columns(
            bed_fasta_mafft_boundary_crop_for_select,
            output_dir,
            final_end,
            'right',
            window_size=50,
        )

        # Concatenate MSA
        # The concat_start and concat_end points correspond with cropped_boundary.start_post and cropped_boundary.end_post
        cropped_boundary_plot_concatenate, concat_start, concat_end = (
            concatenate_alignments(
                cropped_boundary_plot_left,
                cropped_boundary_plot_select_start_end_and_joint,
                cropped_boundary_plot_right,
                input_file_name=cropped_boundary_MSA,
                output_dir=output_dir,
            )
        )
    except Exception as e:
        with open(error_files, 'a') as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(
                f'\nConcatenate MSA for manual inspection plot failed for {seq_name} with error:\n{e}'
            )
            f.write('\n' + tb_content + '\n\n')
        logging.error(
            f'\nConcatenate MSA for manual inspection plot failed for {seq_name} with error:\n{e}\n{tb_content}\n')
        raise Exception from e

    #####################################################################################################
    # Code block: Predict ORFs and PFAM domains, determine sequence direction
    #####################################################################################################
    try:
        # Generate consensus sequence for ORF prediction
        # Use lower threshold to enable  more PFAM results
        orf_cons = con_generater(cropped_boundary_MSA, output_dir, threshold=0.5)

        # Predict ORF, scan for PFAM domains
        orf_domain_plot_object = PlotPfam(
            orf_cons,
            output_dir,
            pfam_database_dir=pfam_dir,
            mini_orf=mini_orf,
            after_tetrimmer=True,
        )

        # run_getorf() function will return 'True' if any ORF is detected. Otherwise, it will return 'False'.
        orf_domain_plot = None
        if_pfam_domain = False
        reverse_complement = False
        if orf_domain_plot_object.run_getorf():
            # run_pfam_scan() will return 'True' if any PFAM domain is found. Otherwise, it will return 'False'.
            if_pfam_domain, pfam_result_file = orf_domain_plot_object.run_pfam_scan()

            if if_pfam_domain:
                # determine_sequence_direction() will return 'True' if the direction of MSA and match are identical.
                # Otherwise, it will return 'False'.
                if determine_sequence_direction(pfam_result_file):
                    # If the direction matches, plot the ORF and PFAM matches
                    orf_domain_plot = orf_domain_plot_object.orf_domain_plot()
                else:
                    # If the direction is wrong, reverse-complement the consensus sequence and the corresponding MSA.
                    # The reverse-complemented file will overwrite the old file.
                    orf_cons = reverse_complement_seq_file(
                        input_file=orf_cons, output_file=orf_cons
                    )

                    # Reverse-complement MSA files
                    cropped_boundary_MSA = reverse_complement_seq_file(
                        input_file=cropped_boundary_MSA,
                        output_file=f'{cropped_boundary_MSA}_rc.fa',
                    )
                    # cropped_boundary_manual_MSA_concatenate will be used for the CIAlign plot
                    cropped_boundary_manual_MSA_concatenate = reverse_complement_seq_file(
                            input_file=cropped_boundary_manual_MSA_concatenate,
                            output_file=f'{cropped_boundary_manual_MSA_concatenate}_rc.fa',
                        )
                    # cropped_boundary_plot_concatenate will be used for the first MSA plot
                    cropped_boundary_plot_concatenate = reverse_complement_seq_file(
                        input_file=cropped_boundary_plot_concatenate,
                        output_file=f'{cropped_boundary_plot_concatenate}_rc.fa',
                    )
                    # bed_fasta_mafft_boundary_crop_for_select is used as raw for the manual curation
                    bed_fasta_mafft_boundary_crop_for_select = reverse_complement_seq_file(
                        input_file=bed_fasta_mafft_boundary_crop_for_select,
                        output_file=f'{bed_fasta_mafft_boundary_crop_for_select}_rc.fa',
                    )

                    # bed_fasta_mafft_with_gap is used as raw for the manual curation in the new version
                    # because this is before the gappy column cleaning and all the TSD information is kept
                    bed_fasta_mafft_with_gap = reverse_complement_seq_file(
                        input_file=bed_fasta_mafft_with_gap,
                        output_file=f'{bed_fasta_mafft_with_gap}_rc.fa',
                    )

                    # out_boundary_msa_for_teaid will be used for TEAid plot for the out boundary MSA
                    out_boundary_msa_for_teaid = reverse_complement_seq_file(
                        input_file=out_boundary_msa_for_teaid,
                        output_file=f'{out_boundary_msa_for_teaid}_rc.fa',
                    )

                    # Reverse complement input sequence, this will be used for dotplot
                    seq_file_reverse_c_path = os.path.join(
                        output_dir, f'{os.path.basename(seq_name)}.fasta_rc'
                    )
                    seq_file_reverse_c = reverse_complement_seq_file(
                        input_file=seq_file, output_file=seq_file_reverse_c_path
                    )

                    reverse_complement = True

                    # Define the new start and end points for cropped_boundary_manual_MSA_concatenate
                    # cropped_boundary_manual_MSA_concatenate will be used for the CIAlign plot
                    # Get MSA sequence length
                    cropped_boundary_manual_MSA_concatenate_align = AlignIO.read(
                        cropped_boundary_manual_MSA_concatenate, 'fasta'
                    )

                    # get_alignment_length() is a build-in function
                    cropped_boundary_manual_MSA_concatenate_length = cropped_boundary_manual_MSA_concatenate_align.get_alignment_length()

                    # Use intermediate number to get the new start and end number
                    concat_start_man_intermediate = (
                        cropped_boundary_manual_MSA_concatenate_length - concat_end_man - 1
                    )
                    concat_end_man_intermediate = (
                        cropped_boundary_manual_MSA_concatenate_length - concat_start_man - 1
                    )
                    concat_start_man = concat_start_man_intermediate
                    concat_end_man = concat_end_man_intermediate

                    # Define the new start and end points for cropped_boundary_plot_concatenate
                    # cropped_boundary_plot_concatenate will be used for the te_trimmer plot
                    # Get MSA sequence length
                    cropped_boundary_plot_concatenate_align = AlignIO.read(
                        cropped_boundary_plot_concatenate, 'fasta'
                    )
                    cropped_boundary_plot_concatenate_length = (
                        cropped_boundary_plot_concatenate_align.get_alignment_length()
                    )

                    # Use intermediate number to get the new start and end points
                    concat_start_intermediate = (
                        cropped_boundary_plot_concatenate_length - concat_end - 1
                    )
                    concat_end_intermediate = (
                        cropped_boundary_plot_concatenate_length - concat_start - 1
                    )
                    concat_start = concat_start_intermediate
                    concat_end = concat_end_intermediate

                    # Based on the new consensus sequence, predict ORFs and PFAM domains again, then plot
                    orf_domain_plot_object = PlotPfam(
                        orf_cons,
                        output_dir,
                        pfam_database_dir=pfam_dir,
                        mini_orf=mini_orf,
                    )
                    if orf_domain_plot_object.run_getorf():
                        # "run_pfam_scan() will return 'True' if any PFAM domains were found. Otherwise, it will return 'False'
                        if_pfam_domain, pfam_result_file = (
                            orf_domain_plot_object.run_pfam_scan()
                        )
                        if if_pfam_domain:
                            orf_domain_plot = orf_domain_plot_object.orf_domain_plot()

            else:
                # If only ORFs but no PFAM domain were found, plot ORFs
                orf_domain_plot = orf_domain_plot_object.orf_domain_plot()
    except Exception as e:
        with open(error_files, 'a') as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(
                f'\nConcatenate MSA for the maunual inspection plot failed for {seq_name} with error:\n{e}'
            )
            f.write('\n' + tb_content + '\n\n')
        logging.error(f'\nORF or PFAM prediction failed for {seq_name} with error:\n{e}\n{tb_content}\n')
        raise Exception from e


    #####################################################################################################
    # Code block: Plot multiple sequence alignment.
    #####################################################################################################
    # Plotting is done after PFAM predictions in case consensus/MSA are in the wrong direction.
    try:
        # Plot MSA, which can easily verify if the start and end crop points are correct
        MSA_plot = process_msa(
            cropped_boundary_plot_concatenate,
            output_dir,
            concat_start,
            concat_end,
            sequence_len,
        )

        #####################################################################################################
        # Code block: Plot entire MSA by CIAlign with start and end points
        #####################################################################################################

        # Define MSA plot output file
        cropped_boundary_manual_MSA_concatenate_plot = (
            f'{cropped_boundary_manual_MSA_concatenate}_plot.pdf'
        )

        # Draw the entire MSA
        cialign.drawMiniAlignment(
            cropped_boundary_manual_MSA_concatenate,
            cropped_boundary_manual_MSA_concatenate_plot,
            concat_start_man,
            concat_end_man,
        )
    except Exception as e:
        with open(error_files, 'a') as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f'\nMSA plot failed for {seq_name} with error:\n{e}')
            f.write('\n' + tb_content + '\n\n')
        logging.error(f'MSA plot failed {seq_name} with error:\n{e}')
        logging.warning(
            'MSA plots are only used to evaluate TEtrimmer and will not affect the final TE consensus library.'
            " For traceback text, please refer to 'error_file.txt' in the 'Multiple_sequence_alignment' folder\n"
        )

    #####################################################################################################
    # Code block: Generate TE-Aid plot for input and output TE consensus sequence, and out boundary MSA
    #####################################################################################################
    try:
        # so.path.abspath(__file__) will return the current executable Python file
        # TE-Aid package is stored in the same directory as this function file.
        TE_aid_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), 'TE-Aid-master'
        )

        # Because terminal repeats were found before, use the previous result
        # The low_copy will only affect it to keep self blast result from TEAid.
        # it is Ture, TE_aid_object.run will check the terminal repeat based on the self blast output of TEAid
        # Otherwise, it uses the terminal repeat result "found_match_crop"

        TE_aid_object = TEAid(
            orf_cons,
            output_dir,
            genome_file,
            blast_database_path,
            mmseqs_database_dir,
            error_file=error_files,
            TE_aid_dir=TE_aid_path,
        )
        TE_aid_plot = TE_aid_object.run(title="after")

        if found_match_crop in ("LTR", "TIR"):
            found_match = found_match_crop
        else:
            # found_match could be False LTR or TIR
            found_match = TE_aid_object.teaid_check_termina_repeat()

        # Run TE_aid to plot the query sequence, if required. Because one query file can return multiple
        # clusters, TE-Aid will check if a TE-Aid plot has been generated before.
        # Set low_
        if plot_query:
            query_file = seq_obj.get_input_fasta()
            TE_aid_object_query = TEAid(
                query_file,
                output_dir,
                genome_file,
                blast_database_path,
                mmseqs_database_dir,
                error_file=error_files,
                TE_aid_dir=TE_aid_path,
            )
            TE_aid_plot_query = TE_aid_object_query.run(title="before")

            found_match_query = TE_aid_object_query.teaid_check_termina_repeat()

            all_blast_hit_n_query, full_blast_query_n = TE_aid_object_query.check_blast_full_n(seq_obj, check_query=True)

            input_te_genome_cov_len = TE_aid_object_query.te_genome_coverage()

            # Update input sequence terminal repeat information
            seq_obj.set_old_terminal_repeat(found_match_query)

            # Update input sequence full length blast number
            seq_obj.set_old_blast_full_n(full_blast_query_n)

            # Update input sequence genome coverage length number
            seq_obj.set_input_genome_cov_len(input_te_genome_cov_len)

        else:
            TE_aid_plot_query = None

        # Do TEAid for out boundary output TE consensus sequence
        out_boundary_msa_for_teaid_cons = con_generater(
            out_boundary_msa_for_teaid,
            output_dir,
            threshold=0.3, ambiguous='N'
        )

        TE_aid_object_out_boundary_msa = TEAid(
            out_boundary_msa_for_teaid_cons,
            output_dir,
            genome_file,
            blast_database_path,
            mmseqs_database_dir,
            error_file=error_files,
            TE_aid_dir=TE_aid_path,
        )

        TE_aid_object_out_boundary_msa_plot = TE_aid_object_out_boundary_msa.run(
            v_x_line_1=start_relate_to_out_boundary_msa_for_teaid,
            v_x_line_2=end_relate_to_out_boundary_msa_for_teaid,
            title="extend"
        )

    except Exception as e:
        with open(error_files, 'a') as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f'\nTE-Aid plot failed for {seq_name} with error:\n{e}')
            f.write('\n' + tb_content + '\n\n')
        logging.error(f'\nTE-Aid plot failed for {seq_name} with error:\n{e}\n' + tb_content + '\n')
        raise Exception from e

    #####################################################################################################
    # Code block: Generate dotplot
    #####################################################################################################
    dotplot_pdf = None
    try:
        # Use reverse complemented input sequence file, when reverse_complement is true
        if reverse_complement:
            dotplot_pdf = dotplot(orf_cons, seq_file_reverse_c, output_dir)
        else:
            dotplot_pdf = dotplot(orf_cons, seq_file, output_dir)
    except Exception as e:
        # dotplot is not mandatory; skip if any error occurred
        with open(error_files, 'a') as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f'\nDotplot failed for {seq_name} with error:\n{e}')
            f.write('\n' + tb_content + '\n\n')

    try:
        if dotplot_pdf is not None:
            # Because dotmatcher can't change output size, scale it up to make it more clear in the merged pdf
            scale_dotplot_pdf = scale_single_page_pdf(
                dotplot_pdf, f'{dotplot_pdf}_su.pdf', scale_ratio=2
            )
            dotplot_pdf = scale_dotplot_pdf
    except Exception as e:
        with open(error_files, 'a') as f:
            # This is not mandatory, skip if any error occurred
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(
                f'\nps file to pdf conversion failed for {seq_name} with error:\n{e}'
            )
            f.write('\n' + tb_content + '\n\n')

    #####################################################################################################
    # Code block: Merge plot files
    #####################################################################################################
    try:
        merged_pdf_path = merge_pdfs(
            output_dir,
            os.path.basename(cropped_boundary_MSA),
            MSA_plot,
            cropped_boundary_manual_MSA_concatenate_plot,
            TE_aid_plot,
            TE_aid_object_out_boundary_msa_plot,
            TE_aid_plot_query,
            orf_domain_plot,
            input_orf_pfam,
            dotplot_pdf,
        )

    except Exception as e:
        with open(error_files, 'a') as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f'\nPlot merging to PDF failed for {seq_name} with error:\n{e}')
            f.write('\n' + tb_content + '\n\n')
        logging.error(f'\nPlot merging to PDF failed for {seq_name} with error:\n{e}\n{tb_content}\n')
        raise Exception from e

    #####################################################################################################
    # Code block: Create directories
    #####################################################################################################
    try:
        # Create a folder in the same directory as output_dir to store annotation files
        parent_output_dir = os.path.dirname(output_dir)

        # Define different levels of proof_curation folder
        perfect_proof = os.path.join(proof_curation_dir, 'Annotations_perfect')
        good_proof = os.path.join(proof_curation_dir, 'Annotations_good')
        intermediate_proof = os.path.join(
            proof_curation_dir, 'Annotations_check_recommended'
        )
        need_check_proof = os.path.join(
            proof_curation_dir, 'Annotations_check_required'
        )

        # Create the directory if it does not exist
        os.makedirs(proof_curation_dir, exist_ok=True)
        os.makedirs(perfect_proof, exist_ok=True)
        os.makedirs(good_proof, exist_ok=True)
        os.makedirs(intermediate_proof, exist_ok=True)
        os.makedirs(need_check_proof, exist_ok=True)

        # Define temporary classified and unknown final consensus file, which will be used for
        # reclassification by RepeatMasker
        final_unknown_con_file = os.path.join(
            classification_dir, 'temp_TEtrimmer_unknown_consensus.fasta'
        )
        final_classified_con_file = os.path.join(
            classification_dir, 'temp_TEtrimmer_classified_consensus.fasta'
        )

        # Generate final consensus sequence
        # Compared with initial consensus, the final consensus sequence can have a different orientation.
        # For this reason, generate consensus again
        final_con = con_generater_no_file(
            cropped_boundary_MSA, threshold=cons_threshold
        )
        sequence = str(final_con).upper()

        #####################################################################################################
        # Code block: Calculate input and output sequence identity and coverage
        #####################################################################################################

        try:
            if reverse_complement:
                input_output_identity, coverage_input_seq, coverage_output_seq = pairwise_seqs_align(
                    seq_file_reverse_c,
                    sequence,
                    seq2_is_file=False
                )
            else:
                input_output_identity, coverage_input_seq, coverage_output_seq = pairwise_seqs_align(
                    seq_file,
                    sequence,
                    seq2_is_file=False
                )

            consi_obj.set_in_out_identity(input_output_identity)
            consi_obj.set_input_coverage(coverage_input_seq)
            consi_obj.set_output_coverage(coverage_output_seq)

        except Exception as e:
            # dotplot is not mandatory; skip if any error occurred
            with open(error_files, 'a') as f:
                # Get the traceback content as a string
                tb_content = traceback.format_exc()
                f.write(f'\nInput and output sequence pairwise alignment failed for {seq_name} with error:\n{e}')
                f.write('\n' + tb_content + '\n\n')
            logging.error(f'\nInput and output sequence pairwise alignment failed for {seq_name} with error:\n{e}'
                          f'\n' + tb_content + '\n\n')

        # Storing consensus sequence length into consi_obj
        consi_obj.set_new_length(len(sequence))

        # Store consensus sequence into consi_obj not necessary
        # consi_obj.set_cons_seq(sequence)

        # Store MSA sequence number into consi_obj
        MSA_for_final_cons = AlignIO.read(cropped_boundary_MSA, 'fasta')
        MSA_for_final_cons_seq_n = len(MSA_for_final_cons)
        consi_obj.set_cons_MSA_n(MSA_for_final_cons_seq_n)

        # Store terminal repeat to consi_obj
        if found_match == 'LTR' or found_match == 'TIR':
            consi_obj.set_new_terminal_repeat(found_match)

        # Store full length sequence number from BLAST search into consi_obj
        # check_blast_full_n is a function in TE_Aid class. TEtrimmer will use the blast result of TE Aid
        all_blast_hit_n_con, blast_full_length_n = TE_aid_object.check_blast_full_n(consi_obj, check_query= False, engine=engine)
        output_genome_cov_len = TE_aid_object.te_genome_coverage()

        # Set output sequence object informations
        consi_obj.set_cons_blast_n(all_blast_hit_n_con)
        consi_obj.set_blast_full_n(blast_full_length_n)
        consi_obj.set_output_genome_cov_len(output_genome_cov_len)

        # Store PFAM predictions to consi_obj
        if if_pfam_domain:
            consi_obj.set_cons_pfam(if_pfam_domain)

    except Exception as e:
        with open(error_files, 'a') as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f'\nUpdate sequence object failed for {seq_name} with error:\n{e}')
            f.write('\n' + tb_content + '\n\n')
        logging.error(f'\nUpdate sequence object failed for {seq_name} with error:\n{e}\n{tb_content}\n')
        raise Exception from e

    #####################################################################################################
    # Code block: Run RepeatClassifier in RepeatModeler to classify TEtrimmer consensus sequences
    #####################################################################################################

    # This classification is different from the final RepeatMasker classification
    # check_unknown returns 'True' if unknowns were detected
    # Classify all elements by RepeatClassify if classify_all is 'True'
    # Rename consensus if classify_unknown is 'True' and if the final consensus length is much longer or
    # shorter than the query sequence or if the TE type is unknown.
    # fast_mode will supress RepeatClassifier step
    # Classification is not mandatory, skip this step if errors are encountered
    try:
        if classify_all or (classify_unknown and (seq_obj.check_unknown() or (left_ex_total + right_ex_total >= 5000))):

            # Define different folders for each sequence
            # the suffix .fasta is important, this ensures that this folder can be deleted later
            classification_seq_folder = os.path.join(
                classification_dir, f'{uniq_seq_name}.fasta'
            )
            os.makedirs(classification_seq_folder, exist_ok=True)

            # Define consensus file path used for classification
            classification_seq_file = os.path.join(
                classification_seq_folder, uniq_seq_name
            )

            with open(classification_seq_file, 'w') as f:
                # RepeatClassifier input cannot be a single sequence, add >Dummy to enable RepeatClassifier run
                f.write(
                    '>'
                    + uniq_seq_name
                    + '\n'
                    + sequence
                    + '\n'
                    + '>Dummy'
                    + '\n'
                    + 'T'
                    + '\n'
                )

            TE_type = classify_single(classification_seq_file)

            # Only update new_TE_type if classify_single was successful
            if TE_type:
                # Set TE_type after RepeatClassifier
                consi_obj.set_new_TE_type(TE_type)

            # Clean classification folder if debug is 'True'
            if not debug:
                remove_files_with_start_pattern(
                    classification_dir, f'{uniq_seq_name}.fasta'
                )
    except Exception as e:
        with open(error_files, 'a') as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f'\nRepeatClassifier classification error for {seq_name}')
            f.write('\n' + tb_content + '\n\n')
        logging.error(f"\nNote: RepeatClassifier doesn't work for {seq_name} with error {e}")
        logging.warning(
            "\nThis won't affect final TE consensus sequences but only the classification. You can choose to ignore this. "
            "For traceback text, please refer to 'error_file.txt' under 'Multiple_sequence_alignment' folder\n"
        )

    # Update final con_TE_type. get_TE_type_for_file will evaluate if TE_type is Unknown. If so, use the
    # original TE classification name
    updated_TE_type = consi_obj.get_TE_type_for_file()
    consi_obj.set_new_TE_type(updated_TE_type)

    ###############################################################################################################
    # Code block: Output evaluation (perfect, good, Reco_check, need_check)
    ###############################################################################################################

    #               Terminal_repeat    Classified    MSA_sequence_number    Blast_full_length_number    if_PFAM
    # Perfect       True               True          >=30                   >=5                         True
    # Good          True               Not_required  >=10                   >=2                         Not_required
    # Reco_check    Not_required       Not_required  >=20                   >=2                         Not_required
    # Need_check    Not_required       Not_required  Not_required           Not_required                Not_required


    # Perfect
    # Good
    # Reco_check
    # Need_check


    try:
        if not extension_enough:
            consi_obj.set_cons_evaluation('Need_ext')
        elif (
            consi_obj.new_TE_terminal_repeat != 'False'
            and consi_obj.new_TE_type != 'NaN'
            and 'unknown' not in consi_obj.new_TE_type.lower()
            and consi_obj.new_TE_MSA_seq_n >= perfect_seq_num
            and consi_obj.new_TE_blast_full_length_n >= 5
            and consi_obj.cons_pfam
        ):
            consi_obj.set_cons_evaluation('Perfect')

        elif (
            consi_obj.new_TE_terminal_repeat != 'False'
            and consi_obj.new_TE_MSA_seq_n >= 10
            and consi_obj.new_TE_blast_full_length_n >= 2
        ):
            consi_obj.set_cons_evaluation('Good')

        elif (
            consi_obj.new_TE_MSA_seq_n >= 20
            and consi_obj.new_TE_blast_full_length_n >= 2
        ):
            consi_obj.set_cons_evaluation('Reco_check')

        else:
            consi_obj.set_cons_evaluation('Need_check')

        #####################################################################################################
        # Code block: Code block: Move file for manual inspection
        #####################################################################################################

        consi_obj.set_proof_curation_file()

        # modify fasta_out_flank_mafft_gap_rm fasta header based on the bed file, this can allow the
        # extension function in the final GUI.
        # For example: change 1(+) to scaffold_1:23256-24757(+)
        # The new header will allow the "Extend" function button in the TEtrimmerGUI
        cropped_boundary_MSA_nm = modify_fasta_headers(
            bed_out_flank_file, cropped_boundary_MSA
        )
        bed_fasta_mafft_boundary_crop_for_select_nm = modify_fasta_headers(
            bed_out_flank_file, bed_fasta_mafft_boundary_crop_for_select
        )

        bed_fasta_mafft_with_gap_nm = modify_fasta_headers(
            bed_out_flank_file, bed_fasta_mafft_with_gap
        )

        # Define file name for inspection file
        file_copy_pattern = [
            (merged_pdf_path, str(consi_obj.proof_pdf)),
            (cropped_boundary_MSA_nm, str(consi_obj.proof_fasta)),
            (bed_fasta_mafft_with_gap_nm, str(consi_obj.proof_raw)),
            (cluster_msa, str(consi_obj.proof_cluster))
        ]

        files_moved_successfully = True

        for pattern, new_name in file_copy_pattern:
            try:
                if consi_obj.cons_evaluation == 'Need_ext':
                    destination_dir = more_extension_dir
                elif consi_obj.cons_evaluation == 'Perfect':
                    destination_dir = perfect_proof
                elif consi_obj.cons_evaluation == 'Good':
                    destination_dir = good_proof
                elif consi_obj.cons_evaluation == 'Reco_check':
                    destination_dir = intermediate_proof
                else:
                    destination_dir = need_check_proof

                # Copy the file to the new location with the new unique name
                if pattern:
                    shutil.copy(pattern, os.path.join(destination_dir, new_name))

            except Exception as e:
                with open(error_files, 'a') as f:
                    # Get the traceback content as a string
                    tb_content = traceback.format_exc()
                    f.write(f'Copy file error for {pattern}\n')
                    f.write(tb_content + '\n\n')
                logging.error(f'Error copying {pattern} to {new_name}: {e}')
                files_moved_successfully = False

        if hmm:  # Generate HMM files
            consi_obj.set_hmm_file()
            hmm_output_file = os.path.join(hmm_dir, consi_obj.hmm_file)
            generate_hmm_from_msa(cropped_boundary_MSA, hmm_output_file, error_files)

        # Classification of unknown consensus TEs will be attempted again later by using successfully classified sequences
        if 'Unknown' in updated_TE_type:
            with open(final_unknown_con_file, 'a') as f:  # 'a' mode for appending
                f.write('>' + uniq_seq_name + '\n' + sequence + '\n')
        else:
            with open(final_classified_con_file, 'a') as f:
                f.write(
                    '>' + uniq_seq_name + '#' + updated_TE_type + '\n' + sequence + '\n'
                )

        # Don't write file to the final_con_file if extension is not enough

        if extension_enough:
            # Write all consensus sequence to final_cons_file.
            with open(final_con_file, 'a') as f:
                f.write(
                    '>' + uniq_seq_name + '#' + updated_TE_type + '\n' + sequence + '\n'
                )

            # Write all consensus sequences to final_cons_file_no_low_copy.
            with open(final_con_file_no_low_copy, 'a') as f:
                f.write(
                    '>' + uniq_seq_name + '#' + updated_TE_type + '\n' + sequence + '\n'
                )

        return files_moved_successfully

    except Exception as e:
        with open(error_files, 'a') as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f'\nMoving of files failed for {seq_name} with error:\n{e}')
            f.write('\n' + tb_content + '\n\n')
        logging.error(f'\nMoving of files failed for {seq_name} with error:\n{e}\n{tb_content}\n')
        raise Exception from e
