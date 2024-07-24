# Standard library imports
import os
import click
import shutil
import traceback
import pandas as pd
from Bio import AlignIO



# Local imports
from functions import generate_hmm_from_msa, extract_fasta, remove_gaps_with_similarity_check, align_sequences, \
    con_generater_no_file, concatenate_alignments, select_window_columns, select_start_end_and_join, \
    con_generater, reverse_complement_seq_file, classify_single, check_terminal_repeat, select_star_to_end, \
    define_crop_end_simi_thr, prcyan, prgre, merge_pdfs, dotplot, scale_single_page_pdf, \
    remove_files_with_start_pattern, find_poly_a_end_position, is_LTR, check_and_update, modify_fasta_headers
from boundaryclass import CropEnd, CropEndByGap, DefineBoundary
from TEaid import TEAid
from orfdomain import PlotPfam, determine_sequence_direction
from MSAcluster import CleanAndSelectColumn, process_msa, clean_and_cluster_MSA
import cialign


def long_bed(input_file, output_dir):
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
    reset_right_filtered_df.loc[reset_right_filtered_df[5] == '+', 1] = reset_right_filtered_df.loc[reset_right_filtered_df[5] == '+', 2] - 1
    reset_right_filtered_df.loc[reset_right_filtered_df[5] == '-', 2] = reset_right_filtered_df.loc[reset_right_filtered_df[5] == '-', 1] + 1

    # Conditionally update values in column 1 and column 2 for left extension
    reset_left_filtered_df = filtered_df.copy()
    reset_left_filtered_df.loc[reset_left_filtered_df[5] == '+', 2] = reset_left_filtered_df.loc[reset_left_filtered_df[5] == '+', 1] + 1
    reset_left_filtered_df.loc[reset_left_filtered_df[5] == '-', 1] = reset_left_filtered_df.loc[reset_left_filtered_df[5] == '-', 2] - 1

    # write new BED file
    reset_right_long_bed = os.path.join(output_dir, f"{os.path.basename(input_file)}_reset_right_long.bed")
    reset_left_long_bed = os.path.join(output_dir, f"{os.path.basename(input_file)}_reset_left_long.bed")
    reset_right_filtered_df.to_csv(reset_right_long_bed, sep='\t', index=False, header=None)
    reset_left_filtered_df.to_csv(reset_left_long_bed, sep='\t', index=False, header=None)

    return reset_left_long_bed, reset_right_long_bed


def extend_end(max_extension, ex_step_size, end, input_file, genome_file, output_dir, crop_end_gap_thr,
               crop_end_gap_win, ext_threshold, define_boundary_win, bed_dic):
    # end: left or right extension
    if_ex = True
    ex_total = 0

    # The following code calculates the majority of the alignment length and finds the longest alignment for extension
    # while ignoring the shorter ones.
    # long_bed function won't do length selection
    # It generates the proper bed file for left and right extension.
    reset_left_long_bed, reset_right_long_bed = long_bed(input_file, output_dir)

    # 100 bp are added to avoid a case where the boundary is at the edge. Therefore, the alignment length is actually
    # ex_step_size + 100 bp
    adjust = 100
    ite = 1
    reset_left_long_bed_copy = reset_left_long_bed
    reset_right_long_bed_copy = reset_right_long_bed
    while if_ex and (ex_total < max_extension):
        # BEDtools will make sure that the extension will not excess the maximum length of the chromosome
        # track extend total
        ex_total += ex_step_size
        
        if end == "left":
            left_ex = ex_total
            right_ex = ex_step_size - ex_total + adjust
            bed_fasta, bed_out_flank_file = extract_fasta(reset_left_long_bed, genome_file, output_dir,
                                                          left_ex, right_ex, nameonly=True)
        if end == "right":
            left_ex = ex_step_size - ex_total + adjust
            right_ex = ex_total
            bed_fasta, bed_out_flank_file = extract_fasta(reset_right_long_bed, genome_file, output_dir,
                                                          left_ex, right_ex, nameonly=True)
        
        # align_sequences() will return extended MSA in absolute file
        bed_fasta_mafft_with_gap = align_sequences(bed_fasta, output_dir)

        if not os.path.isfile(bed_fasta_mafft_with_gap):
            click.echo(f"{input_file} encountered a problem during the MAFFT extension step.")
            return False

        # Remove nucleotide whose proportion is smaller than threshold
        bed_fasta_mafft_with_gap_column_clean_object = CleanAndSelectColumn(bed_fasta_mafft_with_gap, threshold=0.08)
        bed_fasta_mafft_with_gap_column_clean = bed_fasta_mafft_with_gap_column_clean_object.clean_column(output_dir)

        # Remove gaps with similarity check
        bed_fasta_mafft = remove_gaps_with_similarity_check(bed_fasta_mafft_with_gap_column_clean, output_dir,
                                                            gap_threshold=0.6, simi_check_gap_thre=0.4,
                                                            similarity_threshold=0.7,
                                                            min_nucleotide=5
                                                           )

        # bed_fasta_mafft will be false if the MSA column number is smaller than 50 after removing the gaps
        if not bed_fasta_mafft:
            break
        bed_fasta_mafft_object = CropEndByGap(bed_fasta_mafft, gap_threshold=crop_end_gap_thr,
                                              window_size=crop_end_gap_win)
        large_crop_ids, remain_n, remain_ids = bed_fasta_mafft_object.find_large_crops()
        cropped_alignment = bed_fasta_mafft_object.crop_alignment()
        bed_fasta_mafft_cop_end_gap = bed_fasta_mafft_object.write_to_file(output_dir, cropped_alignment)

        # ext_threshold refer to option --ext_thr
        # max_X is float number, means the proportion of X
        bed_boundary = DefineBoundary(bed_fasta_mafft_cop_end_gap, threshold=ext_threshold,
                                      check_window=define_boundary_win, max_X=0.3, if_con_generater=False,
                                      extension_stop_num=300)
        # Read bed_out_flank_file
        bed_out_flank_file_df = pd.read_csv(bed_out_flank_file, sep='\t', header=None)

        # Update bed_dic
        for i in bed_dic:
            id_list = large_crop_ids + remain_ids
            if i in id_list:
                matching_rows = bed_out_flank_file_df[bed_out_flank_file_df[3] == i]
                if end == "left":
                    if bed_dic[i][2] == "+":
                        bed_dic[i][0] = matching_rows.iloc[0, 1]
                    elif bed_dic[i][2] == "-":
                        bed_dic[i][1] = matching_rows.iloc[0, 2]
                if end == "right":
                    if bed_dic[i][2] == "+":
                        bed_dic[i][1] = matching_rows.iloc[0, 2]
                    elif bed_dic[i][2] == "-":
                        bed_dic[i][0] = matching_rows.iloc[0, 1]

        # If the remaining sequences that need further extension are fewer than 5, stop the extension
        if not bed_boundary.if_continue or remain_n < 5:
            break
        elif bed_boundary.if_continue:
            # Subset BED file to only keep sequences that need further extension
            if end == "left":
                if_ex = bed_boundary.left_ext
                df = pd.read_csv(reset_left_long_bed, sep='\t', header=None)

                # Filter out rows where the fourth column is in large_crop_ids
                df_filtered = df[~df[3].isin(large_crop_ids)]
                reset_left_long_bed = reset_left_long_bed_copy + f"_{ite}"

                # Write the filtered DataFrame back to reset_left_long_bed
                df_filtered.to_csv(reset_left_long_bed, sep='\t', index=False, header=None)

            if end == "right":
                if_ex = bed_boundary.right_ext
                df = pd.read_csv(reset_right_long_bed, sep='\t', header=None)

                # Filter out rows where the fourth column is in large_crop_ids
                df_filtered = df[~df[3].isin(large_crop_ids)]
                reset_right_long_bed = reset_right_long_bed_copy + f"_{ite}"

                # Write the filtered DataFrame into reset_left_long_bed
                df_filtered.to_csv(reset_right_long_bed, sep='\t', index=False, header=None)
            ite += 1

    return bed_dic, ex_total


def final_MSA(bed_file, genome_file, output_dir, gap_nul_thr, gap_threshold, ext_threshold,
              define_boundary_win, crop_end_gap_thr, crop_end_gap_win, crop_end_thr, crop_end_win):
    found_match_crop = False

    bed_fasta, bed_out_flank_file = extract_fasta(bed_file, genome_file, output_dir, 0, 0, nameonly=True)

    # align_sequences() will return extended MSA in absolute file
    bed_fasta_mafft_with_gap = align_sequences(bed_fasta, output_dir)

    if not os.path.isfile(bed_fasta_mafft_with_gap):
        click.echo(f"{bed_file} encountered a problem during the final MAFFT extension step.")
        return False

    # Remove nucleotides whose proportion is smaller than threshold
    bed_fasta_mafft_with_gap_column_clean_object = CleanAndSelectColumn(bed_fasta_mafft_with_gap, threshold=0.01)
    bed_fasta_mafft_with_gap_column_clean = bed_fasta_mafft_with_gap_column_clean_object.clean_column(output_dir)

    # Remove gaps by similarity check
    bed_fasta_mafft_gap_sim = remove_gaps_with_similarity_check(bed_fasta_mafft_with_gap_column_clean, output_dir)

    # Crop end before terminal check, use loose threshold to avoid overhead cropping
    bed_fasta_mafft_gap_sim_cp_object = CropEnd(bed_fasta_mafft_gap_sim, threshold=20, window_size=40)
    bed_fasta_mafft_gap_sim_cp_object.crop_alignment()
    bed_fasta_mafft_gap_sim_cp = bed_fasta_mafft_gap_sim_cp_object.write_to_file(output_dir)

    # Generate consensus sequence, use low threshold to reduce number of N's in the consensus sequence
    bed_fasta_mafft_gap_sim_cp_con = con_generater(bed_fasta_mafft_gap_sim_cp, output_dir, threshold=0.45, ambiguous="N")

    # Check terminal repeats
    # Define the output folder to store the temporary BLAST database
    check_terminal_repeat_output = f"{bed_fasta_mafft_gap_sim_cp_con}_tem"
    LTR_boundary, TIR_boundary, found_match_crop = check_terminal_repeat(bed_fasta_mafft_gap_sim_cp_con,
                                                                         check_terminal_repeat_output)

    # Check terminal repeats
    if found_match_crop:
        # Check LTRs first
        if LTR_boundary is not None:
            left = LTR_boundary[0]
            right = LTR_boundary[1]
        else:
            left = TIR_boundary[0]
            right = TIR_boundary[1]

        # Extract MSA based on terminal repeat boundary
        bed_fasta_mafft_gap_sim_selected = select_star_to_end(bed_fasta_mafft_gap_sim_cp, output_dir, left, right)

        # Define the threshold use for cropping end by similarity to make sure the ends of the MSA are not totally removed
        # because the boundary of the MSA has been defined.
        start_pro_mean, end_pro_mean = define_crop_end_simi_thr(bed_fasta_mafft_gap_sim_selected)

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
        bed_fasta_mafft_gap_sim_selected_cp_object_left = CropEnd(
            bed_fasta_mafft_gap_sim_selected, threshold=start_pro_mean, window_size=40, crop_l=True, crop_r=False)
        bed_fasta_mafft_gap_sim_selected_cp_object_left.crop_alignment()
        bed_fasta_mafft_gap_sim_selected_cp_left = bed_fasta_mafft_gap_sim_selected_cp_object_left.write_to_file(output_dir)

        # Crop end by similarity to clean MSA for right side
        bed_fasta_mafft_gap_sim_selected_cp_object_right = CropEnd(
            bed_fasta_mafft_gap_sim_selected_cp_left, threshold=end_pro_mean, window_size=40, crop_l=False, crop_r=True)
        bed_fasta_mafft_gap_sim_selected_cp_object_right.crop_alignment()
        bed_fasta_mafft_gap_sim_selected_cp = bed_fasta_mafft_gap_sim_selected_cp_object_right.write_to_file(output_dir)

        # Remove gaps by similarity again
        bed_fasta_mafft_gap_sim_selected_cp_g, column_mapping = remove_gaps_with_similarity_check(
            bed_fasta_mafft_gap_sim_selected_cp, output_dir, return_map=True)

        # Modify the dictionary to correspond with bed_fasta_mafft_gap_sim
        for key in column_mapping:
            column_mapping[key] += left

        """
        # Crop end by gap to make final MSA look better
        bed_fasta_mafft_gap_sim_selected_cp_g_cpg = CropEndByGap(bed_fasta_mafft_gap_sim_selected_cp_g,
                                                                 gap_threshold=crop_end_gap_thr,
                                                                 window_size=crop_end_gap_win)
        """
        # Set max_X to 1, this will make sure this part won't cut any beginning and end columns. Because
        # the boundary has been defined by terminal repeat. The reason to do DefineBoundary is cropped_boundary
        # and cropped_boundary_MSA are required for the further analysis
        # max_X is the proportion of X, the maximum number is 1
        cropped_boundary = DefineBoundary(bed_fasta_mafft_gap_sim_selected_cp_g, threshold=0.8,
                                          check_window=4, max_X=1)
        cropped_boundary_MSA = cropped_boundary.crop_MSA(output_dir, crop_extension=0)

        cropped_alignment_output_file_g = bed_fasta_mafft_gap_sim_selected_cp_g
        bed_fasta_mafft_boundary_crop_for_select = bed_fasta_mafft_gap_sim

    # If LTR or TIR are not found
    else:
        # Check if poly A can be found from the sequence and return the last A position
        # poly_a will be None if not found
        poly_a = find_poly_a_end_position(bed_fasta_mafft_gap_sim_cp_con, min_length=10)

        bed_boundary = DefineBoundary(bed_fasta_mafft_gap_sim_cp, threshold=ext_threshold,
                                      check_window=define_boundary_win, max_X=0.2,
                                      if_con_generater=False, end_position=poly_a)
        # when if_continue is false, it means the start position is greater than end position
        if bed_boundary.if_continue:
            bed_fasta_mafft_boundary_crop = bed_boundary.crop_MSA(output_dir, crop_extension=300)

            # For LINE elements bed_fasta_mafft_boundary_crop will be changed. Copy alignment to another variable
            bed_fasta_mafft_boundary_crop_for_select = bed_fasta_mafft_boundary_crop

            if "LINE" in bed_file:
                # For highly divergent regions, more gaps can be found. According to this feature, remove
                # high-divergence regions. This function is very useful for dealing with LINE elements.
                cropped_MSA_by_gap = CropEndByGap(bed_fasta_mafft_boundary_crop, gap_threshold=crop_end_gap_thr,
                                                  window_size=crop_end_gap_win)
                bed_fasta_mafft_boundary_crop = cropped_MSA_by_gap.write_to_file(output_dir)

            # Gaps are removed again after processing with CropEnd
            cropped_alignment_output_file_g, column_mapping = crop_end_and_remove_gap(
                bed_fasta_mafft_boundary_crop, output_dir, crop_end_threshold=crop_end_thr,
                window_size=crop_end_win, gap_threshold=0.8)

            # This means the sequence length after gap removal is shorter than 50
            if not cropped_alignment_output_file_g:
                return False

            # CropEnd cannot define the final boundary. Use DefineBoundary again to define start position.
            cropped_boundary = DefineBoundary(cropped_alignment_output_file_g, threshold=0.8,
                                              check_window=4, max_X=0)
            cropped_boundary_MSA = cropped_boundary.crop_MSA(output_dir, crop_extension=0)
        else:
            return False

    return bed_out_flank_file, cropped_boundary_MSA, cropped_alignment_output_file_g, cropped_boundary, \
        column_mapping, bed_fasta_mafft_boundary_crop_for_select, found_match_crop


def crop_end_and_remove_gap(input_file, output_dir, crop_end_threshold=0.8, window_size=20, gap_threshold=0.8):

    # window_size means the checked nucleotide number each time
    # threshold means the nucleotide proportion in window_size must greater than INT
    CropEnd_thres = crop_end_threshold * window_size
    cropped_MSA = CropEnd(input_file, threshold=CropEnd_thres, window_size=window_size)
    cropped_MSA.crop_alignment()
    
    # write_to_file() function will return absolute cropped alignment file
    cropped_MSA_output_file = cropped_MSA.write_to_file(output_dir)

    # Remove gaps again after the CropEnd step
    # "column_mapping" is a dictionary; the key is gap removed MSA index, the value is the corresponding original
    # MSA index
    cropped_alignment_output_file_g, column_mapping = remove_gaps_with_similarity_check(
        cropped_MSA_output_file, output_dir, return_map=True)

    return cropped_alignment_output_file_g, column_mapping


def find_boundary_and_crop(bed_file, genome_file, output_dir, pfam_dir, seq_obj, hmm, classify_all, classify_unknown,
                           error_files, plot_query, classification_dir, final_con_file, final_con_file_no_low_copy,
                           proof_curation_dir, hmm_dir,
                           cons_threshold=0.8, ext_threshold=0.7, ex_step_size=1000,
                           max_extension=7000, gap_threshold=0.4, gap_nul_thr=0.7, crop_end_thr=0.8, crop_end_win=40,
                           crop_end_gap_thr=0.1, crop_end_gap_win=150, start_patterns=None, end_patterns=None,
                           mini_orf=200, define_boundary_win=150, fast_mode=False, engine="blast",
                           input_orf_pfam=False, debug=False, cluster_msa=None):
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
    :param start_patterns: str, patterns to check for start points. Default: None
    :param end_patterns: str, patterns to check for end points. Default: None
    :param min_orf: int, set minimum ORF length for ORF prediction. Default: 200

    """

    # Check if this is a LINE element, if so decrease the ext_threshold number,
    # because LINE have high divergence at the 5' end
    seq_name = seq_obj.name
    seq_file = seq_obj.get_input_fasta()  # Return full file path
    if "LINE" in seq_obj.old_TE_type:
        ext_threshold = ext_threshold - 0.1

    msa_loop_n = 1

    while msa_loop_n <= 2:

        # TEtrimmer will do the second round clustering and extension. The MSA has been extensively extended during
        # the first time of loop, reduce the max_extension to 2000 for the second round, which can prevent to generate
        # very long sequence especially for tandem repeats.
        if msa_loop_n == 2:
            max_extension = 2*ex_step_size

        try:
            # Read BED file and build dictionary to store sequence position
            df = pd.read_csv(bed_file, sep='\t', header=None)

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
            with open(error_files, "a") as f:
                # Get the traceback content as a string
                tb_content = traceback.format_exc()
                f.write(f"\nFind boundary bed_file reading failed for {seq_name} with error:\n{e}")
                f.write('\n' + tb_content + '\n\n')
                prcyan(f"\nStart and end patterns like 'TGT' or 'ACA' definition failed for {seq_name} with error \n{e}")
                prgre("\nThis error will not affect the final result significantly, you can choose to ignore it. "
                      "For traceback text, please refer to 'error_file.txt' in the 'Multiple_sequence_alignment' folder.\n")
                pass

        #####################################################################################################
        # Code block: Define length of extensions on left and right sides of the MSA and define the boundary
        #####################################################################################################
        try:
            left_bed_dic, left_ex_total = extend_end(
                max_extension, ex_step_size, "left", bed_file, genome_file, output_dir, crop_end_gap_thr,
                crop_end_gap_win, ext_threshold, define_boundary_win, bed_dict)

            final_bed_dic, right_ex_total = extend_end(
                max_extension, ex_step_size, "right", bed_file, genome_file, output_dir, crop_end_gap_thr,
                crop_end_gap_win, ext_threshold, define_boundary_win, left_bed_dic)
        except Exception as e:
            with open(error_files, "a") as f:
                # Get the traceback content as a string
                tb_content = traceback.format_exc()
                f.write(f"\nMSA extension failed for {seq_name} with error:\n{e}")
                f.write('\n' + tb_content + '\n\n')
                prcyan(f"\nMSA extension failed for {seq_name} with error:\n{e}")
                prcyan('\n' + tb_content + '\n')
                raise Exception

        try:
            # Update BED file for final MSA
            for key, value in final_bed_dic.items():
                # Find rows where the fourth column matches the key and update position by dictionary value
                df.loc[df[3] == key, [1, 2]] = value[:2]

            # Define BED file name for final MSA
            bed_final_MSA = bed_file + "_fm.bed"
            df.to_csv(bed_final_MSA, sep='\t', index=False, header=False)

            # final_MSA return 'False' if when the start crop point is greater than the end crop point, which means the
            # sequence is too short. If the divergence is too high, CropEnd will return 'False'.
            final_msa_result = final_MSA(bed_final_MSA, genome_file, output_dir, gap_nul_thr, gap_threshold, ext_threshold,
                                         define_boundary_win, crop_end_gap_thr, crop_end_gap_win, crop_end_thr, crop_end_win)

            # Check if final_msa_result returned 'False', indicating an error or specific condition
            if not final_msa_result:
                # Handle the error or specific condition
                return False
            else:
                # Unpack the returned values if the function executed normally
                bed_out_flank_file, cropped_boundary_MSA, cropped_alignment_output_file_g, cropped_boundary, \
                    column_mapping, bed_fasta_mafft_boundary_crop_for_select, found_match_crop = final_msa_result

        except Exception as e:
            with open(error_files, "a") as f:
                # Return the traceback content as a string
                tb_content = traceback.format_exc()
                f.write(f"Boundary definition error {seq_name}\n")
                f.write('\n' + tb_content + '\n\n')
                prcyan(f"\nBoundary definition failed for {seq_name} with error:\n{e}")
                prcyan('\n' + tb_content + '\n')
                raise Exception

        #####################################################################################################
        # Code block: Check the consistency of the final MSA
        #####################################################################################################
        # Plotting is done after PFAM predictions in case consensus/MSA are in the wrong direction.
        try:
            if msa_loop_n <= 1:
                final_msa_consistency = clean_and_cluster_MSA(cropped_boundary_MSA, bed_out_flank_file,
                                                              output_dir, input_msa=cropped_boundary_MSA,
                                                              clean_column_threshold=0.01,
                                                              min_length_num=10, cluster_num=2, cluster_col_thr=500)

                # len(final_msa_consistency) == 1 means clustering is not necessary
                if not final_msa_consistency:
                    break
                else:
                    cluster_bed_files_list, fasta_out_flank_mafft_gap_rm = final_msa_consistency

                # If the clustered MSA contains only 1 element, check if the BED file has the same line number as the original file
                if len(cluster_bed_files_list) == 1:
                    df_new = pd.read_csv(cluster_bed_files_list[0], sep='\t', header=None)

                    if len(df) == len(df_new):
                        break
                    else:
                        bed_file = cluster_bed_files_list[0]
                else:
                    # Only use the top 1 cluster for further analysis
                    bed_file = cluster_bed_files_list[0]
        except Exception as e:
            with open(error_files, "a") as f:
                # Get the traceback content as a string
                tb_content = traceback.format_exc()
                f.write(f"Boundary definition clustering error {seq_name}\n")
                f.write('\n' + tb_content + '\n\n')
                prcyan(f"\nBoundary definition clustering failed for {seq_name} with error:\n{e}")
                prcyan('\n' + tb_content + '\n')
                raise Exception

        msa_loop_n += 1

        # Use smaller extension step for the next round extension
        ex_step_size = 700

    #####################################################################################################
    # Code block: Check if the final MSA contains too many instances of the ambiguous letter "N"
    #####################################################################################################
    try:
        initial_cons = con_generater_no_file(cropped_boundary_MSA, threshold=0.55, ambiguous="N")

        # Calculate the proportion of 'N' in the sequence
        n_proportion = initial_cons.count("N") / len(initial_cons)

        # Check if the proportion is greater than 40%
        # If 40% of this consensus sequence is "N", stop analysis for this MSA
        if n_proportion > 0.4:
            return False

    #####################################################################################################
    # Code block: For LTR elements, check if the cropped MSA starts with given patterns like TGA ACA
    #####################################################################################################

        # If both start and end patterns are 'None', and found_match_crop is 'False' (no terminal repeat found), 
        # skip this block, because terminal repeats can precisely define the start and end points of TEs.
        if start_patterns is not None or end_patterns is not None and not found_match_crop:
            if is_LTR(cropped_alignment_output_file_g):  # Check if file name contains "LTR"
                # Generate consensus sequences
                consensus_seq = con_generater_no_file(cropped_alignment_output_file_g, threshold=0.7, ambiguous="X")

                # Convert consensus_seq to list
                consensus_seq = list(consensus_seq)

                # Four variables will be returned
                start_matched, end_matched, check_start, check_end = check_and_update(
                    consensus_seq, start=cropped_boundary.start_post, end=cropped_boundary.end_post,
                    start_patterns=start_patterns, end_patterns=end_patterns)

                # If the new start or end positions are different from previous MSA, replace with the new 
                # start or end position in the cropped_boundary object and generate the new MSA file
                if cropped_boundary.start_post != check_start or cropped_boundary.end_post != check_end:
                    cropped_boundary.start_post = check_start
                    cropped_boundary.end_post = check_end
                    # Generate the new MSA file based on new start and end positions
                    cropped_boundary_MSA = cropped_boundary.crop_MSA(output_dir, crop_extension=0)
    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"\nStart and end pattern definition failed for {seq_name} with error \n{e}\n")
            f.write('\n' + tb_content + '\n\n')
            prcyan(f"\nStart and end patterns like 'TGT' or 'ACA' definition failed for {seq_name} with error \n{e}")
            prgre("\nThis error will not affect the final result significantly, you can choose to ignore it. "
                  "For traceback text, please refer to 'error_file.txt' in the 'Multiple_sequence_alignment' folder.\n")

    #####################################################################################################
    # Code block: Generate MSA for CIAlign plot
    #####################################################################################################
    try:
        # Add 300 columns to the left of start point
        cropped_boundary_manual_MSA_left = select_window_columns(
            bed_fasta_mafft_boundary_crop_for_select, output_dir, column_mapping[cropped_boundary.start_post],
            "left", window_size=300)

        # Add 300 columns to the right of the end point
        cropped_boundary_manual_MSA_right = select_window_columns(
            bed_fasta_mafft_boundary_crop_for_select, output_dir, column_mapping[cropped_boundary.end_post],
            "right", window_size=300)

        # Concatenate MSA for manual curation
        # The concat_start and concat_end points correspond with cropped_boundary.start_post and cropped_boundary.end_post
        # Because the input file of "concatenate_alignment()" is an alignment object and not a file path, read "cropped_boundary_MSA"
        # as an alignment object.
        input_file_name = str(os.path.basename(bed_fasta_mafft_boundary_crop_for_select)) + "_proof_anno"
        cropped_boundary_MSA_alignment = AlignIO.read(cropped_boundary_MSA, "fasta")
        cropped_boundary_manual_MSA_concatenate, concat_start_man, concat_end_man = concatenate_alignments(
            cropped_boundary_manual_MSA_left, cropped_boundary_MSA_alignment, cropped_boundary_manual_MSA_right,
            input_file_name=input_file_name, output_dir=output_dir
        )

    #####################################################################################################
    # Code block: Concatenate beginning and end part of MSA and plotting
    #####################################################################################################

        # "sequence_len" represent the length of the final cropped MSA
        # Extract the beginning and end columns of the cropped MSA, then join them by "----------" (pseudo-gaps).
        # Return MSA length, which will be used for plotting
        cropped_boundary_plot_select_start_end_and_joint, sequence_len = select_start_end_and_join(
            cropped_alignment_output_file_g, output_dir, cropped_boundary.start_post, cropped_boundary.end_post
        )

        # "column_mapping" is a dictionary, the key represents the cropped nucleotide position and the value represents the
        # original MSA nucleotide position.
        # Add 50 columns to the left of the start point
        cropped_boundary_plot_left = select_window_columns(
            bed_fasta_mafft_boundary_crop_for_select, output_dir, column_mapping[cropped_boundary.start_post], "left",
            window_size=50
        )

        # Add 50 columns to the right of the end point
        cropped_boundary_plot_right = select_window_columns(
            bed_fasta_mafft_boundary_crop_for_select, output_dir, column_mapping[cropped_boundary.end_post], "right",
            window_size=50
        )

        # Concatenate MSA
        # The concat_start and concat_end points correspond with cropped_boundary.start_post and cropped_boundary.end_post
        cropped_boundary_plot_concatenate, concat_start, concat_end = concatenate_alignments(
            cropped_boundary_plot_left, cropped_boundary_plot_select_start_end_and_joint, cropped_boundary_plot_right,
            input_file_name=cropped_boundary_MSA, output_dir=output_dir
        )
    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"\nConcatenate MSA for manual inspection plot failed for {seq_name} with error:\n{e}")
            f.write('\n' + tb_content + '\n\n')
            prcyan(f"\nConcatenate MSA for manual inspection plot failed for {seq_name} with error:\n{e}")
            prcyan('\n' + tb_content + '\n')
            raise Exception

    #####################################################################################################
    # Code block: Predict ORFs and PFAM domains, determine sequence direction
    #####################################################################################################
    try:
        # Generate consensus sequence for ORF prediction
        # Use lower threshold to enable  more PFAM results
        orf_cons = con_generater(cropped_boundary_MSA, output_dir, threshold=0.5)

        # Predict ORF, scan for PFAM domains
        orf_domain_plot_object = PlotPfam(orf_cons, output_dir, pfam_database_dir=pfam_dir, mini_orf=mini_orf,
                                          after_tetrimmer=True)

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
                        input_file=orf_cons, output_file=orf_cons)

                    # Reverse-complement MSA files
                    cropped_boundary_MSA = reverse_complement_seq_file(
                        input_file=cropped_boundary_MSA, output_file=cropped_boundary_MSA
                    )
                    cropped_boundary_manual_MSA_concatenate = reverse_complement_seq_file(
                        input_file=cropped_boundary_manual_MSA_concatenate,
                        output_file=cropped_boundary_manual_MSA_concatenate
                    )

                    cropped_boundary_plot_concatenate = reverse_complement_seq_file(
                        input_file=cropped_boundary_plot_concatenate,
                        output_file=cropped_boundary_plot_concatenate
                    )
                    # Reverse complement input sequence, this will be used for dotplot
                    seq_file_reverse_c_path = os.path.join(output_dir, f"{os.path.basename(seq_name)}.fasta_rc")
                    seq_file_reverse_c = reverse_complement_seq_file(input_file=seq_file,
                                                                     output_file=seq_file_reverse_c_path)
                    reverse_complement = True

                    # Define the new start and end points for cropped_boundary_manual_MSA_concatenate
                    # cropped_boundary_manual_MSA_concatenate will be used for the CIAlign plot
                    # Get MSA sequence length
                    cropped_boundary_manual_MSA_concatenate_align = AlignIO.read(cropped_boundary_manual_MSA_concatenate, "fasta")

                    # get_alignment_length() is a build-in function
                    cropped_boundary_manual_MSA_concatenate_length = cropped_boundary_manual_MSA_concatenate_align.get_alignment_length()

                    # Use intermediate number to get the new start and end number
                    concat_start_man_intermediate = cropped_boundary_manual_MSA_concatenate_length - concat_end_man - 1
                    concat_end_man_intermediate = cropped_boundary_manual_MSA_concatenate_length - concat_start_man - 1
                    concat_start_man = concat_start_man_intermediate
                    concat_end_man = concat_end_man_intermediate

                    # Define the new start and end points for cropped_boundary_plot_concatenate
                    # cropped_boundary_plot_concatenate will be used for the te_trimmer plot
                    # Get MSA sequence length
                    cropped_boundary_plot_concatenate_align = AlignIO.read(cropped_boundary_plot_concatenate, "fasta")
                    cropped_boundary_plot_concatenate_length = cropped_boundary_plot_concatenate_align.get_alignment_length()

                    # Use intermediate number to get the new start and end points
                    concat_start_intermediate = cropped_boundary_plot_concatenate_length - concat_end - 1
                    concat_end_intermediate = cropped_boundary_plot_concatenate_length - concat_start - 1
                    concat_start = concat_start_intermediate
                    concat_end = concat_end_intermediate

                    # Based on the new consensus sequence, predict ORFs and PFAM domains again, then plot
                    orf_domain_plot_object = PlotPfam(orf_cons, output_dir, pfam_database_dir=pfam_dir, mini_orf=mini_orf)
                    if orf_domain_plot_object.run_getorf():
                        # "run_pfam_scan() will return 'True' if any PFAM domains were found. Otherwise, it will return 'False'
                        if_pfam_domain, pfam_result_file = orf_domain_plot_object.run_pfam_scan()
                        if if_pfam_domain:
                            orf_domain_plot = orf_domain_plot_object.orf_domain_plot()

            else:
                # If only ORFs but no PFAM domain were found, plot ORFs
                orf_domain_plot = orf_domain_plot_object.orf_domain_plot()
    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"\nConcatenate MSA for the maunual inspection plot failed for {seq_name} with error:\n{e}")
            f.write('\n' + tb_content + '\n\n')
            prcyan(f"\nORF or PFAM prediction failed for {seq_name} with error:\n{e}")
            prcyan('\n' + tb_content + '\n')
            raise Exception

    #####################################################################################################
    # Code block: Plot multiple sequence alignment. 
    #####################################################################################################
    # Plotting is done after PFAM predictions in case consensus/MSA are in the wrong direction.
    try:
        # Plot MSA, which can easily verify if the start and end crop points are correct
        MSA_plot = process_msa(cropped_boundary_plot_concatenate, output_dir, concat_start, concat_end, sequence_len)

    #####################################################################################################
    # Code block: Plot entire MSA by CIAlign with start and end points
    #####################################################################################################

        # Define MSA plot output file
        cropped_boundary_manual_MSA_concatenate_plot = f"{cropped_boundary_manual_MSA_concatenate}_plot.pdf"

        # Convert MSA to array
        cropped_boundary_manual_MSA_concatenate_array, nams = cialign.FastaToArray(cropped_boundary_manual_MSA_concatenate)

        # Draw the entire MSA
        cialign.drawMiniAlignment(cropped_boundary_manual_MSA_concatenate_array, nams,
                                  cropped_boundary_manual_MSA_concatenate_plot, concat_start_man, concat_end_man)
    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"\nMSA plot failed for {seq_name} with error:\n{e}")
            f.write('\n' + tb_content + '\n\n')
            prcyan(f"\nMSA plot failed {seq_name} with error:\n{e}")
            prgre("\nMSA plots are only used to evaluate TEtrimmer and will not affect the final TE consensus library."
                  " For traceback text, please refer to 'error_file.txt' in the 'Multiple_sequence_alignment' folder\n")

    #####################################################################################################
    # Code block: Generate TE-Aid plot
    #####################################################################################################
    try:
        # so.path.abspath(__file__) will return the current executable Python file
        # TE-Aid package is stored in the same directory as this function file.
        TE_aid_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "TE-Aid-master")

        # Because terminal repeats were found before, use the previous result
        # Run TE_aid. If low_copy is 'True'
        # The low_copy will only affect it to keep self blast result from TEAid.
        # it is Ture, TE_aid_object.run will check the terminal repeat based on the self blast output of TEAid
        # Otherwise, it use the terminal repeat result "found_match_crop"
        if not found_match_crop:
            TE_aid_object = TEAid(orf_cons, output_dir, genome_file, error_file=error_files, TE_aid_dir=TE_aid_path)
            TE_aid_plot, found_match = TE_aid_object.run(low_copy=True)
        else:
            # low_copy=False will prevent self-BLAST and instead check for terminal repeats
            TE_aid_object = TEAid(orf_cons, output_dir, genome_file, error_file=error_files, TE_aid_dir=TE_aid_path)

            # found_match returns 'False' in this case
            TE_aid_plot, found_match = TE_aid_object.run(low_copy=False)

            # Assign found_match_crop ("LTR" or "TIR") to found_match for further analysis
            found_match = found_match_crop

        # Run TE_aid to plot the query sequence, if required. Because one query file can return multiple
        # clusters, TE-Aid will check if a TE-Aid plot has been generated before.
        if plot_query:
            query_file = seq_obj.get_input_fasta()
            TE_aid_object_query = TEAid(query_file, output_dir, genome_file, error_file=error_files,
                                        TE_aid_dir=TE_aid_path)
            TE_aid_plot_query, found_match_query = TE_aid_object_query.run(low_copy=False, label=False)
        else:
            TE_aid_plot_query = None
    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"\nTE-Aid plot failed for {seq_name} with error:\n{e}")
            f.write('\n' + tb_content + '\n\n')
            prcyan(f"\nTE-Aid plot failed for {seq_name} with error:\n{e}")
            prcyan('\n' + tb_content + '\n')
            raise Exception

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
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"\nDotplot failed for {seq_name} with error:\n{e}")
            f.write('\n' + tb_content + '\n\n')

    try:
        if dotplot_pdf is not None:
            # Because dotmatcher can't change output size, scale it up to make it more clear in the merged pdf
            scale_dotplot_pdf = scale_single_page_pdf(dotplot_pdf, f"{dotplot_pdf}_su.pdf", scale_ratio=2)
            dotplot_pdf = scale_dotplot_pdf
    except Exception:
        with open(error_files, "a") as f:
            # This is not mandatory, skip if any error occurred
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"\nps file to pdf conversion failed for {seq_name} with error:\n{e}")
            f.write('\n' + tb_content + '\n\n')
            
    #####################################################################################################
    # Code block: Merge plot files
    #####################################################################################################
    try:
        merged_pdf_path = merge_pdfs(output_dir, os.path.basename(cropped_boundary_MSA),
                                     MSA_plot, cropped_boundary_manual_MSA_concatenate_plot,
                                     TE_aid_plot, TE_aid_plot_query, orf_domain_plot, input_orf_pfam, dotplot_pdf)

    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"\nPlot merging to PDF failed for {seq_name} with error:\n{e}")
            f.write('\n' + tb_content + '\n\n')
            prcyan(f"\nPlot merging to PDF failed for {seq_name} with error:\n{e}")
            prcyan('\n' + tb_content + '\n')
            raise Exception

    #####################################################################################################
    # Code block: Update sequence object (seq_obj)
    #####################################################################################################
    try:
        # Create a folder in the same directory as output_dir to store annotation files
        parent_output_dir = os.path.dirname(output_dir)

        # Define different levels of proof_curation folder
        perfect_proof = os.path.join(proof_curation_dir, "Annotations_perfect")
        good_proof = os.path.join(proof_curation_dir, "Annotations_good")
        intermediate_proof = os.path.join(proof_curation_dir, "Annotations_check_recommended")
        need_check_proof = os.path.join(proof_curation_dir, "Annotations_check_required")

        # Create the directory if it does not exist
        os.makedirs(proof_curation_dir, exist_ok=True)
        os.makedirs(perfect_proof, exist_ok=True)
        os.makedirs(good_proof, exist_ok=True)
        os.makedirs(intermediate_proof, exist_ok=True)
        os.makedirs(need_check_proof, exist_ok=True)

        # Define temporary classified and unknown final consensus file, which will be used for
        # reclassification by RepeatMasker
        final_unknown_con_file = os.path.join(classification_dir, "temp_TEtrimmer_unknown_consensus.fasta")
        final_classified_con_file = os.path.join(classification_dir, "temp_TEtrimmer_classified_consensus.fasta")

        # Define unique sequence names
        # Because seq_obj.create_consi_obj() is generated after the unique name definition, len(seq_obj.consi_obj_list)
        # is the appended number for the unique name.
        consi_n = len(seq_obj.consi_obj_list)

        if consi_n > 0:
            uniq_seq_name = f"{seq_name}_{consi_n:02}"
        else:
            uniq_seq_name = seq_name

        # Create consensus object
        consi_obj = seq_obj.create_consi_obj(uniq_seq_name)

        # Generate final consensus sequence
        # Compared with initial consensus, the final consensus sequence can have a different orientation.
        # For this reason, generate consensus again
        final_con = con_generater_no_file(cropped_boundary_MSA, threshold=cons_threshold)
        sequence = str(final_con).upper()

        # Storing consensus sequence length into consi_obj
        consi_obj.set_new_length(len(sequence))

        # Store consensus sequence into consi_obj not necessary
        # consi_obj.set_cons_seq(sequence)

        # Store MSA sequence number into consi_obj
        MSA_for_final_cons = AlignIO.read(cropped_boundary_MSA, "fasta")
        MSA_for_final_cons_seq_n = len(MSA_for_final_cons)
        consi_obj.set_cons_MSA_n(MSA_for_final_cons_seq_n)

        # Store terminal repeat to consi_obj
        if found_match == "LTR" or found_match == "TIR":
            consi_obj.set_new_terminal_repeat(found_match)

        # Store full length sequence number from BLAST search into consi_obj
        # check_blast_full_n is a function in TE_Aid class. TEtrimmer will use the blast result of TE Aid
        blast_full_length_n = TE_aid_object.check_blast_full_n(consi_obj, engine=engine)
        consi_obj.set_blast_full_n(blast_full_length_n)

        # Store PFAM predictions to consi_obj
        if if_pfam_domain:
            consi_obj.set_cons_pfam(if_pfam_domain)

    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"\nUpdate sequence object failed for {seq_name} with error:\n{e}")
            f.write('\n' + tb_content + '\n\n')
            prcyan(f"\nUpdate sequence object failed for {seq_name} with error:\n{e}")
            prcyan('\n' + tb_content + '\n')
            raise Exception

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
        if classify_all or (
                classify_unknown and (seq_obj.check_unknown() or (left_ex_total + right_ex_total >= 4000))):
            # Define different folders for each sequence
            # the suffix .fasta is important, this ensures that this folder can be deleted later
            classification_seq_folder = os.path.join(classification_dir, f"{uniq_seq_name}.fasta")
            os.makedirs(classification_seq_folder, exist_ok=True)

            # Define consensus file path used for classification
            classification_seq_file = os.path.join(classification_seq_folder, uniq_seq_name)

            with open(classification_seq_file, "w") as f:
                # RepeatClassifier input cannot be a single sequence, add >Dummy to enable RepeatClassifier run
                f.write(">" + uniq_seq_name + "\n" + sequence + "\n" + ">Dummy" + "\n" + "T" + "\n")

            TE_type = classify_single(classification_seq_file)

            # Only update new_TE_type if classify_single was successful
            if TE_type:
                # Set TE_type after RepeatClassifier
                consi_obj.set_new_TE_type(TE_type)

            # Clean classification folder if debug is 'True'
            if not debug:
                remove_files_with_start_pattern(classification_dir, f"{uniq_seq_name}.fasta")
    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"\nRepeatClassifier classification error for {seq_name}")
            f.write('\n' + tb_content + '\n\n')
        prcyan(f"\nNote: RepeatClassifier doesn't work for {seq_name} with error {e}")
        prgre("\nThis won't affect final TE consensus sequences but only the classification. You can choose to ignore this. "
              "For traceback text, please refer to 'error_file.txt' under 'Multiple_sequence_alignment' folder\n")

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
    try:
        if (consi_obj.new_TE_terminal_repeat != "False" and
                consi_obj.new_TE_type != "NaN" and "unknown" not in consi_obj.new_TE_type.lower() and
                consi_obj.new_TE_MSA_seq_n >= 30 and
                consi_obj.new_TE_blast_full_length_n >= 5 and
                consi_obj.cons_pfam):
            consi_obj.set_cons_evaluation("Perfect")

        elif (consi_obj.new_TE_terminal_repeat != "False" and
              consi_obj.new_TE_MSA_seq_n >= 10 and
              consi_obj.new_TE_blast_full_length_n >= 2):
            consi_obj.set_cons_evaluation("Good")

        elif (consi_obj.new_TE_MSA_seq_n >= 20 and
              consi_obj.new_TE_blast_full_length_n >= 2):
            consi_obj.set_cons_evaluation("Reco_check")

        else:
            consi_obj.set_cons_evaluation("Need_check")

    #####################################################################################################
    # Code block: Code block: Move file for manual inspection
    #####################################################################################################

        consi_obj.set_proof_curation_file()

        # modify fasta_out_flank_mafft_gap_rm fasta header based on the bed file, this can allow the
        # extension function in the final GUI.
        # For example: change 1(+) to scaffold_1:23256-24757(+)
        cropped_boundary_MSA_nm = modify_fasta_headers(bed_out_flank_file, cropped_boundary_MSA)
        bed_fasta_mafft_boundary_crop_for_select_nm = modify_fasta_headers(bed_out_flank_file,
                                                                           bed_fasta_mafft_boundary_crop_for_select)

        # Define file name for inspection file
        file_copy_pattern = [
            (merged_pdf_path, str(consi_obj.proof_pdf)),
            (cropped_boundary_MSA_nm, str(consi_obj.proof_fasta)),
            (bed_fasta_mafft_boundary_crop_for_select_nm, str(consi_obj.proof_raw)),
            (cluster_msa, str(consi_obj.proof_cluster))
        ]

        files_moved_successfully = True

        for pattern, new_name in file_copy_pattern:
            try:
                if consi_obj.cons_evaluation == "Perfect":
                    destination_dir = perfect_proof
                elif consi_obj.cons_evaluation == "Good":
                    destination_dir = good_proof
                elif consi_obj.cons_evaluation == "Reco_check":
                    destination_dir = intermediate_proof
                else:
                    destination_dir = need_check_proof

                # Copy the file to the new location with the new unique name
                if pattern:
                    shutil.copy(pattern, os.path.join(destination_dir, new_name))

            except Exception as e:
                with open(error_files, "a") as f:
                    # Get the traceback content as a string
                    tb_content = traceback.format_exc()
                    f.write(f"Copy file error for {pattern}\n")
                    f.write(tb_content + '\n\n')
                click.echo(f"Error copying {pattern} to {new_name}: {e}")
                files_moved_successfully = False

        if hmm:  # Generate HMM files
            consi_obj.set_hmm_file()
            hmm_output_file = os.path.join(hmm_dir, consi_obj.hmm_file)
            generate_hmm_from_msa(cropped_boundary_MSA, hmm_output_file, error_files)

        # Classification of unknown consensus TEs will be attempted again later by using successfully classified sequences
        if "Unknown" in updated_TE_type:
            with open(final_unknown_con_file, "a") as f:  # 'a' mode for appending
                f.write(">" + uniq_seq_name + "\n" + sequence + "\n")
        else:
            with open(final_classified_con_file, "a") as f:
                f.write(">" + uniq_seq_name + "#" + updated_TE_type + "\n" + sequence + "\n")

        # Write all consensus sequence to final_cons_file.
        with open(final_con_file, "a") as f:
            f.write(">" + uniq_seq_name + "#" + updated_TE_type + "\n" + sequence + "\n")

        # Write all consensus sequences to final_cons_file_no_low_copy.
        with open(final_con_file_no_low_copy, "a") as f:
            f.write(">" + uniq_seq_name + "#" + updated_TE_type + "\n" + sequence + "\n")

        return files_moved_successfully

    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"\nMoving of files failed for {seq_name} with error:\n{e}")
            f.write('\n' + tb_content + '\n\n')
            prcyan(f"\nMoving of files failed for {seq_name} with error:\n{e}")
            prcyan('\n' + tb_content + '\n')
            raise Exception




