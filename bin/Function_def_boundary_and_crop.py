from Class_define_boundary import DefineBoundary
from Class_crop_end import CropEnd
from Class_crop_end_by_gap import CropEndByGap
from Class_blast_extension_mafft import SequenceManipulator
from Class_select_ditinct_columns import CleanAndSelectColumn
from Class_MSA_plot import MSAPainter
from Class_check_start_end import StartEndChecker
from Class_TE_aid import TEAid
from Class_orf_domain_prediction import PlotPfam, determine_sequence_direction
from Function_clean_and_clauster_MSA import clean_and_cluster_MSA
from Bio import AlignIO
import Cialign_plot
import os
import click
import shutil
import pandas as pd
from PyPDF2 import PdfMerger
from Classification import classify_single
import traceback


def check_self_alignment(seq_obj, seq_file, output_dir, genome_file, blast_hits_count, blast_out_file):

    blast_full_length_n = check_blast_low_copy(seq_obj, blast_out_file, identity=90, coverage=0.9, min_hit_length=100)

    # At least 2 lines need to meet the requirement
    if blast_full_length_n >= 2:

        check_blast = True

        # Check self-alignment of terminal repeats
        TE_aid_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "TE-Aid-master")

        # Run TE_aid
        TE_aid_object = TEAid(seq_file, output_dir, genome_file, TE_aid_dir=TE_aid_path)
        TE_aid_plot, found_match = TE_aid_object.run(low_copy=True)
        seq_obj.update_blast_hit_n(blast_hits_count)

        # Convert found_match to True when LTR or TIR is found
        if found_match == "LTR" or found_match == "TIR":
            found_match = True

        # Update_low_copy return true when check_80 and found_match are both true
        check_low_copy = seq_obj.update_low_copy(check_blast, found_match)

    else:
        check_low_copy = False
    return check_low_copy, blast_full_length_n


def check_blast_low_copy(seq_obj, blast_out_file, identity=90, coverage=0.9, min_hit_length=100, te_aid_blast=False):

    # The TE Aid blast output have a header
    if te_aid_blast:
        df = pd.read_csv(blast_out_file, sep="\s+")
    else:
        df = pd.read_csv(blast_out_file, sep="\s+", header=None)

    identity_condition = df[2] > identity
    coverage_condition = df[3] / seq_obj.old_length > coverage
    length_condition = df[3] > min_hit_length

    # Filter the DataFrame
    filtered_df = df[identity_condition & coverage_condition & length_condition]

    blast_full_length_n = filtered_df.shape[0]

    return blast_full_length_n


def crop_end_and_clean_column(input_file, output_dir, crop_end_threshold=16, window_size=20, gap_threshold=0.8):

    # Window_size means the checked nucleotide number each time
    # Threshold means the sum of nucleotide proportion in window_size must greater than that
    cropped_MSA = CropEnd(input_file, threshold=crop_end_threshold, window_size=window_size)
    cropped_MSA.pro_calculation()
    cropped_MSA.find_positions()
    cropped_MSA.crop_alignment()

    # write_to_file() function will return absolute cropped alignment file
    cropped_MSA_output_file = cropped_MSA.write_to_file(output_dir)

    # Remove gaps again after crop end step
    cropped_alignment_output_file_with_gap = SequenceManipulator()

    # "column_mapping" is a dictionary the key is gap removed MSA index, the value is the corresponded original
    # MSA index
    cropped_alignment_output_file_no_gap, column_mapping = cropped_alignment_output_file_with_gap.remove_gaps(
        cropped_MSA_output_file, output_dir, threshold=gap_threshold, min_nucleotide=5)

    return cropped_alignment_output_file_no_gap, column_mapping


# Because for each query file, several subgroups of multiple sequence alignment will be generated,
# This function will name them differently
def get_unique_filename(base_path, original_name):
    counter = 0
    name, ext = os.path.splitext(original_name)

    # If the extension is .anno_fasta, change it to .anno.fasta
    # This will make sure the right order of the files
    if ext == ".anno_fasta":
        ext = ".anno.fasta"

    # Initialize with the original name
    new_name = f"{name}{ext}"

    # While the constructed name already exists, increment the counter and reconstruct the name
    while os.path.exists(os.path.join(base_path, new_name)):
        counter += 1
        new_name = f"{name}_{counter:02}{ext}"  # Adjusted to use zero-padding

    return new_name


def find_boundary_and_crop(bed_file, genome_file, output_dir, pfam_dir, seq_obj, hmm,
                           cons_threshold=0.8, ext_threshold=0.7, ex_step_size=1000, max_extension=7000,
                           gap_threshold=0.4, gap_nul_thr=0.7, crop_end_thr=16, crop_end_win=20,
                           crop_end_gap_thr=0.1, crop_end_gap_win=150, start_patterns=None, end_patterns=None,
                           mini_orf=200, define_boundary_win=150):
    """
    :param bed_file: str, bed file directory
    :param genome_file: str, genome directory
    :param output_dir: str, output file directory
    :param pfam_dir: str, pfam database directory
    :param cons_threshold: num 0 to 1 default: 0.8, threshold used for final consensus sequence generation
    :param ext_threshold: num 0 to 1 default: 0.7, threshold used for define the extension extent. The higher number
    represent more stringent extension threshold (short extension in total), the smaller number means it become easier
    to have a final longer extension for each sider of the sequence.
    :param ex_step_size: num default: 1000, number of nucleotide will be added to left or right side each time
    :param max_extension: num default: 7000, the maximum extension number for right and left side
    :param crop_end_thr: num default: 16, when the ambiguous letter more than this number, delete this window
    :param crop_end_win: num default: 20, window size used for crop end process
    :param gap_threshold: num default: 0.4, columns with greater gap proportion than gap_threshold and the most common
    nucleotide proportion in this column is less than gap_nul_thr will be removed
    :param gap_nul_thr: num default: 0.7, set nucleotide proportion to decide if remove this column
    :param crop_end_gap_thr: num default: 0.05, set gap threshold for crop end by gap
    :param crop_end_gap_win: num default: 300, set window size used to crop end by gap
    :param start_patterns: str default: None, patterns to check for start points
    :param end_patterns: str default: None, patterns to check for end points
    :param min_orf: num default: 200, set minimum orf length for orf prediction

    """
    # Check if this is a LINE element, if so decrease the ext_threshold number,
    # because LINE have higher divergence at 5' end
    seq_name = seq_obj.name
    if"LINE" in seq_obj.old_TE_type:
        ext_threshold = ext_threshold - 0.2
    if_left_ex = True
    if_right_ex = True

    #####################################################################################################
    # Code block: Define left and right sides extension number
    #####################################################################################################

    left_ex = ex_step_size
    right_ex = ex_step_size
    final_MSA_consistent = False
    intact_loop_times = 1

    # intact_loop_times is the iteration numbers, <3 means iterate 2 times.
    while (if_left_ex or if_right_ex) and (left_ex <= max_extension and right_ex <= max_extension) \
            and not final_MSA_consistent and intact_loop_times < 4:

        bed = SequenceManipulator()
        # bedtools will makesure the extension won't excess the maximum length of that chromosome
        bed_fasta, bed_out_flank_file = bed.extract_fasta(bed_file, genome_file, output_dir, left_ex, right_ex)

        # align_sequences() will return extended MSA absolute file
        bed_fasta_mafft_with_gap = bed.align_sequences(bed_fasta, output_dir)

        if not os.path.isfile(bed_fasta_mafft_with_gap):
            click.echo(f"{bed_file} has problem during mafft extension step")
            return False

        # Remove nucleotide whose proportion is smaller than threshold
        bed_fasta_mafft_with_gap_column_clean_object = CleanAndSelectColumn(bed_fasta_mafft_with_gap, threshold=0.08)
        bed_fasta_mafft_with_gap_column_clean = bed_fasta_mafft_with_gap_column_clean_object.clean_column(output_dir)

        # Remove gaps with similarity check
        bed_fasta_mafft = bed.remove_gaps_with_similarity_check(bed_fasta_mafft_with_gap_column_clean, output_dir,
                                                                gap_threshold=0.6, simi_check_gap_thre=0.4,
                                                                similarity_threshold=0.7,
                                                                min_nucleotide=5
                                                                )

        bed_fasta_mafft_object = CropEndByGap(bed_fasta_mafft, gap_threshold=crop_end_gap_thr,
                                              window_size=crop_end_gap_win)

        bed_fasta_mafft_cop_end_gap = bed_fasta_mafft_object.write_to_file(output_dir)

        # Threshold to generate consensus sequence
        bed_boundary = DefineBoundary(bed_fasta_mafft_cop_end_gap, threshold=ext_threshold,
                                      check_window=define_boundary_win, max_X=0.3, if_con_generater=False)

        if bed_boundary.if_continue:
            if bed_boundary.left_ext:  # boundary.left_ext will become true when more extension is required
                left_ex += ex_step_size
            else:
                if_left_ex = bed_boundary.left_ext

            if bed_boundary.right_ext:
                right_ex += ex_step_size
            else:
                if_right_ex = bed_boundary.right_ext
        elif not bed_boundary.if_continue:

            return "Short_sequence"

        #####################################################################################################
        # Code block: Remove gaps and define the start and end position
        #####################################################################################################

        if (not if_left_ex and not if_right_ex) or left_ex == max_extension or right_ex == max_extension:

            # Remove gap block with similarity check
            # gap_threshold=0.8 means if gap proportion is greater than 80% this column will be regarded as gap column directly
            # without further nucleotide similarity check
            bed_fasta_mafft_gap_block_sim = bed.remove_gaps_block_with_similarity_check(
                bed_fasta_mafft_with_gap_column_clean, output_dir, gap_threshold=0.8, simi_check_gap_thre=gap_threshold,
                similarity_threshold=gap_nul_thr, conservation_threshold=0.5)

            bed_boundary = DefineBoundary(bed_fasta_mafft_gap_block_sim, threshold=ext_threshold,
                                          check_window=define_boundary_win, max_X=0.2, if_con_generater=False)
            bed_fasta_mafft_boundary_crop = bed_boundary.crop_MSA(output_dir, crop_extension=300)

            # Because for LINE element bed_fasta_mafft_boundary_crop will be changed. Copy it to the other variable
            bed_fasta_mafft_boundary_crop_for_select = bed_fasta_mafft_boundary_crop

            if "LINE" in bed_file:
                # For the high divergence region, more gaps can be found. According to this feature, remove high divergence
                # region this function is very useful for dealing with LINE elements
                cropped_MSA_by_gap = CropEndByGap(bed_fasta_mafft_boundary_crop, gap_threshold=crop_end_gap_thr,
                                                  window_size=crop_end_gap_win)

                bed_fasta_mafft_boundary_crop = cropped_MSA_by_gap.write_to_file(output_dir)

            # Gaps are removed again after crop end process
            cropped_alignment_output_file_no_gap, column_mapping = crop_end_and_clean_column(bed_fasta_mafft_boundary_crop,
                                                                                               output_dir,
                                                                                               crop_end_threshold=crop_end_thr,
                                                                                               window_size=crop_end_win,
                                                                                               gap_threshold=0.8)

            # Crop end can't define the final boundary, use DefineBoundary again to define start position
            cropped_boundary = DefineBoundary(cropped_alignment_output_file_no_gap, threshold=0.8, check_window=3, max_X=0)
            cropped_boundary_MSA = cropped_boundary.crop_MSA(output_dir, crop_extension=0)

            #####################################################################################################
            # Code block: Check the consistency of the final MSA
            #####################################################################################################
            final_MSA_consistency = clean_and_cluster_MSA(cropped_boundary_MSA, bed_out_flank_file, output_dir,
                                                          clean_column_threshold=0.08,
                                                          min_length_num=10, cluster_num=2, cluster_col_thr=250
                                                          )

            # False means the sequences number in each cluster is smaller than minimum number, normally 10
            # True means not necessary to do further cluster
            if final_MSA_consistency is False or final_MSA_consistency is True:

                final_MSA_consistent = True



            # else means that further clustering is required
            else:

                # Only use the first cluster for further analysis, which contains the most sequences
                new_bed_file = final_MSA_consistency[0]
                bed_file = new_bed_file

                intact_loop_times = intact_loop_times + 1

                # Reset extension parameters to enable whole while loop
                if_left_ex = True
                if_right_ex = True

                left_ex = 0
                right_ex = 0

    #####################################################################################################
    # Code block: Check if the final MSA contains too many ambiguous letter "N"
    #####################################################################################################

    initial_cons_object = SequenceManipulator()
    initial_cons = initial_cons_object.con_generater_no_file(cropped_boundary_MSA, threshold=0.7, ambiguous="N")

    # Calculate the proportion of 'N' in the sequence
    n_proportion = initial_cons.count("N") / len(initial_cons)

    # Check if the proportion is greater than 30%
    # If 30% of this consensus sequence is "N", stop analysis for this MSA
    if n_proportion > 0.3:
        return False

    #####################################################################################################
    # Code block: For LTR element, check if the cropped MSA starts with give patterns like TGA ACA
    #####################################################################################################

    # When both start and end patterns are None, skip this block
    if start_patterns is not None or end_patterns is not None:

        # Create StartEndChecker object
        check_start_end_object = StartEndChecker(cropped_alignment_output_file_no_gap,
                                                 start=cropped_boundary.start_post, end=cropped_boundary.end_post,
                                                 start_patterns=start_patterns, end_patterns=end_patterns,
                                                 threshold=0.7)

        if check_start_end_object.check_LTR():  # Check if file name contain "LTR"
            check_start_end_object.con_generater()  # Generate consensus sequences

            # Four variables will be returned,
            start_matched, end_matched, check_start, check_end = check_start_end_object.check_and_update()

            # When the new start or end positions are different with previous, convey the new start or end position
            # to cropped_boundary object and generate the new MSA file
            if cropped_boundary.start_post != check_start or cropped_boundary.end_post != check_end:
                cropped_boundary.start_post = check_start
                cropped_boundary.end_post = check_end
                # Generate the new MSA file based on new start and end positions
                cropped_boundary_MSA = cropped_boundary.crop_MSA(output_dir, crop_extension=0)

    #####################################################################################################
    # Code block: Generate MSA for CIAlign plot
    #####################################################################################################

    cropped_boundary_manual_MSA_object = SequenceManipulator()

    # Get 300 columns left the start point
    cropped_boundary_manual_MSA_left = cropped_boundary_manual_MSA_object.select_window_columns(
        bed_fasta_mafft_boundary_crop_for_select, output_dir, column_mapping[cropped_boundary.start_post], "left", window_size=300
    )

    # Get 300 columns right the end point
    cropped_boundary_manual_MSA_right = cropped_boundary_manual_MSA_object.select_window_columns(
        bed_fasta_mafft_boundary_crop_for_select, output_dir, column_mapping[cropped_boundary.end_post], "right", window_size=300
    )

    # Concatenate MSA for manual curation
    # The concat_start and concat_end points correspond with cropped_boundary.start_post and cropped_boundary.end_post
    # Because the input file of "concatenate_alignment()" is alignment object not file path, read "cropped_boundary_MSA"
    # to an alignment object.
    input_file_name = str(os.path.basename(bed_fasta_mafft_boundary_crop)) + "_proof_anno"
    cropped_boundary_MSA_alignment = AlignIO.read(cropped_boundary_MSA, "fasta")
    cropped_boundary_manual_MSA_concatenate, concat_start_man, concat_end_man = cropped_boundary_manual_MSA_object.concatenate_alignments(
        cropped_boundary_manual_MSA_left, cropped_boundary_MSA_alignment, cropped_boundary_manual_MSA_right,
        input_file_name=input_file_name, output_dir=output_dir
    )

    #####################################################################################################
    # Code block: Concatenate beginning and end part of MSA and plotting
    #####################################################################################################

    # "sequence_len" represent the length of final cropped MSA
    # Extract the beginning and end columns of cropped MSA then join them by "----------".
    # Return MSA length, which will be used for plotting
    cropped_boundary_plot_object = SequenceManipulator()
    cropped_boundary_plot_select_start_end_and_joint, sequence_len = cropped_boundary_plot_object.select_start_end_and_join(
        cropped_alignment_output_file_no_gap, output_dir, cropped_boundary.start_post, cropped_boundary.end_post
    )

    # "column_mapping" is a dictionary, the key represent cropped nucleotide position. The value represent the
    # original MSA nucleotide position.
    # Get 50 columns left the start point
    cropped_boundary_plot_left = cropped_boundary_plot_object.select_window_columns(
        bed_fasta_mafft_boundary_crop_for_select, output_dir, column_mapping[cropped_boundary.start_post], "left",
        window_size=50
    )

    # Get 50 columns right the end point
    cropped_boundary_plot_right = cropped_boundary_plot_object.select_window_columns(
        bed_fasta_mafft_boundary_crop_for_select, output_dir, column_mapping[cropped_boundary.end_post], "right",
        window_size=50
    )

    # Concatenate MSA
    # The concat_start and concat_end points correspond with cropped_boundary.start_post and cropped_boundary.end_post
    cropped_boundary_plot_concatenate, concat_start, concat_end = cropped_boundary_plot_object.concatenate_alignments(
        cropped_boundary_plot_left, cropped_boundary_plot_select_start_end_and_joint, cropped_boundary_plot_right,
        input_file_name=cropped_boundary_MSA, output_dir=output_dir
    )

    #####################################################################################################
    # Code block: Predict ORF and Pfam domain, determine sequence direction
    #####################################################################################################

    # Generate consensus sequence for ORF prediction
    orf_cons_object = SequenceManipulator()

    # Use lower threshold to enable give more pfam results
    orf_cons = orf_cons_object.con_generater(cropped_boundary_MSA, output_dir, threshold=0.5)

    # Predict orf, scan with Pfam
    orf_domain_plot_object = PlotPfam(orf_cons, output_dir, pfam_database_dir=pfam_dir, mini_orf=mini_orf)

    # "run_getorf()" function will return True when any ORF are detected. Otherwise, it will be False
    orf_domain_plot = None
    if orf_domain_plot_object.run_getorf():

        # "run_pfam_scan()" will return True when pfam domains are found, otherwise it will return False
        if_pfam_domain, pfam_result_file = orf_domain_plot_object.run_pfam_scan()
        if if_pfam_domain:

            # "determine_sequence_direction()" will return True when the direction is right, otherwise, it will be False
            if determine_sequence_direction(pfam_result_file):

                # When the direction is right, plot the orf and pfam directly
                orf_domain_plot = orf_domain_plot_object.orf_domain_plot()
            else:
                # When the direction is wrong, reverse complement the corresponded sequence and MSA
                # Reverse complement consensus sequence
                determine_direction_object = SequenceManipulator()

                # The reverse complemented file will overwrite the old file
                orf_cons = determine_direction_object.reverse_complement_seq_file(
                    input_file=orf_cons, output_file=orf_cons)

                # Reverse complement MSA files
                cropped_boundary_MSA = determine_direction_object.reverse_complement_seq_file(
                    input_file=cropped_boundary_MSA, output_file=cropped_boundary_MSA
                )
                cropped_boundary_manual_MSA_concatenate = determine_direction_object.reverse_complement_seq_file(
                    input_file=cropped_boundary_manual_MSA_concatenate,
                    output_file=cropped_boundary_manual_MSA_concatenate
                )

                cropped_boundary_plot_concatenate = determine_direction_object.reverse_complement_seq_file(
                    input_file=cropped_boundary_plot_concatenate,
                    output_file=cropped_boundary_plot_concatenate
                )

                # Define the new start and end points for cropped_boundary_manual_MSA_concatenate
                # cropped_boundary_manual_MSA_concatenate will be used for CIAlign plot
                # Get MSA sequence length
                cropped_boundary_manual_MSA_concatenate_align = AlignIO.read(cropped_boundary_manual_MSA_concatenate, "fasta")

                # get_alignment_length() is a build in function
                cropped_boundary_manual_MSA_concatenate_length = cropped_boundary_manual_MSA_concatenate_align.get_alignment_length()

                # Use intermediate number to get the new start and end number
                concat_start_man_intermediate = cropped_boundary_manual_MSA_concatenate_length - concat_end_man - 1
                concat_end_man_intermediate = cropped_boundary_manual_MSA_concatenate_length - concat_start_man - 1
                concat_start_man = concat_start_man_intermediate
                concat_end_man = concat_end_man_intermediate

                # Define the new start and end points for cropped_boundary_plot_concatenate
                # cropped_boundary_plot_concatenate will be used for te_trimmer plot
                # Get MSA sequence length
                cropped_boundary_plot_concatenate_align = AlignIO.read(cropped_boundary_plot_concatenate, "fasta")
                cropped_boundary_plot_concatenate_length = cropped_boundary_plot_concatenate_align.get_alignment_length()

                # Use intermediate number to get the new start and end number
                concat_start_intermediate = cropped_boundary_plot_concatenate_length - concat_end - 1
                concat_end_intermediate = cropped_boundary_plot_concatenate_length - concat_start - 1
                concat_start = concat_start_intermediate
                concat_end = concat_end_intermediate

                # Based on the new consensus sequence, predict orf and pfam domain again. At the end plot
                orf_domain_plot_object = PlotPfam(orf_cons, output_dir, pfam_database_dir=pfam_dir, mini_orf=mini_orf)
                if orf_domain_plot_object.run_getorf():
                    # "run_pfam_scan()" will return when pfam domains are found, otherwise it will return False
                    if_pfam_domain, pfam_result_file = orf_domain_plot_object.run_pfam_scan()
                    if if_pfam_domain:
                        orf_domain_plot = orf_domain_plot_object.orf_domain_plot()

        else:
            # When only ORF is found but Pfam domain isn't found, plot ORF
            orf_domain_plot = orf_domain_plot_object.orf_domain_plot()

    #####################################################################################################
    # Code block: Plot multiple sequence alignment
    #####################################################################################################

    # Plot MSA, which can easily verify if the start and end crop points are right
    MSA_plot_object = MSAPainter(cropped_boundary_plot_concatenate, output_dir, sequence_len)
    MSA_plot = MSA_plot_object.process(concat_start, concat_end)

    #####################################################################################################
    # Code block: Plot whole MSA by CIAlign style with start and end point
    #####################################################################################################

    # Define while MSA plot output file
    cropped_boundary_manual_MSA_concatenate_plot = f"{cropped_boundary_manual_MSA_concatenate}_plot.pdf"

    # Convert MSA to array
    cropped_boundary_manual_MSA_concatenate_array, nams = Cialign_plot.FastaToArray(cropped_boundary_manual_MSA_concatenate)

    # Draw the whole MSA
    Cialign_plot.drawMiniAlignment(cropped_boundary_manual_MSA_concatenate_array,
                                   nams, cropped_boundary_manual_MSA_concatenate_plot, concat_start_man, concat_end_man)

    #####################################################################################################
    # Code block: Generate TE-Aid plot
    #####################################################################################################

    # so.path.abspath(__file__) will return the current executable python file
    # TE-Aid package is stored at the same directory as this function file.
    TE_aid_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "TE-Aid-master")

    # Run TE_aid
    TE_aid_object = TEAid(orf_cons, output_dir, genome_file, TE_aid_dir=TE_aid_path)
    TE_aid_plot, found_match = TE_aid_object.run(low_copy=True)

    #####################################################################################################
    # Code block: Merge plot files
    #####################################################################################################

    # Code block: Merge plot files
    merger = PdfMerger()

    # List with your PDF file paths. Ensure they are not None before adding them to the list.
    pdf_files = [pdf for pdf in [MSA_plot, cropped_boundary_manual_MSA_concatenate_plot,
                                 TE_aid_plot, orf_domain_plot] if pdf is not None]

    valid_pdf_count = 0  # Counter to keep track of valid PDFs added

    # Iterate over the list of the file paths
    for pdf_file in pdf_files:
        # Check if the file exists before appending
        if os.path.exists(pdf_file):
            # Append PDF files
            merger.append(pdf_file)
            valid_pdf_count += 1
        else:
            click.echo(f"Note: {os.path.basename(pdf_file)} does not exist and will not be merged. This won't affect "
                       f"the final result")

    # Write out the merged PDF file only if there's at least one valid PDF appended
    if valid_pdf_count > 0:
        merged_pdf_path = os.path.join(output_dir, f"{os.path.basename(cropped_boundary_MSA)}_me_plot.pdf")
        merger.write(merged_pdf_path)

    merger.close()

    #####################################################################################################
    # Code block: Move files
    #####################################################################################################

    # Create a folder at the same directory with output_dir to store proof annotation files
    parent_output_dir = os.path.dirname(output_dir)

    # Define final consensus file
    final_con_file = os.path.join(parent_output_dir, "TE_Trimmer_consensus.fasta")

    # Define proof_annotation folder path
    proof_annotation_dir = os.path.join(parent_output_dir, "TE_Trimmer_for_proof_annotation")

    # Define low copy folder
    low_copy_dir = os.path.join(proof_annotation_dir, "Low_copy_TE")
    
    # Construct the path for the new folder
    classification_dir = os.path.join(parent_output_dir, "Classification")

    # Create the directory if it doesn't exist
    if not os.path.exists(proof_annotation_dir):
        os.makedirs(proof_annotation_dir)

    if not os.path.exists(low_copy_dir):
        os.makedirs(low_copy_dir)

    if not os.path.exists(classification_dir):
        os.makedirs(classification_dir)

    # Define proof annotation file name
    file_copy_pattern = [
        (merged_pdf_path, f"{seq_name}.pdf"),
        (cropped_boundary_MSA, f"{seq_name}.fasta"),  # Changed to append .fasta
        (bed_fasta_mafft_boundary_crop_for_select, f"{seq_name}.anno_fasta")
    ]

    files_moved_successfully = True

    for pattern, new_name in file_copy_pattern:

        # Determine the final unique name before copying, add 01, 02, 03 .... to the end of name if required
        unique_new_name = get_unique_filename(proof_annotation_dir, new_name)

        try:
            # Copy the file to the new location with the unique new name
            shutil.copy(pattern, os.path.join(proof_annotation_dir, unique_new_name))

        except Exception as e:
            files_moved_successfully = False

        if pattern == cropped_boundary_MSA:

            # Generate consensus sequence
            try:
                final_con_object = SequenceManipulator()
                if hmm:  # Generate hmm files when the user want it
                    hmm_dir = os.path.join(os.path.dirname(output_dir), "HMM_files")
                    final_con_object.generate_hmm_from_msa(cropped_boundary_MSA, hmm_dir)

                # Generate final consensus sequence
                # Comparing initial consensus, the final consensus sequence might have different orientation.
                # For this reason, generate consensus again
                final_con = final_con_object.con_generater_no_file(cropped_boundary_MSA, threshold=cons_threshold)
                header = ">" + os.path.splitext(unique_new_name)[0]  # Use the filename without extension
                sequence = str(final_con).upper()

                # Create consensus_object and update length, new consensus name will be created in seq_obj
                # For example self.consensus_name, self.new_length
                consi_obj = seq_obj.create_consi_obj(os.path.splitext(unique_new_name)[0])
                consi_obj.set_new_lenth(len(sequence))

                # Store consensus sequence input consi_obj
                #consi_obj.set_cons_seq(sequence)

                # Store MSA sequence number into consi_obj
                MSA_for_final_cons = AlignIO.read(cropped_boundary_MSA, "fasta")
                MSA_for_final_cons_seq_n = len(MSA_for_final_cons)
                consi_obj.set_cons_MSA_n(MSA_for_final_cons_seq_n)

                # Store terminal repeat to consi_obj
                if found_match == "LTR" or found_match == "TIR":

                    consi_obj.set_new_flank_repeat(found_match)

                # Store blast full length sequence number into consi_obj
                # check_blast_full_n is a function in TE_Aid class. TE Trimmer will use the blast result of TE Aid
                blast_full_length_n = TE_aid_object.check_blast_full_n(consi_obj)
                consi_obj.set_blast_full_n(blast_full_length_n)

                #####################################################################################################
                # Code block: Run RepeatClassifier in RepeatModeler to classify TE_trimmer consensus sequences
                #####################################################################################################
                # Write single consensus sequence to a separate place for RepeatClassifier
                # only classify unknown seq
                # check_unknown returns true if unknown is detected
                # Rename consensus when the final consensus length is much longer or shorter than the query sequence
                if seq_obj.check_unknown() or abs(consi_obj.new_length - seq_obj.old_length) >= 1000:

                    # Define different folder for each sequence
                    sequence_con_file = os.path.join(classification_dir, str(unique_new_name))
                    os.makedirs(sequence_con_file)

                    # Define consensus name path used for classification
                    n_sequence_con_file = os.path.join(sequence_con_file, str(unique_new_name))

                    with open(n_sequence_con_file, "w") as f:

                        # RepeatClassifier input cannot be single sequence, add >Add to enable to run RepeatClassifier
                        f.write(header + "\n" + sequence + "\n" + ">Add" + "\n" + "T" + "\n")

                    TE_type = classify_single(n_sequence_con_file)

                    # Only update new_TE_type when classify_single is successfully
                    if TE_type:
                        # Update TE_type after RepeatClassifier
                        consi_obj.set_new_TE_type(TE_type)

                TE_type = consi_obj.get_TE_type_for_file()
                with open(final_con_file, "a") as f:  # 'a' mode for appending
                    f.write(header + "#" + TE_type + "\n" + sequence + "\n")
                    
            except Exception as e:
                click.echo(f"Error Generate consensus sequence and create consensus_object: {str(e)}")
                click.echo(traceback.format_exc())
                files_moved_successfully = False
        
    return files_moved_successfully
