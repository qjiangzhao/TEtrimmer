from Class_Define_boundary import DefineBoundary
from Class_crop_end import CropEnd
from Class_crop_end_by_gap import CropEndByGap
from Class_blast_extension_mafft import SequenceManipulator
from Class_select_ditinct_columns import CleanAndSelectColumn
from Class_MSA_plot import MSAPainter
from Class_check_start_end import StartEndChecker
from Class_TE_aid import TEAid
from Class_pfam_scan_and_plot import PlotPfam
from Bio import AlignIO
import os
import re
import shutil
from PyPDF2 import PdfMerger


def crop_end_and_clean_column(input_file, output_dir, crop_end_threshold=16, window_size=20, gap_threshold=0.8):

    # Window_size means the checked nucleotide number each time
    # Threshold means the sum of nucleotide proportion in window_size must greater than that
    cropped_MSA = CropEnd(input_file, threshold=crop_end_threshold, window_size=window_size)

    # write_to_file() function will return absolute cropped alignment file
    cropped_MSA_output_file = cropped_MSA.write_to_file(output_dir)

    # Remove gaps again after crop end step
    cropped_alignment_output_file_with_gap = SequenceManipulator()

    # "column_mapping" is a dictionary the key is gap removed MSA index, the value is the corresponded original MSA index
    cropped_alignment_output_file_no_gap, column_mapping = cropped_alignment_output_file_with_gap.remove_gaps(
        cropped_MSA_output_file, output_dir, threshold=gap_threshold, min_nucleotide=5)

    return cropped_alignment_output_file_no_gap, column_mapping


def find_boundary_and_crop(bed_file, genome_file, output_dir, pfam_dir, seq_name, cons_threshold=0.8,
                           ext_threshold=0.7, ex_step_size=1000, max_extension=7000,
                           gap_threshold=0.4, gap_nul_thr=0.7, crop_end_thr=16, crop_end_win=20,
                           crop_end_gap_thr=0.05, crop_end_gap_win=300, start_patterns=None, end_patterns=None,
                           mini_orf=200):
    """
    :param bed_file: str, bed file directory
    :param genome_file: str, genome directory
    :param output_dir: str, output file directory
    :param pfam_dir: str, pfam database directory
    :param con_threshold: num 0 to 1 default: 0.8, threshold used for final consensus sequence generation
    :param ext_threshold: num 0 to 1 default: 0.7, threshold used for define the extension extent. The higher number
    represent more strigent extension threshold (short extension in total), the smaller number means it become easier
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
    if"LINE" in bed_file:
        ext_threshold = ext_threshold - 0.2
    if_left_ex = True
    if_right_ex = True

    #####################################################################################################
    # Code block: Define left and right sides extension number
    #####################################################################################################

    left_ex = ex_step_size
    right_ex = ex_step_size
    while (if_left_ex or if_right_ex) and (left_ex <= max_extension and right_ex <= max_extension):
        bed = SequenceManipulator()
        # bedtools will makesure the extension won't excess the maximum length of that chromosome
        bed_fasta = bed.extract_fasta(bed_file, genome_file, output_dir, left_ex, right_ex)

        # align_sequences() will return extended MSA absolute file
        bed_fasta_mafft_with_gap = bed.align_sequences(bed_fasta, output_dir)

        # Remove nucleotide whose proportion is smaller than threshold
        bed_fasta_mafft_with_gap_column_clean_object = CleanAndSelectColumn(bed_fasta_mafft_with_gap, threshold=0.08)
        bed_fasta_mafft_with_gap_column_clean = bed_fasta_mafft_with_gap_column_clean_object.clean_column(output_dir)

        # Remove gaps with similarity check
        bed_fasta_mafft = bed.remove_gaps_with_similarity_check(bed_fasta_mafft_with_gap_column_clean, output_dir,
                                                                   gap_threshold=0.6, simi_check_gap_thre=0.4,
                                                                   similarity_threshold=0.7,
                                                                   min_nucleotide=5
                                                                   )

        bed_fasta_mafft_object = CropEndByGap(bed_fasta_mafft, gap_threshold=0.05, window_size=300)

        bed_fasta_mafft_cop_end_gap = bed_fasta_mafft_object.write_to_file(output_dir)

        # Threshold here means the threshold to generate consensus sequence
        bed_boundary = DefineBoundary(bed_fasta_mafft_cop_end_gap, threshold=ext_threshold,
                                      check_window=200, max_X=0.2, if_con_generater=False)

        if bed_boundary.left_ext:  # boundary.left_ext will become true when more extension is required
            left_ex += ex_step_size
        else:
            if_left_ex = bed_boundary.left_ext

        if bed_boundary.right_ext:
            right_ex += ex_step_size
        else:
            if_right_ex = bed_boundary.right_ext

    #####################################################################################################
    # Code block: Remove gaps and define the start and end position
    #####################################################################################################

    # Remove gap block with similarity check
    bed_fasta_mafft_gap_block_sim = bed.remove_gaps_block_with_similarity_check(
        bed_fasta_mafft_with_gap_column_clean, output_dir, gap_threshold=0.8, simi_check_gap_thre=gap_threshold,
        similarity_threshold=gap_nul_thr, conservation_threshold=0.5
    )

    bed_boundary = DefineBoundary(bed_fasta_mafft_gap_block_sim, threshold=ext_threshold, check_window=200, max_X=0.2)
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
    # Code block: Generate MSA for manual curation
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
        bed_fasta_mafft_boundary_crop_for_select, output_dir, column_mapping[cropped_boundary.start_post], "left", window_size=50
    )

    # Get 50 columns right the end point
    cropped_boundary_plot_right = cropped_boundary_plot_object.select_window_columns(
        bed_fasta_mafft_boundary_crop_for_select, output_dir, column_mapping[cropped_boundary.end_post], "right", window_size=50
    )

    # Concatenate MSA
    # The concat_start and concat_end points correspond with cropped_boundary.start_post and cropped_boundary.end_post
    cropped_boundary_plot_concatenate, concat_start, concat_end = cropped_boundary_plot_object.concatenate_alignments(
        cropped_boundary_plot_left, cropped_boundary_plot_select_start_end_and_joint, cropped_boundary_plot_right,
        input_file_name=cropped_boundary_MSA, output_dir=output_dir
    )

    # Plot MSA, which can easily verify if the start and end crop points are right
    MSA_plot_object = MSAPainter(cropped_boundary_plot_concatenate, output_dir, sequence_len)
    MSA_plot = MSA_plot_object.process(concat_start, concat_end)

    #####################################################################################################
    # Code block: Check by TE-Aid and combine MSA plot with TE-Aid plot
    #####################################################################################################

    # Generate consensus sequence for TE_aid
    TE_aid_consensus_object = SequenceManipulator()
    # The low threshold will make sure most of columns will generate a consensus nucleotide that is import for orf prediction
    TE_aid_consensus = TE_aid_consensus_object.con_generater(cropped_boundary_MSA, output_dir, threshold=0.5)

    # so.path.abspath(__file__) will return the current executable python file
    # TE-Aid package is stored at the same directory as this function file.
    TE_aid_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "TE-Aid-master")
    # Run TE_aid
    TE_aid_object = TEAid(TE_aid_consensus, output_dir, genome_file, TE_aid_dir=TE_aid_path)
    TE_aid_plot = TE_aid_object.run()

    #####################################################################################################
    # Code block: Predict orf, scan with Pfam, and plot orf and domain features
    #####################################################################################################

    # Code block: Predict orf, scan with Pfam, and plot orf and domain features
    orf_domain_plot_object = PlotPfam(TE_aid_consensus, output_dir, pfam_database_dir=pfam_dir, mini_orf=mini_orf)
    orf_domain_plot = orf_domain_plot_object.orf_domain_plot() if orf_domain_plot_object.orf_result else None

    #####################################################################################################
    # Code block: Merge plot files
    #####################################################################################################

    # Code block: Merge plot files
    merger = PdfMerger()

    # List with your PDF file paths. Ensure they are not None before adding them to the list.
    pdf_files = [pdf for pdf in [MSA_plot, orf_domain_plot, TE_aid_plot] if pdf is not None]

    valid_pdf_count = 0  # Counter to keep track of valid PDFs added

    # Iterate over the list of the file paths
    for pdf_file in pdf_files:
        # Check if the file exists before appending
        if os.path.exists(pdf_file):
            # Append PDF files
            merger.append(pdf_file)
            valid_pdf_count += 1
        else:
            print(f"Warning: {pdf_file} does not exist and will not be merged.")

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

    # Construct the path for the new folder
    proof_annotation_dir = os.path.join(parent_output_dir, "TE_Trimmer_for_proof_annotation")

    # Create the directory if it doesn't exist
    if not os.path.exists(proof_annotation_dir):
        os.makedirs(proof_annotation_dir)

    # Because for each query file, several subgroups of multiple sequence alignment will be generated,
    # This function will name them differently
    def get_unique_filename(base_path, original_name):
        counter = 0
        name, ext = os.path.splitext(original_name)

        # If the extension is .anno_fasta, change it to .anno.fasta
        if ext == ".anno_fasta":
            ext = ".anno.fasta"

        # Initialize with the original name
        new_name = f"{name}{ext}"

        # While the constructed name already exists, increment the counter and reconstruct the name
        while os.path.exists(os.path.join(base_path, new_name)):
            counter += 1
            new_name = f"{name}_{counter:02}{ext}"  # Adjusted to use zero-padding

        return new_name

    # Eliminate .fasta from seq_name, this can make sure the right file order when do proof annotation
    seq_name_name, seq_name_ext = os.path.splitext(seq_name)

    file_copy_pattern = [
        (rf"^{seq_name}.*me_plot.pdf$", f"{seq_name_name}.pdf"),
        (rf"^{seq_name}.*gap_rm.fa_bou_crop.fa$", f"{seq_name_name}.fasta"),  # Changed to append .fasta
        (rf"^{seq_name}.*_proof_anno_me.fa$", f"{seq_name_name}.anno_fasta")
    ]

    final_msa_pattern = rf"^{seq_name}.*gap_rm.fa_bou_crop.fa$"

    files = os.listdir(output_dir)
    for file in files:
        for pattern, new_name in file_copy_pattern:
            if re.match(pattern, file):
                # Decide on the final unique name before copying
                unique_new_name = get_unique_filename(proof_annotation_dir, new_name)

                # Copy the file to the new location with the unique new name
                shutil.copy(os.path.join(output_dir, file), os.path.join(proof_annotation_dir, unique_new_name))
                os.remove(os.path.join(output_dir, file))

                # Generate consensus sequence
                if re.match(final_msa_pattern, file):
                    final_con_object = SequenceManipulator()
                    final_con = final_con_object.con_generater_no_file(
                        os.path.join(proof_annotation_dir, unique_new_name), threshold=cons_threshold)
                    header = ">" + os.path.splitext(unique_new_name)[0]  # Use the filename without extension
                    sequence = str(final_con).upper()
                    with open(final_con_file, "a") as f:  # 'a' mode for appending
                        f.write(header + "\n" + sequence + "\n")

    return cropped_boundary
