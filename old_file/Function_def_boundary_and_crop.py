from Class_Define_boundary import DefineBoundary
from Class_crop_end import CropEnd
from Class_crop_end_by_gap import CropEndByGap
from Class_blast_extension_mafft import SequenceManipulator
from Class_select_ditinct_columns import CleanAndSelectColumn
from Class_MSA_plot import MSAPainter
from Class_check_TGT_ACA import StartEndChecker
from Class_TE_aid import TEAid
from Class_pfam_scan_and_plot import PlotPfam
from Bio import AlignIO
import os
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

def find_boundary_and_crop(bed_file, genome_file, output_dir, pfam_dir, cons_threshold=0.7,
                           left_ex=1000, right_ex=1000, crop_end_thr=16, crop_end_win=20):
    # Check if this is a LINE element, if so decrease the cons_threshold number,
    # because LINE have higher divergence at 5' end
    if"LINE" in bed_file:
        cons_threshold = cons_threshold - 0.2
    if_left_ex = True
    if_right_ex = True
    bed_fasta_mafft = None
    bed_boundary = None

    #####################################################################################################
    # Code block: Define left and right sides extension number
    #####################################################################################################

    # The maximum extension is 7000
    while (if_left_ex or if_right_ex) and (left_ex <= 7000 and right_ex <= 7000):
        bed = SequenceManipulator()
        bed_fasta = bed.extract_fasta(bed_file, genome_file, output_dir, left_ex, right_ex)
        # align_sequences() will return extended MSA absolute file
        bed_fasta_mafft_with_gap = bed.align_sequences(bed_fasta, output_dir)
        # Remove nucleotide whose proportion is smaller than threshold
        bed_fasta_mafft_with_gap_column_clean_object = CleanAndSelectColumn(bed_fasta_mafft_with_gap, threshold=0.08)
        bed_fasta_mafft_with_gap_column_clean = bed_fasta_mafft_with_gap_column_clean_object.clean_column(output_dir)

        bed_fasta_mafft = bed.remove_gaps_with_similarity_check(bed_fasta_mafft_with_gap_column_clean, output_dir,
                                                                   gap_threshold=0.6, simi_check_gap_thre=0.4,
                                                                   similarity_threshold=0.7,
                                                                   min_nucleotide=5
                                                                   )

        bed_fasta_mafft_object = CropEndByGap(bed_fasta_mafft, gap_threshold=0.05, window_size=300)

        bed_fasta_mafft_cop_end_gap = bed_fasta_mafft_object.write_to_file(output_dir)

        # Threshold here means the threshold to generate consensus sequence
        bed_boundary = DefineBoundary(bed_fasta_mafft_cop_end_gap, threshold=cons_threshold,
                                      check_window=200, max_X=40, if_con_generater=False)

        if bed_boundary.left_ext:  # boundary.left_ext will become true when more extension is required
            left_ex += 1000
        else:
            if_left_ex = bed_boundary.left_ext

        if bed_boundary.right_ext:
            right_ex += 1000
        else:
            if_right_ex = bed_boundary.right_ext

    # crop bed_fasta_mafft according to the start and end position
    # return absolute path for boundary cropped sequence

    #####################################################################################################
    # Code block: Remove gaps and define the start and end position
    #####################################################################################################

    bed_fasta_mafft_gap_block_sim = bed.remove_gaps_block_with_similarity_check(
        bed_fasta_mafft_with_gap_column_clean, output_dir, gap_threshold=0.8, similarity_threshold=0.7,
        conservation_threshold=0.5
    )

    bed_boundary = DefineBoundary(bed_fasta_mafft_gap_block_sim, threshold=cons_threshold, check_window=200, max_X=40)
    bed_fasta_mafft_boundary_crop = bed_boundary.crop_MSA(output_dir, crop_extension=300)

    # Because for LINE element bed_fasta_mafft_boundary_crop will be changed. Copy it to the other variable
    bed_fasta_mafft_boundary_crop_for_select = bed_fasta_mafft_boundary_crop

    if "LINE" in bed_file:
        # For the high divergence region, more gaps can be found. According to this feature, remove high divergence region
        # this function is very useful for dealing with LINE elements
        cropped_MSA_by_gap = CropEndByGap(bed_fasta_mafft_boundary_crop, gap_threshold=0.05, window_size=300)
        bed_fasta_mafft_boundary_crop = cropped_MSA_by_gap.write_to_file(output_dir)

    cropped_alignment_output_file_no_gap, column_mapping = crop_end_and_clean_column(bed_fasta_mafft_boundary_crop,
                                                                                       output_dir,
                                                                                       crop_end_threshold=crop_end_thr,
                                                                                       window_size=crop_end_win,
                                                                                       gap_threshold=0.8)

    # Crop end can't define the final boundary, use DefineBoundary again to define start position
    #cropped_boundary = DefineBoundary(cropped_alignment_output_file_no_gap, threshold=0.8, check_window=6, max_X=1)
    #cropped_boundary_MSA = cropped_boundary.crop_MSA(output_dir, crop_extension=0)

    cropped_boundary = DefineBoundary(cropped_alignment_output_file_no_gap, threshold=0.8, check_window=3, max_X=0)
    cropped_boundary_MSA = cropped_boundary.crop_MSA(output_dir, crop_extension=0)

    #####################################################################################################
    # Code block: For LTR element, check if the cropped MSA starts with TGT and ends with ACA
    #####################################################################################################

    # Create StartEndChecker object
    check_start_end_object = StartEndChecker(cropped_alignment_output_file_no_gap,
                                             start=cropped_boundary.start_post,
                                             end=cropped_boundary.end_post,
                                             start_pattern="TGT", end_pattern="ACA", threshold=0.7)

    if check_start_end_object.check_LTR():  # Check if file name contain "LTR"
        check_start_end_object.con_generater()  # Generate consensus sequences

        # Three variables will be returned, if "check_result" is false, that means "TGT" and "ACA" are not found
        check_result, check_start, check_end = check_start_end_object.check_and_update(other_end_pattern=["aga", "ata"])

        if check_result:
            # When the input start is equal the output start position, that means the original sequence start with "TGT"
            # same like the end position. If one of them are different, convey the new start or end position to
            # cropped_boundary object and generate the new MSA file
            if cropped_boundary.start_post != check_start or cropped_boundary.end_post != check_end:
                cropped_boundary.start_post = check_start
                cropped_boundary.end_post = check_end
                # Generate the new MSA file based on new start and end positions
                cropped_boundary_MSA = cropped_boundary.crop_MSA(output_dir, crop_extension=0)

    #####################################################################################################
    # Code block: Generate MSA for manual curation when the start and end point are wrong
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
    # to a alignment object.
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

    orf_doamin_plot_object = PlotPfam(TE_aid_consensus, output_dir, pfam_database_dir=pfam_dir)
    orf_domain_plot = orf_doamin_plot_object.orf_domain_plot()

    #####################################################################################################
    # Code block: Merge plot files
    #####################################################################################################

    # Create an instance of PdfFileMerger class
    merger = PdfMerger()

    # List with your two PDF file paths
    pdf_files = [MSA_plot, orf_domain_plot, TE_aid_plot]

    # Iterate over the list of the file paths
    for pdf_file in pdf_files:
        # Append PDF files
        merger.append(pdf_file)

    # Write out the merged PDF file
    merged_pdf_path = os.path.join(output_dir, f"{os.path.basename(cropped_boundary_MSA)}_me_plot.pdf")
    merger.write(merged_pdf_path)
    merger.close()

    return
