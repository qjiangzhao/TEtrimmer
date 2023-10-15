from Class_select_ditinct_columns import CleanAndSelectColumn
from Class_group_MSA import MultipleSequenceAlignmentCluster
import Function_blast_extension_mafft


def clean_and_cluster_MSA(input_file, bed_file, output_dir, div_column_thr=0.8, clean_column_threshold=0.08,
                          min_length_num=10, cluster_num=2, cluster_col_thr=500, muscle_ite_times=4, fast_mode=False):
    """
    This function will cluster multiple sequence alignment file
    :param input_file: str, The direct fasta file derived from bed file
    :param bed_file: The bed file used to generate pattern alignment
    :param output_dir: Output directory
    :param gap_threshold: num (0-1) default 0.8, columns with gap percentage higher than "gap_threshold" will be removed
    :param clean_column_threshold: num (0-1) default 0.08, nucleotide percentage (gap not count)
    lower than threshold will be converted to "-"
    :param min_length_num: num default 10, the minimum line number for each cluster
    :param cluster_num: num default 2, the maximum cluster number for each MSA
    :return: A list of subset pattern alignment and bed files
    """

    # Align_sequences will return the absolute file path of alignment file
    if fast_mode:
        muscle_ite_times = 2
    try:
        fasta_out_flank_mafft_file = Function_blast_extension_mafft.muscle_align(input_file, output_dir,
                                                                                 ite_times=muscle_ite_times)
    except Exception as e:
        fasta_out_flank_mafft_file = False
        pass

    # When muscle goes wrong, use mafft
    if not fasta_out_flank_mafft_file:
        fasta_out_flank_mafft_file = Function_blast_extension_mafft.align_sequences(input_file, output_dir)

    # Remove gaps. Return absolute path for gap removed alignment file
    fasta_out_flank_mafft_file_gap_filter = Function_blast_extension_mafft.remove_gaps_with_similarity_check(
        fasta_out_flank_mafft_file, output_dir, gap_threshold=0.8, simi_check_gap_thre=0.4,
        similarity_threshold=0.85, min_nucleotide=5)

    # Extract columns that contain different alignment patter, use this for group separation.
    # This threshold will be used to replace nucleotides that are less than threshold with a gap character for each column
    pattern_alignment = CleanAndSelectColumn(fasta_out_flank_mafft_file_gap_filter, threshold=clean_column_threshold)

    # Clean_column() function will help MSA cluster and return the absolute path of column cleaned alignment file
    pattern_alignment.clean_column(output_dir)

    # Select_divergent_column() function will return a boolean value, true represents need cluster step
    if pattern_alignment.select_divergent_column(cluster_col_thr=cluster_col_thr, dis_col_threshold=div_column_thr):

        # write_alignment_filtered() function return pattern_alignment absolute path
        pattern_alignment = pattern_alignment.write_alignment_filtered(output_dir)
        """
        MultipleSequenceAlignmentCluster() function requires pattern alignment. It also require bed file for subsetting 
        bed file.
        "min_lines" will define the minimum line numbers for each cluster
        "max_cluster" will define the maximum cluster numbers that will be returned
        this function will return a list contain all cluster file absolute path
        """
        cluster_MSA_object = MultipleSequenceAlignmentCluster(pattern_alignment, output_dir,
                                                              min_cluster_size=min_length_num,
                                                              max_cluster=cluster_num)

        # Test if silhouette_scores is high enough to perform cluster
        if cluster_MSA_object.if_cluster:

            """
            Test if cluster number is 0, if so skip this sequence
            If cluster is 0, that means no cluster has line numbers greater than "min_lines". In this case, 
            it will be hard to still use multiple sequence alignment method to define consensus sequence.
            that isn't to say this won't be a TE, but with less copy numbers. Low copy TE will also be checked later
            """
            if len(cluster_MSA_object.filtered_cluster_records) == 0:

                return False

            else:
                # Subset bed file and return a list contain all clustered bed absolute files
                cluster_bed_files_list = cluster_MSA_object.subset_bed_file(bed_file)
                return cluster_bed_files_list
        else:
            return True

    else:
        return True





