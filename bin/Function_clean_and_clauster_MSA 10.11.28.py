from Class_select_ditinct_columns import CleanAndSelectColumn
from Class_group_MSA import MultipleSequenceAlignmentCluster
from Class_blast_extension_mafft import SequenceManipulator

def clean_and_cluster_MSA(input_file, bed_file, output_dir, gap_threshold=0.8, clean_column_threshold=0.08,
                          min_length_num=10, cluster_num=2):
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

    seq_blast = SequenceManipulator()
    # Align_sequences will return the absolute file path of alignment file
    fasta_out_flank_mafft_file = seq_blast.align_sequences(input_file, output_dir)

    # Remove gaps. Return absolute path for gap removed alignment file
    fasta_out_flank_mafft_file_gap_filter, column_mapping = seq_blast.remove_gaps(
        fasta_out_flank_mafft_file, output_dir, threshold=gap_threshold)

    # Extract columns that contain different alignment patter, use this for group separation.
    # This threshold will be used to rplace nucleotides that are less than threshold with a gap character for each column
    pattern_alignment = CleanAndSelectColumn(fasta_out_flank_mafft_file_gap_filter, threshold=clean_column_threshold)

    # Clean_column() function will help MSA cluster and return the absolute path of column cleaned alignment file
    pattern_alignment.clean_column(output_dir)

    # Select_divergent_column() function will return a boolean value, true represents need cluster step
    if pattern_alignment.select_divergent_column():
        # write_alignment_filtered() function return pattern_alignment absolute path
        pattern_alignment = pattern_alignment.write_alignment_filtered(output_dir)
        """
        MultipleSequenceAlignmentCluster() function requires pattern alignment. It also require bed file for subseting 
        bed file.
        "min_lines" will define the minimum line numbers for each cluster
        "max_cluster" will define the maximum cluster numbers that will be returned
        this function will return a list contain all cluster file absolute path
        """
        cluster_MSA = MultipleSequenceAlignmentCluster(pattern_alignment, output_dir, min_lines=min_length_num,
                                                       max_cluster=cluster_num)

        """
        Test if cluster number is 0, if so skip this sequence
        If cluster is 0, that means no cluster has line numbers greater than "min_lines". In this case, it will be hard
        to still use multiple sequence alignment method to define consensus sequence.
        that isn't to say this won't be a TE, but with less copy numbers.
        """
        if len(cluster_MSA.filtered_cluster_records) == 0:

            return False
        else:
            # Subset pattern alignment file, return a list contain all pattern alignment clusters
            cluster_pattern_alignment_list = cluster_MSA.subset_alignment_dis()
            # Subset bed file and return a list contain all clustered bed absolute files
            cluster_bed_files_list = cluster_MSA.subset_bed_file(bed_file)

            return cluster_pattern_alignment_list, cluster_bed_files_list
    else:
        return True





