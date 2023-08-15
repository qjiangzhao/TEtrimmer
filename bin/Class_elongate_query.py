from Class_blast_extension_mafft import SequenceManipulator
from Class_bed_filter import BEDFile
from Class_crop_end_by_gap import CropEndByGap
from Class_crop_end import CropEnd
from Class_group_MSA import MultipleSequenceAlignmentCluster
from Class_select_ditinct_columns import CleanAndSelectColumn
import os
import traceback
import click


class SequenceElongation:
    """
    When query sequence is too short, it won't enable efficient MSA cluster. This class can elongate query sequence.
    when they are too short
    """

    def __init__(self, input_file, genome_file, output_dir, skipped_file, min_seq_num=10):
        self.input_file = input_file
        self.genome_file = genome_file
        self.output_dir = output_dir
        self.skipped_file = skipped_file
        self.min_seq_num = min_seq_num

    def clean_and_cluster_elongate_MSA(self, input_file, output_dir, min_seq_num=10, clean_column_threshold=0.08):

        # Extract columns that contain different alignment patter, use this for group separation.
        # This threshold will be used to replace nucleotides that are less than threshold with a gap character
        # for each column
        pattern_alignment_object = CleanAndSelectColumn(input_file, threshold=clean_column_threshold)

        # Clean_column() function will help MSA cluster and return the absolute path of column cleaned alignment file
        pattern_alignment_object.clean_column(output_dir)

        # Select_divergent_column() function will return a boolean value, true represents need cluster step
        if pattern_alignment_object.select_divergent_column():

            # write_alignment_filtered() function return pattern_alignment absolute path
            pattern_alignment = pattern_alignment_object.write_alignment_filtered(output_dir)

            cluster_MSA_object = MultipleSequenceAlignmentCluster(pattern_alignment, output_dir, min_lines=min_seq_num,
                                                                  max_cluster=1)

            if len(cluster_MSA_object.filtered_cluster_records) == 0:
                return False
            else:

                # Subset input alignment file, return a list contain all pattern alignment clusters
                cluster_intact_alignment_list = cluster_MSA_object.subset_alignment_intact(input_file)

                return cluster_intact_alignment_list
        else:
            return True

    def elongate_query_seq(self):

        input_file_name = os.path.basename(self.input_file)

        try:
            seq_blast = SequenceManipulator()  # this class contain many functions to do blast and filter blast result

            # run blast for each single fasta file and return a bed file absolute path
            bed_out_file_dup = seq_blast.blast(self.input_file, self.genome_file, self.output_dir)

        except Exception as e:
            click.echo(f"Error while running blast for sequence: {input_file_name} elongation. Error: {str(e)}")
            return False

        try:
            # check if blast hit number is equal 0, then skip this sequence
            if seq_blast.blast_hits_count == 0:
                with open(self.skipped_file, "a") as f:
                    f.write(input_file_name + "\telongation_blast_equal_0\n")
                click.echo(f"{input_file_name} is skipped due to elongation blast hit number is 0\n")

                return False

            # check if blast hit number is smaller than 10
            elif seq_blast.blast_hits_count != 0 and seq_blast.blast_hits_count < 10:
                with open(self.skipped_file, "a") as f:
                    f.write(input_file_name + f"\telongation_blast_smaller_than_{self.min_seq_num}\n")
                click.echo(f"{input_file_name} is skipped due to elongation blast hit number is "
                           f"smaller than {self.min_seq_num}\n")
                return False  # when blast hit number is smaller than 10, code will execute next fasta file

        except Exception as e:
            click.echo(f"Error while checking uniqueness for sequence: {input_file_name} elongation. Error: {str(e)}")
            return False

        try:
            # remove duplicated lines
            bed_out_file = seq_blast.check_bed_uniqueness(self.output_dir, bed_out_file_dup)

            # test bed_out_file line number and extract top longest lines
            # return bed_out_filter_file absolute path
            bed_out_filter = BEDFile(bed_out_file)

            # for process_lines() function. threshold represent the maximum number to keep for MSA
            # top_longest_lines_count means the number of sequences with top length
            # for example if threshold = 100, top_longest_lines_count = 50, then 50 sequences will be
            # randomly chose from the rest of sequences
            bed_out_filter_file = bed_out_filter.process_lines(self.output_dir, threshold=100,
                                                               top_longest_lines_count=100)

            # extract fasta from bed_out_filter_file
            # return fasta_out_flank_file absolute path
            # because have to group MSA the first round extend for left and right side are both 0
            fasta_out_flank_file = seq_blast.extract_fasta(bed_out_filter_file, self.genome_file, self.output_dir,
                                                           left_ex=1500, right_ex=1500)

            # Do multiple sequence alignment
            fasta_out_flank_file_MSA = seq_blast.align_sequences(fasta_out_flank_file, self.output_dir)

            # Remove gap block with similarity check
            fasta_out_flank_file_MSA_gap_block_sim = seq_blast.remove_gaps_block_with_similarity_check(
                fasta_out_flank_file_MSA, self.output_dir, gap_threshold=0.8, simi_check_gap_thre=0.4,
                similarity_threshold=0.7, conservation_threshold=0.6, min_nucleotide=5
            )

            # Crop end by gap
            crop_end_by_gap_object = CropEndByGap(fasta_out_flank_file_MSA_gap_block_sim,
                                                  gap_threshold=0.05, window_size=300)
            crop_end_by_gap = crop_end_by_gap_object.write_to_file(self.output_dir)

            # Crop end by divergence
            crop_end_by_divergence_object = CropEnd(crop_end_by_gap, threshold=16, window_size=20)
            crop_end_by_divergence = crop_end_by_divergence_object.write_to_file(self.output_dir)

            # Remove gaps
            crop_end_by_divergence_remove_gap = seq_blast.remove_gaps_with_similarity_check(
                crop_end_by_divergence, self.output_dir, gap_threshold=0.8, simi_check_gap_thre=0.4,
                similarity_threshold=0.7, min_nucleotide=5
            )

            # Return False when cluster number is 0, otherwise, return True when divergent column number is smaller
            # than 100 otherwise it will return the subset intact MSA file absolute path
            cluster_MSA_result = self.clean_and_cluster_elongate_MSA(crop_end_by_divergence_remove_gap, self.output_dir,
                                                                     min_seq_num=self.min_seq_num,
                                                                     clean_column_threshold=0.08)
        except Exception as e:
            click.echo(
                f"Error during processing lines, extracting fasta, or clustering MSA for"
                f" sequence: {input_file_name} elongation Error: {str(e)}")
            traceback.print_exc()
            return False

        try:
            if cluster_MSA_result is False:  # When it is False means all clusters' sequence number is less than 10
                with open(self.skipped_file, "a") as f:
                    f.write(input_file_name + f"\telongation_cluster_blast_smaller_than_{self.min_seq_num}\n")
                click.echo(f"{input_file_name} is skipped due to sequence number in elongation cluster is "
                           f"smaller than {self.min_seq_num}\n")
                return False

            # When it is True means not necessary to do clustering. Use the input MSA to generate consensus sequence
            elif cluster_MSA_result is True:
                output_con = seq_blast.con_generater(crop_end_by_divergence_remove_gap, self.output_dir, threshold=0.8,
                                                     ambiguous="N")
            else:
                # Generate consensus sequence
                output_con = seq_blast.con_generater(str(cluster_MSA_result[0]), self.output_dir, threshold=0.8,
                                                     ambiguous="N")

            # Change consensus sequence name, otherwise it will be too long for the further analysis
            new_output_con = os.path.join(self.output_dir, input_file_name)
            os.rename(output_con, new_output_con)

            return new_output_con

        except Exception as e:
            click.echo(f"Error during consensus generation for sequence: {input_file_name} elongation Error: {str(e)}")
            traceback.print_exc()
