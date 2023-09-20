import os
import concurrent.futures
import shutil
import traceback
from Bio import SeqIO
from datetime import datetime, timedelta
from Class_separate_fasta import FastaSequenceSeparator
from Class_blast_extension_mafft import SequenceManipulator
from Class_bed_filter import BEDFile
from Function_def_boundary_and_crop import find_boundary_and_crop
from Function_clean_and_clauster_MSA import clean_and_cluster_MSA
from Class_elongate_query import SequenceElongation
from Class_orf_domain_prediction import prepare_pfam_database
import click
import json


# Define function to check progress file, which will be used for continue analysis
def check_progress_file(progress_file_path):
    try:
        # If the progress file exists, open the file and read the list of completed sequences
        with open(progress_file_path, 'r') as f:
            local_completed_sequences = f.readlines()

        # Strip newline characters and remove duplicates
        local_completed_sequences = list(
            set([sequence.strip().split('\t')[0] for sequence in local_completed_sequences]))

    except Exception as e:
        raise Exception(
            f"An error occurred while reading the progress file {progress_file_path}: {e}")
    return local_completed_sequences


#####################################################################################################
# Code block: Define analyze_sequence function
#####################################################################################################


def analyze_sequence_helper(params):
    return analyze_sequence(*params)


def analyze_sequence(seq_name, single_file_dir, genome_file, MSA_dir, min_blast_len, min_seq_num, max_msa_lines,
                     top_mas_lines, max_cluster_num, min_el, min_el_dna, min_el_sine, cons_thr, ext_thr, ex_step,
                     max_extension, gap_thr, gap_nul_thr, crop_end_thr, crop_end_win, crop_end_gap_thr,
                     crop_end_gap_win, start_patterns, end_patterns, output_dir, pfam_dir, mini_orf, elongation_dir,
                     single_fasta_n, hmm, check_extension_win, keep_intermediate):
    #####################################################################################################
    # Code block: Elongate query sequence when it is too short
    #####################################################################################################

    try:
        # Get query fasta file path
        seq_file = os.path.join(single_file_dir, seq_name)

        # Define the skipped_file, all single sequence that is skipped by TE Trimmer will be stored inside this file
        skipped_file = os.path.join(output_dir, "Skipped_sequence_name.txt")

        # Define the progress file
        progress_file = os.path.join(output_dir, "Finished_sequence_name.txt")

        # If sequence length is less than min_length_elongation. Elongation process will be executed.
        # Because short query sequence won't enable efficient MSA cluster
        seq_record = SeqIO.read(seq_file, "fasta")
        seq_file_len = len(seq_record.seq)
        min_length_elongation = min_el

        # TODO add more elements for different length
        elongation_ext_n = 1500
        # Due to DNA element will be much shorter than LTR and LINE elements, set min_length_elongation to 500.
        if "DNA" in seq_name:
            min_length_elongation = min_el_dna
            elongation_ext_n = 800
            ex_step = 500
            max_extension = 7000
            min_blast_len = 150
            crop_end_gap_win = 100
            check_extension_win = 50

        # The average length of SINE element is around 500 bp, set smaller min_length_elongation number.
        if "SINE" in seq_name:
            min_length_elongation = min_el_sine
            elongation_ext_n = 200
            ex_step = 200
            max_extension = 1400
            min_blast_len = 80
            crop_end_gap_win = 50
            check_extension_win = 50

        if "Helitron" in seq_name:
            min_length_elongation = 500
            elongation_ext_n = 500
            ex_step = 500
            max_extension = 7000
            min_blast_len = 150
            crop_end_gap_win = 100
            check_extension_win = 50

        if "MITE" in seq_name:
            min_length_elongation = 50
            elongation_ext_n = 200
            ex_step = 100
            max_extension = 500
            min_blast_len = 50
            crop_end_gap_win = 40
            check_extension_win = 50

        if seq_file_len < min_length_elongation:

            seq_elongation_object = SequenceElongation(seq_file, genome_file, elongation_dir, skipped_file,
                                                       min_seq_num=min_seq_num,
                                                       crop_end_gap_thr=crop_end_gap_thr,
                                                       crop_end_gap_win=crop_end_gap_win,
                                                       ext_n=elongation_ext_n)
            seq_elongation = seq_elongation_object.elongate_query_seq()

            if seq_elongation:
                seq_file = seq_elongation

            # When seq_elongation is false, skip this sequence.
            elif not seq_elongation:
                return

        # run blast for each single fasta file and return a bed file absolute path
        seq_blast = SequenceManipulator()  # this class contain many functions to do blast and filter blast result
        bed_out_file_dup = seq_blast.blast(seq_file, genome_file, MSA_dir, min_length=min_blast_len)

    except Exception as e:
        click.echo(f"Error while running blast for sequence: {seq_name}. Error: {str(e)}")
        return

    try:
        # check if blast hit number is equal 0, then skip this sequence
        if seq_blast.blast_hits_count == 0:
            with open(skipped_file, "a") as f:
                f.write(seq_name + "\tno_blast\tblast_equal_0\n")
            click.echo(f"\n{seq_name} is skipped due to blast hit number is 0\n")

            if not keep_intermediate:
                remove_file_object = SequenceManipulator()
                remove_file_object.remove_files_with_start_pattern(MSA_dir, seq_name)
                remove_file_object.remove_files_with_start_pattern(elongation_dir, seq_name)
            return

        # check if blast hit number is smaller than "min_seq_num", not include "min_seq_num"
        elif seq_blast.blast_hits_count != 0 and seq_blast.blast_hits_count < min_seq_num:
            with open(skipped_file, "a") as f:
                f.write(seq_name + f"\tblast_less\tblast_smaller_than_{min_seq_num}\n")
            click.echo(f"\n{seq_name} is skipped due to blast hit number is smaller than {min_seq_num}\n")

            if not keep_intermediate:
                remove_file_object = SequenceManipulator()
                remove_file_object.remove_files_with_start_pattern(MSA_dir, seq_name)
                remove_file_object.remove_files_with_start_pattern(elongation_dir, seq_name)
            return  # when blast hit number is smaller than 10, code will execute next fasta file
            # with less hits numbers can also refer to young TE, but less blast hit number will hamper to generate
            # consensus sequence by multiple sequence alignment method.
            # TODO check if blast hit lengths occupied 80% of query if so, check if query have LTR and TE domain
            # TODO if it fits critium, keep it to final TE library but label with single copy TE
            # TODO write this fasta name to a file and do 4000 extension

    except Exception as e:
        click.echo(f"Error while checking uniqueness for sequence: {seq_name}. Error: {str(e)}")
        return

    try:
        # remove duplicated lines
        bed_out_file = seq_blast.check_bed_uniqueness(MSA_dir, bed_out_file_dup)

        # test bed_out_file line number and extract top longest lines
        # return bed_out_filter_file absolute path
        bed_out_filter = BEDFile(bed_out_file)
        # for process_lines() function. threshold represent the maximum number to keep for MSA
        # top_longest_lines_count means the number of sequences with top length
        # for example if threshold = 100, top_longest_lines_count = 50, then 50 sequences will be
        # randomly chose from the rest of sequences
        # top_mas_lines has to be equal or smaller than max_mas_lines
        bed_out_filter_file = bed_out_filter.process_lines(MSA_dir, threshold=max_msa_lines,
                                                           top_longest_lines_count=top_mas_lines)

        # extract fast from bed_out_filter_file
        # return fasta_out_flank_file absolute path
        # because have to group MSA the first round extend for left and right side are both 0
        fasta_out_flank_file, bed_out_flank_file = seq_blast.extract_fasta(
            bed_out_filter_file, genome_file, MSA_dir, left_ex=0, right_ex=0)

        # Return False when cluster number is 0. Return True when divergent column number is smaller than 100
        # Otherwise it will return the subset bed and alignment file
        cluster_MSA_result = clean_and_cluster_MSA(fasta_out_flank_file, bed_out_filter_file, MSA_dir,
                                                   clean_column_threshold=0.08,
                                                   min_length_num=min_seq_num, cluster_num=max_cluster_num,
                                                   cluster_col_thr=100)
    except Exception as e:
        click.echo(
            f"Error during processing lines, extracting fasta, or clustering MSA for sequence: {seq_name}. Error: {str(e)}")
        traceback.print_exc()
        return

    try:
        # cluster false means no cluster, TE Trimmer will skip this sequence.
        if cluster_MSA_result is False:
            with open(skipped_file, "a") as f:
                f.write(seq_name + f"\tblast_less\tcluster_blast_smaller_than_{min_seq_num}\n")
            click.echo(f"\n{seq_name} is skipped due to sequence number in each cluster is smaller than {min_seq_num}\n")

            if not keep_intermediate:
                remove_file_object = SequenceManipulator()
                remove_file_object.remove_files_with_start_pattern(MSA_dir, seq_name)
                remove_file_object.remove_files_with_start_pattern(elongation_dir, seq_name)
            return

        # cluster True means not necessary to cluster MSA
        elif cluster_MSA_result is True:
            find_boundary_result = find_boundary_and_crop(
                    bed_out_filter_file, genome_file, MSA_dir, pfam_dir, seq_name, hmm,
                    cons_threshold=cons_thr, ext_threshold=ext_thr,
                    ex_step_size=ex_step, max_extension=max_extension,
                    gap_threshold=gap_thr, gap_nul_thr=gap_nul_thr,
                    crop_end_thr=crop_end_thr, crop_end_win=crop_end_win,
                    crop_end_gap_thr=crop_end_gap_thr, crop_end_gap_win=crop_end_gap_win,
                    start_patterns=start_patterns, end_patterns=end_patterns, mini_orf=mini_orf,
                    define_boundary_win=check_extension_win
            )

            if find_boundary_result == "Short_sequence":
                with open(skipped_file, "a") as f:
                    f.write(seq_name + f"\tshort\telement is too short for TE Trimmer\n")
                click.echo(f"\n{seq_name} is skipped due to too short length\n")

                if not keep_intermediate:
                    remove_file_object = SequenceManipulator()
                    remove_file_object.remove_files_with_start_pattern(MSA_dir, seq_name)
                    remove_file_object.remove_files_with_start_pattern(elongation_dir, seq_name)
                return False
            elif not find_boundary_result:  # This means the errors happen in the function
                return

        else:
            cluster_bed_files_list = cluster_MSA_result

            # cluster_pattern_alignment_list have the same index with cluster_bed_files_list
            all_inner_skipped = True
            for i in range(len(cluster_bed_files_list)):

                # Based on the bed file list, extract fasta file
                inner_fasta_out_flank_file, inner_bed_out_flank_file = seq_blast.extract_fasta(
                    cluster_bed_files_list[i], genome_file, MSA_dir, left_ex=0, right_ex=0)

                inner_cluster_MSA_result = clean_and_cluster_MSA(inner_fasta_out_flank_file,
                                                                 cluster_bed_files_list[i], MSA_dir,
                                                                 clean_column_threshold=0.08,
                                                                 min_length_num=min_seq_num,
                                                                 cluster_num=max_cluster_num,
                                                                 cluster_col_thr=100)
                # inner_cluster_MSA_result is false means this cluster sequence number is too less
                if inner_cluster_MSA_result is False:
                    continue
                elif inner_cluster_MSA_result is True:  # Means don't need to cluster

                    inner_find_boundary_result = find_boundary_and_crop(
                            cluster_bed_files_list[i], genome_file, MSA_dir, pfam_dir, seq_name,
                            hmm, cons_threshold=cons_thr, ext_threshold=ext_thr,
                            ex_step_size=ex_step, max_extension=max_extension,
                            gap_threshold=gap_thr, gap_nul_thr=gap_nul_thr,
                            crop_end_thr=crop_end_thr, crop_end_win=crop_end_win,
                            crop_end_gap_thr=crop_end_gap_thr, crop_end_gap_win=crop_end_gap_win,
                            start_patterns=start_patterns, end_patterns=end_patterns,
                            mini_orf=mini_orf, define_boundary_win=check_extension_win)

                    if inner_find_boundary_result == "Short_sequence":
                        continue
                    elif inner_find_boundary_result:
                        all_inner_skipped = False
                    elif not inner_find_boundary_result:  # This means the errors happen in the function
                        return

                else:  # If need more cluster, clean_and_cluster_MSA will return two values
                    inner_cluster_bed_files_list = inner_cluster_MSA_result

                    for j in range(len(inner_cluster_bed_files_list)):
                        inner_inner_find_boundary_result = find_boundary_and_crop(
                                inner_cluster_bed_files_list[j], genome_file, MSA_dir,
                                pfam_dir, seq_name, hmm,
                                cons_threshold=cons_thr, ext_threshold=ext_thr,
                                ex_step_size=ex_step, max_extension=max_extension,
                                gap_threshold=gap_thr, gap_nul_thr=gap_nul_thr,
                                crop_end_thr=crop_end_thr, crop_end_win=crop_end_win,
                                crop_end_gap_thr=crop_end_gap_thr, crop_end_gap_win=crop_end_gap_win,
                                start_patterns=start_patterns, end_patterns=end_patterns,
                                mini_orf=mini_orf, define_boundary_win=check_extension_win)
                        if inner_inner_find_boundary_result == "Short_sequence":
                            continue
                        elif inner_inner_find_boundary_result:
                            all_inner_skipped = False
                        elif not inner_inner_find_boundary_result:  # This means the errors happen in the function
                            return

            # Check the flag after the loop. If all inner clusters were skipped, write to skipped_file
            if all_inner_skipped:
                with open(skipped_file, "a") as f:
                    f.write(seq_name + f"\tblast_less\tcluster_blast_smaller_than_{min_seq_num} or too short\n")
                    click.echo(
                        f"\n{seq_name} is skipped due to sequence number in second round each cluster is "
                        f"smaller than {min_seq_num} or the sequence is too short\n")

                # If all this sequence is skipped remove all files contain this name
                if not keep_intermediate:
                    remove_file_object = SequenceManipulator()
                    remove_file_object.remove_files_with_start_pattern(MSA_dir, seq_name)
                    remove_file_object.remove_files_with_start_pattern(elongation_dir, seq_name)
                return

    except Exception as e:
        click.echo(f"Error during boundary finding and cropping for sequence: {seq_name}. Error: {str(e)}")
        traceback.print_exc()
        return

    # After all processing is done, write the name of the file to the progress file
    with open(progress_file, "a") as f:
        f.write(seq_name + "\n")

    click.echo(f"Finished {seq_name}")

    # If all this sequence is finished remove all files contain this name
    if not keep_intermediate:
        remove_file_object = SequenceManipulator()
        remove_file_object.remove_files_with_start_pattern(MSA_dir, seq_name)
        remove_file_object.remove_files_with_start_pattern(elongation_dir, seq_name)
    # Check the sequences numbers in Finished_sequence_name.txt and Skipped_sequence_name.txt
    # and give how many sequences are left
    # Read and count sequences from Finished_sequence_name.txt
    with open(progress_file, 'r') as file:
        progress_lines = file.readlines()
        progress_file_count = len(progress_lines)

    # Read and count sequences from Skipped_sequence_name.txt
    with open(skipped_file, 'r') as file:
        skipped_lines = file.readlines()
        skipped_file_count = len(skipped_lines)

    # Calculate the total count
    processed_count = progress_file_count + skipped_file_count

    # Calculate sequences number that hasn't been processed by TE Trimmer
    rest_sequence = single_fasta_n - processed_count

    click.echo(f"{progress_file_count} sequences have been successfully processed, "
               f"{skipped_file_count} sequences were skipped, "
               f"{rest_sequence} sequences need to be processed\n")

#####################################################################################################
# Code block: Import json species_config file and define the default parameters
#####################################################################################################


species_config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'species_config.json')

# Load the JSON configuration file
with open(species_config_path, "r") as config_file:
    species_config = json.load(config_file)


@click.command(context_settings=dict(max_content_width=100))
@click.option('--input_file', '-i', required=True, type=str,
              help='Input fasta file. TE consensus sequences')
@click.option('--output_dir', '-o', default=os.getcwd(), type=str,
              help='Output directory. Default: current directory')
@click.option('--genome_file', '-g', required=True, type=str,
              help='Provide the genome file path')
@click.option('--species', '-s', required=True, default='fungi', type=click.Choice(species_config.keys()),
              help='Select the species for which you want to run TE Trimmer')
@click.option('--continue_analysis', default=False, is_flag=True,
              help='If --continue_analysis is provided on the command line, TE Trimmer will continue the analysis '
                   'based on the existing data. Otherwise, it will overlap the existing files. Default: False')
@click.option('--cd_hit_merge', default=False, is_flag=True,
              help='If --cd_hit_merge, the input file will be merged first')
@click.option('--genome_anno', default=False, is_flag=True,
              help='If to perform final genome TE annotation by RepeatMasker')
@click.option('--hmm', default=False, is_flag=True,
              help='If to generate hmm files')
@click.option('--keep_intermediate', default=False, is_flag=True,
              help='Use this option if you want to keep the raw files. ATTENTION: Many files will be generated')
@click.option('--pfam_dir', default=None, type=str,
              help="Pfam database directory. Leave this option when you don't have Pfam database, "
                   "TE Trimmer will download automatically")
@click.option('--cons_thr', type=float,
              help='Threshold used for the final consensus sequence generation. Default: 0.8')
@click.option('--max_msa_lines', type=int,
              help='Set the maximum sequences number for multiple sequence alignment. Default: 100')
@click.option('--top_mas_lines', type=int,
              help='When the sequence number of multiple sequence alignment (MSA) is greater than '
                   '"max_mas_lines". It will order sequences by length and choose '
                   '"top_msa_lines" number of sequences. Then randomly choose the rest number of sequences '
                   'Default: 100')
@click.option('--min_seq_num', type=int,
              help='The minimum sequence number for each multiple sequence alignment. Default: 10')
@click.option('--min_blast_len', type=int,
              help='The minimum hit sequence length for blast. Default: 150')
@click.option('--max_cluster_num', type=int,
              help='The maximum cluster number for each multiple sequence alignment. Each multiple sequence alignment '
                   'can be divided into different clusters. TE Trimmer will sort cluster by sequence number and choose'
                   'the top --max_cluster_num of clusters for the further analysis. Default: 2')
@click.option('--min_el', type=int,
              help='If the query sequence length is smaller than the given number, TE Trimmer will elongate the query '
                   'sequence. This can help identify the intact TE. Default: 2000')
@click.option('--min_el_dna', type=int,
              help='Same like --min_el option. But set different given numbers for DNA elements. Default: 500')
@click.option('--min_el_sine', type=int,
              help='Same like --min_el option. But set different given number for SINE element. Default: 100')
@click.option('--ext_thr', type=float,
              help="threshold used for define the extension extent. The smaller number means it become easier "
                   "to have a final longer extension for each side of the sequence. Default: 0.7")
@click.option('--ex_step', type=int,
              help='Number of nucleotides will be added to the left or right side of multiple sequence alignment. '
                   'TE_Trimmer will iteratively add --ex_step number of nucleotide until finding the boundary. '
                   'Default: 1000')
@click.option('--max_extension', type=int,
              help='The maximum extension number for the right and left side. For example, if --ex_step is 1000, '
                   'it can only add seven times to the MSA left or right side. Default: 7000')
@click.option('--gap_thr', type=float,
              help='If columns have a larger gap proportion than --gap_thr and the most common nucleotide proportion '
                   'in this column is less than --gap_nul_thr, this column will be removed. Default: 0.4')
@click.option('--gap_nul_thr', type=float,
              help='Set nucleotide proportion to decide if remove this column. Coupled with --gap_nul_thr option. '
                   'Default: 0.7')
@click.option('--crop_end_win', type=int,
              help='Window size used for crop end process. Coupled with --crop_end_thr option. Default: 20')
@click.option('--crop_end_thr', type=int,
              help='Crop end function will convert each nucleotide in MSA into proportion number. This function will '
                   'check from the beginning and end of each sequence from MSA by iteratively choosing a slide window '
                   'and sum up the proportion numbers. It will stop until the summary of proportion is larger than '
                   '--crop_end_thr. The nucleotide do not match --crop_end_thr will be converted to -. The recommended '
                   'number is 0.8 * --crop_end_win. Default: 16')
@click.option('--crop_end_gap_win', type=int,
              help='Define window size used to crop end by gap, coupled with --crop_end_gap_thr option. Default: 150')
@click.option('--crop_end_gap_thr', type=float,
              help='Crop end by gap function will check from the beginning and end of each sequence from MSA by '
                   'iteratively choosing a slide window. It will stop until the gap proportion in this slide window is '
                   'smaller than --crop_end_gap_thr. The nucleotide do not match the requirement will be converted to -'
                   ' Default: 0.1')
@click.option('--start_patterns', type=str,
              help='LTR elements will always start with fixed pattern. TE Trimmer will check if it starts with those '
                   'patterns. If not, it will seek around the start point, if the pattern is found, the start point '
                   'will be converted to there. Note: if you want to give multiple start patterns, separate them by comma. '
                   'Like: TG,TA,TC (No space between them). The order of the given patterns is matter. Default: TG')
@click.option('--end_patterns', type=str,
              help='LTR elements will always end with fixed pattern. TE Trimmer will check if it end with those '
                   'patterns. If not, it will seek around the end point, if the pattern is found, the start point will '
                   'be converted to there. Note: if you want to give multiple end patterns, separate them by comma. Like: '
                   'CA,TA,GA (No space between them). The order of the given patterns is matter. Default: AC')
@click.option('--mini_orf', type=int,
              help='Set the minimum ORF length that will be predicted by TE Trimmer. Default: 200')
@click.option('--check_extension_win', type=str,
              help='Define check windows size for extension. Deafault: 150')
@click.option('--num_threads', '-t', default=10, type=int,
              help='Threads numbers used for TE Trimmer. Default: 10')
# add the key parameters into click options
def main(input_file, genome_file, output_dir, continue_analysis, pfam_dir, min_blast_len, num_threads, max_msa_lines,
         top_mas_lines, min_seq_num, max_cluster_num, min_el, min_el_dna, min_el_sine, cons_thr, ext_thr, ex_step,
         max_extension, gap_thr, gap_nul_thr, crop_end_thr, crop_end_win, crop_end_gap_thr, crop_end_gap_win,
         start_patterns, end_patterns, mini_orf, species, check_extension_win, cd_hit_merge, genome_anno, hmm,
         keep_intermediate):
    """
        TE Trimmer is software that can replace transposable element (TE) manual curation. Two mandatory arguments are
        required. They are TE consensus file and genome file. TE Trimmer can automatically define proper extension size
        and find the right boundary.

        python ./path_to_TE_Trimmer_bin/main.py --input_file "your_TE_consensus_file_path" --genome_file "your_genome_file_path"

    """

    start_time = datetime.now()
    print(f"\nTE Trimmer started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    click.echo("TE Trimmer is running...... \n"
               "WARNING: This is a beta version, please send bugs to jqian@bio1.rwth-aachen.de ;)\n")

    #####################################################################################################
    # Code block: Define the default options according to the given species
    #####################################################################################################

    # Load species-specific default values from the JSON config
    default_values = species_config.get(species, {})

    # Override default values with command-line options if provided

    if cons_thr is None:
        cons_thr = default_values.get("cons_thr")

    if max_msa_lines is None:
        max_msa_lines = default_values.get("max_msa_lines")

    if top_mas_lines is None:
        top_mas_lines = default_values.get("top_mas_lines")

    if min_seq_num is None:
        min_seq_num = default_values.get("min_seq_num")

    if min_blast_len is None:
        min_blast_len = default_values.get("min_blast_len")

    if max_cluster_num is None:
        max_cluster_num = default_values.get("max_cluster_num")

    if min_el is None:
        min_el = default_values.get("min_el")

    if min_el_dna is None:
        min_el_dna = default_values.get("min_el_dna")

    if min_el_sine is None:
        min_el_sine = default_values.get("min_el_sine")

    if ext_thr is None:
        ext_thr = default_values.get("ext_thr")

    if ex_step is None:
        ex_step = default_values.get("ex_step")

    if max_extension is None:
        max_extension = default_values.get("max_extension")

    if gap_thr is None:
        gap_thr = default_values.get("gap_thr")

    if gap_nul_thr is None:
        gap_nul_thr = default_values.get("gap_nul_thr")

    if crop_end_thr is None:
        crop_end_thr = default_values.get("crop_end_thr")

    if crop_end_win is None:
        crop_end_win = default_values.get("crop_end_win")

    if crop_end_gap_thr is None:
        crop_end_gap_thr = default_values.get("crop_end_gap_thr")

    if crop_end_gap_win is None:
        crop_end_gap_win = default_values.get("crop_end_gap_win")

    if start_patterns is None:
        start_patterns = default_values.get("start_patterns")

    if end_patterns is None:
        end_patterns = default_values.get("end_patterns")

    if mini_orf is None:
        mini_orf = default_values.get("mini_orf")

    if check_extension_win is None:
        check_extension_win = default_values.get("check_extension_win")

    #####################################################################################################
    # Code block: Define input file, output directory, genome, check blast database
    #####################################################################################################

    # bin_py_path contains all classes and bash code
    # so.path.abspath(__file__) will return the current executable python file
    bin_py_path = os.path.dirname(os.path.abspath(__file__))

    # check if path exist otherwise create one
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    output_dir = os.path.abspath(output_dir)  # get absolute path

    # Check if output directory is empty when --continue_analysis is False
    if os.listdir(output_dir) and not continue_analysis:

        click.echo(
            f"WARNING: The output directory {output_dir} is not empty. Please empty your output directory or "
            f"choose another empty directory\n")

        # Stop the whole program when the output directory is not empty
        return

    # make a new folder for single fasta sequence
    single_file_dir = os.path.join(output_dir, "Single_fasta_files")
    if not os.path.exists(single_file_dir):
        os.mkdir(single_file_dir)

    # make a new folder for MSA
    MSA_dir = os.path.join(output_dir, "Multiple_sequence_alignment")
    if not os.path.exists(MSA_dir):
        os.mkdir(MSA_dir)

    # make a new folder for HMM file
    if hmm:
        hmm_dir = os.path.join(output_dir, "HMM_files")
        if not os.path.exists(hmm_dir):
            os.mkdir(hmm_dir)

    # make a new folder for elongation files
    elongation_dir = os.path.join(output_dir, "Elongation_folder")
    if not os.path.exists(elongation_dir):
        os.mkdir(elongation_dir)

    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"The fasta file {input_file} does not exist.")
    input_file = os.path.abspath(input_file)  # get input absolute path

    # check if genome file exist
    if not os.path.isfile(genome_file):
        raise FileNotFoundError(f"The genome fasta file {genome_file} does not exist.")
    genome_file = os.path.abspath(genome_file)  # get genome absolute path

    # Define the skipped_file, all single sequence that is skipped by TE Trimmer will be stored inside this file
    skipped_file = os.path.join(output_dir, "Skipped_sequence_name.txt")
    # Check and create skipped_file if it doesn't exist
    if not os.path.exists(skipped_file):
        with open(skipped_file, 'a'):
            pass

    # Define the progress_file, finished seq IDs will be stored here
    progress_file = os.path.join(output_dir, "Finished_sequence_name.txt")
    # Check and create progress_file if it doesn't exist
    if not os.path.exists(progress_file):
        with open(progress_file, 'a'):
            pass

    # If pfam database isn't provided, create pfam database at TE Trimmer software folder,
    # the database will be downloaded into there.
    # If the pfam_dir is given, but the database can't be found there, TE Trimmer will download pfam database there
    # and generate index file
    if pfam_dir is None:
        pfam_dir = os.path.join(os.path.dirname(bin_py_path), "pfam_database")
        if not os.path.exists(pfam_dir):
            os.mkdir(pfam_dir)

    prepare_pfam_database(pfam_dir)

    #####################################################################################################
    # Code block: Merge input file and generate single fasta file
    #####################################################################################################

    # Generate single files when continue_analysis is false
    if not continue_analysis:

        # Do cd-hit-est merge when cd_hit_merge is true and continue_analysis is false
        if cd_hit_merge:
            click.echo("\nTE Trimmer is merging input sequences, this might take some time.\n")
            cd_hit_merge_output = os.path.join(output_dir, f"{input_file}_cd_hit.fa")
            cd_hit_merge_object = SequenceManipulator()
            # Set lower identity threshold for the query, this can increase sensitive
            cd_hit_merge_object.cd_hit_est(input_file, cd_hit_merge_output, identity_thr=0.9,
                                           aL=0.9, aS=0.9, s=0.9, thread=num_threads)

            # Convert input_file to merged input_file
            input_file = cd_hit_merge_output
            click.echo("Merge finished.\n")

        # separate fasta to single files, if fasta header contain "/" or " " or "#" convert them to "_"
        separate_fasta = FastaSequenceSeparator(input_file, single_file_dir)
        separate_fasta.separate_sequences()  # call this function to separate to single fasta files

        # Calculate the total sequence number in single_file_dir
        single_fasta_n = len([f for f in os.listdir(single_file_dir) if f.endswith('.fasta')])
        click.echo(f"{single_fasta_n} sequences are detected from the input file.")

        # create new object to check blast database availability
        test_blast_db = SequenceManipulator()

        # check if blast database and genome length files are available, otherwise create them at the
        # same directory of genome file
        test_blast_db.check_database(genome_file)
        sequence_names = [f for f in os.listdir(single_file_dir) if f.endswith(".fasta")]

    else:

        # Check if progress_file and skipped_file are both empty
        if not os.listdir(single_file_dir):
            click.echo("\nWARNING: TE Trimmer can't do continue analysis, please make sure the output directory is same"
                       "with your previous analysis.")
            return

        else:
            click.echo("\nTE Trimmer will continue to analyze based on previous results.\n")
            single_fasta_n = len([f for f in os.listdir(single_file_dir) if f.endswith('.fasta')])

            # Check which sequences have already been processed
            progress_sequences = check_progress_file(progress_file)
            skipped_sequences = check_progress_file(skipped_file)
            complete_sequences = progress_sequences + skipped_sequences

            # Filter out already complete sequences from the total sequences
            sequence_names = [f for f in os.listdir(single_file_dir) if
                              f.endswith(".fasta") and f not in complete_sequences]

    #####################################################################################################
    # Code block: Enable multiple threads
    #####################################################################################################
    analyze_sequence_params = [
        (seq_name, single_file_dir, genome_file, MSA_dir, min_blast_len, min_seq_num, max_msa_lines,
         top_mas_lines, max_cluster_num, min_el, min_el_dna, min_el_sine, cons_thr, ext_thr, ex_step,
         max_extension, gap_thr, gap_nul_thr, crop_end_thr, crop_end_win, crop_end_gap_thr, crop_end_gap_win,
         start_patterns, end_patterns, output_dir, pfam_dir, mini_orf, elongation_dir, single_fasta_n, hmm,
         check_extension_win, keep_intermediate
         ) for seq_name in sequence_names]

    # Using a ProcessPoolExecutor to run the function in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        executor.map(analyze_sequence_helper, analyze_sequence_params)

    #####################################################################################################
    # Code block: Check processed sequence number and rename final consensus file
    #####################################################################################################

    final_con_file = os.path.join(output_dir, "TE_Trimmer_consensus.fasta")

    cd_hit_merge_output_final = None
    if os.path.exists(final_con_file):

        # Read the file into memory
        with open(final_con_file, 'r') as file:
            file_contents = file.read()

        # Make the replacements
        file_contents = file_contents.replace('___', '#').replace('__', '/')

        # Write the modified contents back to the file
        with open(final_con_file, 'w') as file:
            file.write(file_contents)

        # Do cd-hit-est merge when finish all sequence
        cd_hit_merge_output_final = os.path.join(output_dir, "TE_Trimmer_consensus_merged.fasta")
        cd_hit_merge_object = SequenceManipulator()
        # According to 80-80-80 rule to filter final consensus sequences
        cd_hit_merge_object.cd_hit_est(final_con_file, cd_hit_merge_output_final, identity_thr=0.8,
                                       aL=0.8, aS=0.8, s=1, thread=num_threads)

        #

    # At the end of the program, check if all sequences have been processed
    with open(progress_file, 'r') as file:
        progress_lines = file.readlines()
        progress_file_count = len(progress_lines)

    # Read and count sequences from Skipped_sequence_name.txt
    with open(skipped_file, 'r') as file:
        skipped_lines = file.readlines()
        skipped_file_count = len(skipped_lines)

    # Calculate the total count
    processed_count = progress_file_count + skipped_file_count

    if processed_count == single_fasta_n:
        click.echo(f"All sequences have been processed! {skipped_file_count} sequences are skipped")

        if not keep_intermediate:
            # Remove all single files when all the sequences are processed
            shutil.rmtree(single_file_dir)

    else:
        remaining = single_fasta_n - processed_count
        click.echo(f"{remaining} sequences have not been processed.")

    # Delete MSA_dir and elongation_dir if they are empty
    if not os.listdir(MSA_dir):
        os.rmdir(MSA_dir)
    if not os.listdir(elongation_dir):
        os.rmdir(elongation_dir)

    # If 95% of the query sequences are processed, RepeatMasker is allowed to be performed
    if processed_count >= single_fasta_n*0.9:
        # Run RepeatMasker
        if genome_anno and cd_hit_merge_output_final:
            click.echo("TE Trimmer is performing whole genome TE annotation by RepeatMasker")
            # make a new folder for RepeatMasker output
            repeatmasker_dir = os.path.join(output_dir, "RepeatMasker_result")
            if not os.path.exists(repeatmasker_dir):
                os.mkdir(repeatmasker_dir)
            genome_anno_object = SequenceManipulator()
            genome_anno_result = \
                genome_anno_object.repeatmasker(genome_file, cd_hit_merge_output_final, repeatmasker_dir, thread=num_threads)
            if genome_anno_result:
                click.echo("Finished whole genome TE annotation by RepeatMasker")

    else:
        click.echo("Less than 90% of the query sequences processed, TE Trimmer can't perform whole genome TE annotation")

    end_time = datetime.now()
    duration = end_time - start_time

    # Remove microseconds from the duration
    duration_without_microseconds = timedelta(days=duration.days, seconds=duration.seconds)

    print(f"\nTE Trimmer finished at {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    print(f"TE Trimmer runtime was {duration_without_microseconds}")


# The following is necessary to make the script executable, i.e., python myscript.py.
if __name__ == '__main__':
    main()
