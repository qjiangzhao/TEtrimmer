import os
import concurrent.futures
import traceback
from Bio import SeqIO
from Class_separate_fasta import FastaSequenceSeparator
from Class_blast_extension_mafft import SequenceManipulator
from Class_bed_filter import BEDFile
from Function_def_boundary_and_crop import find_boundary_and_crop
from Function_clean_and_clauster_MSA import clean_and_cluster_MSA
from Class_elongate_query import SequenceElongation
from Class_pfam_scan_and_plot import prepare_pfam_database
from Class_cd_hit_est import CDHitEst
import click


class DirectoryNotEmptyError(Exception):
    pass


def analyze_sequence_helper(params):
    return analyze_sequence(*params)


def analyze_sequence(seq_name, Single_files_dir, genome_file, MSA_dir, min_blast_len, min_seq_num, max_msa_lines,
                     top_mas_lines, max_cluster_num, min_el, min_el_dna, min_el_sine, cons_thr, ex_step,
                     max_extension, gap_thr, gap_nul_thr, crop_end_thr, crop_end_win, crop_end_gap_thr,
                     crop_end_gap_win, start_patterns, end_patterns, output_dir, pfam_dir, mini_orf, elongation_dir,
                     single_fasta_n):
    #####################################################################################################
    # Code block: Elongate query sequence when they are too short
    #####################################################################################################

    try:
        # Get query fasta file path
        seq_file = os.path.join(Single_files_dir, seq_name)

        # Define the skipped_file, all single sequence that is skipped by TE Trimmer will be stored inside this file
        skipped_file = os.path.join(output_dir, "Skipped_file.txt")
        # Check and create skipped_file if it doesn't exist
        if not os.path.exists(skipped_file):
            with open(skipped_file, 'a'):
                pass

        # If sequence length is less than min_length_elongation. Elongation process will be executed.
        # Because short query sequence won't enable efficient MSA cluster
        seq_record = SeqIO.read(seq_file, "fasta")
        seq_file_len = len(seq_record.seq)
        min_length_elongation = min_el

        # TODO add more elements for different length
        # Due to DNA element will be much shorter than LTR and LINE elements, set min_length_elongation to 500.
        if "DNA" in seq_name:
            min_length_elongation = min_el_dna

        # The average length of SINE element is around 500 bp, set smaller min_length_elongation number.
        if "SINE" in seq_name:
            min_length_elongation = min_el_sine

        if seq_file_len < min_length_elongation:

            seq_elongation_object = SequenceElongation(seq_file, genome_file, elongation_dir)
            seq_elongation = seq_elongation_object.elongate_query_seq()

            if seq_elongation:
                seq_file = seq_elongation

            # When seq_elongation is false, skip this sequence
            elif not seq_elongation:
                with open(skipped_file, "a") as f:
                    f.write(seq_name + "\n")
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
                f.write(seq_name + "\n")
            click.echo(f"{seq_name} is skipped due to blast hit number is 0\n")
            return

        # check if blast hit number is smaller than 10
        elif seq_blast.blast_hits_count != 0 and seq_blast.blast_hits_count <= min_seq_num:
            with open(skipped_file, "a") as f:
                f.write(seq_name + "\n")
            click.echo(f"{seq_name} is skipped due to blast hit number is smaller than {min_seq_num}\n")
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
        fasta_out_flank_file = seq_blast.extract_fasta(bed_out_filter_file, genome_file, MSA_dir, left_ex=0,
                                                       right_ex=0)

        # Return False when cluster number is 0. Return True when divergent column number is smaller than 100
        # Otherwise it will return the subset bed and alignment file
        cluster_MSA_result = clean_and_cluster_MSA(fasta_out_flank_file, bed_out_filter_file, MSA_dir,
                                                   gap_threshold=0.8, clean_column_threshold=0.08,
                                                   min_length_num=min_seq_num, cluster_num=max_cluster_num)
    except Exception as e:
        click.echo(
            f"Error during processing lines, extracting fasta, or clustering MSA for sequence: {seq_name}. Error: {str(e)}")
        traceback.print_exc()
        return

    try:
        # cluster false means no cluster, TE Trimmer will skip this sequence.
        if cluster_MSA_result is False:
            with open(skipped_file, "a") as f:
                f.write(seq_name + "\n")
            click.echo(f"{seq_name} is skipped due to sequence number in each cluster is smaller than {min_seq_num}\n")
            return
        elif cluster_MSA_result is True:
            find_boundary_and_crop(bed_out_filter_file, genome_file, MSA_dir,
                                   pfam_dir, seq_name, cons_threshold=cons_thr,
                                   ex_step_size=ex_step, max_extension=max_extension,
                                   gap_threshold=gap_thr, gap_nul_thr=gap_nul_thr,
                                   crop_end_thr=crop_end_thr, crop_end_win=crop_end_win,
                                   crop_end_gap_thr=crop_end_gap_thr, crop_end_gap_win=crop_end_gap_win,
                                   start_patterns=start_patterns, end_patterns=end_patterns, mini_orf=mini_orf)
        else:
            cluster_pattern_alignment_list, cluster_bed_files_list = cluster_MSA_result

            # cluster_pattern_alignment_list have the same index with cluster_bed_files_list
            for i in range(len(cluster_pattern_alignment_list)):
                inner_cluster_MSA_result = clean_and_cluster_MSA(cluster_pattern_alignment_list[i],
                                                                 cluster_bed_files_list[i], MSA_dir,
                                                                 gap_threshold=0.8, clean_column_threshold=0.15,
                                                                 min_length_num=min_seq_num,
                                                                 cluster_num=max_cluster_num)
                if inner_cluster_MSA_result is False:
                    continue
                elif inner_cluster_MSA_result is True:
                    find_boundary_and_crop(cluster_bed_files_list[i], genome_file, MSA_dir,
                                           pfam_dir, seq_name, cons_threshold=cons_thr,
                                           ex_step_size=ex_step, max_extension=max_extension,
                                           gap_threshold=gap_thr, gap_nul_thr=gap_nul_thr,
                                           crop_end_thr=crop_end_thr, crop_end_win=crop_end_win,
                                           crop_end_gap_thr=crop_end_gap_thr, crop_end_gap_win=crop_end_gap_win,
                                           start_patterns=start_patterns, end_patterns=end_patterns,
                                           mini_orf=mini_orf)
                else:
                    inner_cluster_pattern_alignment_list, inner_cluster_bed_files_list = inner_cluster_MSA_result

                    for j in range(len(inner_cluster_bed_files_list)):
                        find_boundary_and_crop(inner_cluster_bed_files_list[j], genome_file, MSA_dir,
                                               pfam_dir, seq_name, cons_threshold=cons_thr,
                                               ex_step_size=ex_step, max_extension=max_extension,
                                               gap_threshold=gap_thr, gap_nul_thr=gap_nul_thr,
                                               crop_end_thr=crop_end_thr, crop_end_win=crop_end_win,
                                               crop_end_gap_thr=crop_end_gap_thr, crop_end_gap_win=crop_end_gap_win,
                                               start_patterns=start_patterns, end_patterns=end_patterns,
                                               mini_orf=mini_orf)
    except Exception as e:
        click.echo(f"Error during boundary finding and cropping for sequence: {seq_name}. Error: {str(e)}")
        traceback.print_exc()

    # After all processing is done, write the name of the file to the progress file
    progress_file = os.path.join(output_dir, "progress_file.txt")
    with open(progress_file, "a") as f:
        f.write(seq_name + "\n")

    click.echo(f"Finished {seq_name}")

    # Check the sequences numbers in progress_file.txt and Skipped_file.txt and give how many sequences left
    # Read and count sequences from progress_file.txt
    with open(progress_file, 'r') as file:
        progress_lines = file.readlines()
        progress_file_count = len(progress_lines)

    # Read and count sequences from skipped_file.txt
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


@click.command(context_settings=dict(max_content_width=100))
@click.option('--input_file', '-i', required=True, type=str,
              help='Input fasta file. TE consensus sequences')
@click.option('--output_dir', '-o', default=os.getcwd(), type=str,
              help='Output directory. Default: current directory')
@click.option('--genome_file', '-g', required=True, type=str,
              help='Provide the genome file path')
@click.option('--continue_analysis', default=False, is_flag=True,
              help='If --continue_analysis is provided on the command line, TE Trimmer will continue the analysis '
                   'based on the existing data. Otherwise, it will overlap the existing files')
@click.option('--pfam_dir', default=None, type=str,
              help="Pfam database directory. Leave this option when you don't have Pfam database, "
                   "TE Trimmer will download automatically")
@click.option('--max_msa_lines', default=100, type=int,
              help='Set the maximum sequences number for multiple sequence alignment. Default: 100')
@click.option('--min_seq_num', default=10, type=int,
              help='The minimum sequence number for each multiple sequence alignment. Default: 10')
@click.option('--min_blast_len', default=200, type=int,
              help='Blast sequence length lower than --min_blast_len will be ignored')
@click.option('--top_mas_lines', default=100, type=int,
              help='When the sequence number of multiple sequence alignment (MSA) is greater than '
                   '"max_mas_lines". It will order sequences by length and choose '
                   '"top_msa_lines" number of sequences. Then randomly choose the rest number of sequences '
                   'Default: 100')
@click.option('--max_cluster_num', default=2, type=int,
              help='The maximum cluster number for each multiple sequence alignment. Each multiple sequence alignment '
                   'can be divided into different groups. The group is called cluster')
@click.option('--min_el', default=2000, type=int,
              help='If the query sequence length is smaller than the given number, 1500bp will be added to the left '
                   'and right sides of this sequence. Meanwhile, redundant elongation will be cleaned')
@click.option('--min_el_dna', default=500, type=int,
              help='Same like --min_el option. But set different given numbers for DNA elements. Default: 500')
@click.option('--min_el_sine', default=200, type=int,
              help='Same like --min_el option. But set different given number for SINE element. Default: 200')
@click.option('--cons_thr', default=0.7, type=float,
              help='Threshold used for consensus sequence generation for multiple sequence alignment')
@click.option('--ex_step', default=1000, type=int,
              help='Number of nucleotides will be added to the left or right side of multiple sequence alignment. '
                   'TE_Trimmer will iteratively add --ex_step number of nucleotide until finding the boundary')
@click.option('--max_extension', default=7000, type=int,
              help='The maximum extension number for the right and left side. For example, if --ex_step is 1000, '
                   'it can only add seven times to the MSA left or right side')
@click.option('--gap_thr', default=0.4, type=float,
              help='If columns have a larger gap proportion than --gap_thr and the most common nucleotide proportion in '
                   'this column is less than --gap_nul_thr, this column will be removed')
@click.option('--gap_nul_thr', default=0.7, type=float,
              help='Set nucleotide proportion to decide if remove this column. Coupled with --gap_nul_thr option')
@click.option('--crop_end_thr', default=16, type=int,
              help='Crop end function will convert each nucleotide in MSA into proportion number. This function will '
                   'check from the beginning and end of each sequence from MSA by iteratively choosing a slide window '
                   'and sum up the proportion numbers. It will stop until the summary of proportion is larger than '
                   '--crop_end_thr. The nucleotide do not match --crop_end_thr will be converted to -')
@click.option('--crop_end_win', default=20, type=int,
              help='Window size used for crop end process. Coupled with --crop_end_thr option')
@click.option('--crop_end_gap_thr', default=0.05, type=float,
              help='Crop end by gap function will check from the beginning and end of each sequence from MSA by '
                   'iteratively choosing a slide window. I will stop until the gap proportion in this slide window is '
                   'smaller than --crop_end_gap_thr. The nucleotide do not match the requirement will be converted to -')
@click.option('--crop_end_gap_win', default=300, type=int,
              help='Define window size used to crop end by gap, coupled with --crop_end_gap_thr option')
@click.option('--start_patterns', default='TG', type=str,
              help='LTR elements will always start with fixed pattern. TE Trimmer will check if it starts with those '
                   'patterns. If not, it will seek around the start point, if the pattern is found, the start point '
                   'will be converted to there. Note: if you want to give multiple start patterns, separate them by comma. '
                   'Like: TG,TA,TC (No space between them). The order of the given patterns is important')
@click.option('--end_patterns', default='CA', type=str,
              help='LTR elements will always end with fixed pattern. TE Trimmer will check if it end with those '
                   'patterns. If not, it will seek around the end point, if the pattern is found, the start point will '
                   'be converted to there. Note: if you want to give multiple end patterns, separate them by comma. Like: '
                   'CA,TA,GA (No space between them). The order of the given patterns is important')
@click.option('--mini_orf', default=200, type=int,
              help='Set the minimum ORF length that will be predicted by TE Trimmer')
@click.option('--num_threads', '-t', default=10, type=int,
              help='Threads numbers used for TE Trimmer. Default: 10')
# add the key parameters into click options
def main(input_file, genome_file, output_dir, continue_analysis, pfam_dir, min_blast_len, num_threads, max_msa_lines,
         top_mas_lines, min_seq_num, max_cluster_num, min_el, min_el_dna, min_el_sine, cons_thr, ex_step,
         max_extension, gap_thr, gap_nul_thr, crop_end_thr, crop_end_win, crop_end_gap_thr, crop_end_gap_win,
         start_patterns, end_patterns, mini_orf):
    """
        TE Trimmer is software that can replace transposable element (TE) manual curation. Two mandatory arguments are
        required. They are TE consensus file and genome file. TE Trimmer can automatically define proper extension size
        and find the right boundary.

        python ./path_to_TE_Trimmer_bin/main.py --input_file "your_TE_consensus_file_path" --genome_file "your_genome_file_path"

    """
    click.echo("TE Trimmer is running......")

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
    try:
        if os.listdir(output_dir) and not continue_analysis:
            raise DirectoryNotEmptyError

    except DirectoryNotEmptyError:
        click.echo(
            f"The directory {output_dir} is not empty. Please choose a different output directory")
        return

    # make a new folder for single fasta sequence
    Single_files_dir = os.path.join(output_dir, "Single_fasta_files")
    if not os.path.exists(Single_files_dir):
        os.mkdir(Single_files_dir)

    # make a new folder for MSA
    MSA_dir = os.path.join(output_dir, "Multiple_sequence_alignment")
    if not os.path.exists(MSA_dir):
        os.mkdir(MSA_dir)

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

    # Create a folder at the at directory with genome for blast database
    genome_file_path = os.path.dirname(genome_file)
    blast_database_dir = os.path.join(genome_file_path, "TE_Trimmer_blast_db")
    if not os.path.exists(blast_database_dir):
        os.mkdir(blast_database_dir)

    # If pfam database isn't provided, create pfam database at TE Trimmer software folder,
    # the database will be downloaded into there.
    # If the pfam_dir is given, but the database can't be found there, TE Trimmer will download pfam database there
    # and generate index file
    if pfam_dir is None:
        pfam_dir = os.path.join(os.path.dirname(bin_py_path), "pfam_database")
        if not os.path.exists(pfam_dir):
            os.mkdir(pfam_dir)

    prepare_pfam_database(pfam_dir)

    # separate fasta to single files, if fasta header contain "/" or " " or "#" convert them to "_"
    separate_fasta = FastaSequenceSeparator(input_file, Single_files_dir)
    separate_fasta.separate_sequences()  # call this function to separate to single fasta files

    # Calculate the total sequence number in Single_files_dir
    single_fasta_n = len([f for f in os.listdir(Single_files_dir) if f.endswith('.fasta')])
    click.echo(f"{single_fasta_n} sequences are detected from the input file")

    # create new object just use to check blast database availability
    test_blast_db = SequenceManipulator()

    # check if blast database and genome length files are available, otherwise create genome length file at the
    # same directory of genome and create blast database at blast_database_dir
    test_blast_db.check_database(genome_file)

    #####################################################################################################
    # Code block: Start analysis based on the existed files
    #####################################################################################################

    def check_progress_file(progress_file_path):
        completed_sequences = []

        try:
            # If the progress file exists, open the file and read the list of completed sequences
            with open(progress_file_path, 'r') as progress_file:
                completed_sequences = progress_file.readlines()

            # Strip newline characters and remove duplicates
            completed_sequences = list(set([sequence.strip() for sequence in completed_sequences]))
        except FileNotFoundError:
            click.echo("progress_file.txt isn't found, please change the -o option to the right output directory.")
            raise
        except Exception as e:
            click.echo(f"An error occurred while reading the progress file: {e}")

        return completed_sequences

    progress_file = os.path.join(output_dir, "progress_file.txt")
    skipped_file = os.path.join(output_dir, "Skipped_file.txt")

    if continue_analysis:

        # Check which sequences have already been processed
        completed_sequences = check_progress_file(progress_file)
        skipped_sequences = check_progress_file(skipped_file)

        processed_sequences = completed_sequences + skipped_sequences

        # Filter out already processed sequences from the sequence names to process
        sequence_names = [f for f in os.listdir(Single_files_dir) if
                          f.endswith(".fasta") and f not in processed_sequences]
    else:
        # Delete the old progress file when the new round TE Trimmer is executed
        if os.path.exists(progress_file):
            os.remove(progress_file)
        sequence_names = [f for f in os.listdir(Single_files_dir) if f.endswith(".fasta")]

    #####################################################################################################
    # Code block: Enable multiple threads
    #####################################################################################################
    analyze_sequence_params = [
        (seq_name, Single_files_dir, genome_file, MSA_dir, min_blast_len, min_seq_num, max_msa_lines,
         top_mas_lines, max_cluster_num, min_el, min_el_dna, min_el_sine, cons_thr, ex_step,
         max_extension, gap_thr, gap_nul_thr, crop_end_thr, crop_end_win, crop_end_gap_thr, crop_end_gap_win,
         start_patterns, end_patterns, output_dir, pfam_dir, mini_orf, elongation_dir, single_fasta_n
         ) for seq_name in sequence_names]

    # Using a ProcessPoolExecutor to run the function in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        executor.map(analyze_sequence_helper, analyze_sequence_params)

    #####################################################################################################
    # Code block: Check processed sequence number
    #####################################################################################################

    # At the end of the program, check if all sequences have been processed
    with open(progress_file, 'r') as file:
        progress_lines = file.readlines()
        progress_file_count = len(progress_lines)

    # Read and count sequences from skipped_file.txt
    with open(skipped_file, 'r') as file:
        skipped_lines = file.readlines()
        skipped_file_count = len(skipped_lines)

    # Calculate the total count
    processed_count = progress_file_count + skipped_file_count

    if processed_count == single_fasta_n:
        click.echo(f"All sequences have been processed! {skipped_file_count} sequences are skipped")

    else:
        remaining = single_fasta_n - processed_count
        click.echo(f"{remaining} sequences have not been processed.")

    #####################################################################################################
    # Code block: Run cd-hit-est to merge final consensus sequence and prepare for proof annotation
    #####################################################################################################





# The following is necessary to make the script executable, i.e., python myscript.py.
if __name__ == '__main__':
    main()
