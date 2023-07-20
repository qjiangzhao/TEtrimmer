import os
import concurrent.futures
import traceback
from Bio import SeqIO
from Class_separate_fasta import FastaSequenceSeparator
from Class_blast_extension_mafft import SequenceManipulator
from Class_bed_filter import BEDFile
from Function_def_boundary_and_crop import find_boundary_and_crop
from Function_clean_and_clauster_MSA import clean_and_cluster_MSA
from Class_move_files import FileCopier
from Class_elongate_query import SequenceElongation

#####################################################################################################
# Code block: Define input file, output directory, genome, check blast database
#####################################################################################################

# bin_py_path contains all classes and bash code
# so.path.abspath(__file__) will return the current executable python file
bin_py_path = os.path.dirname(os.path.abspath(__file__))

output_dir = "/work/ur376715/TE_analysis/Software_develpment_for_TE_curation/version_6/output_dir"
# check if path exist otherwise create one
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)
output_dir = os.path.abspath(output_dir)  # get absolute path

# make a new foler for single fasta sequence
Single_files_dir = os.path.join(output_dir, "Single_fasta_files")
if not os.path.exists(Single_files_dir):
    os.mkdir(Single_files_dir)

# make a new foler for MSA
MSA_dir = os.path.join(output_dir, "Multiple_sequence_alignment")
if not os.path.exists(MSA_dir):
    os.mkdir(MSA_dir)

input_fasta_file = "/work/ur376715/TE_analysis/Software_develpment_for_TE_curation/version_6/EDTA_sequence.fasta"

if not os.path.isfile(input_fasta_file):
    raise FileNotFoundError(f"The fasta file {input_fasta_file} does not exist.")
input_fasta_dir = os.path.dirname(os.path.abspath(input_fasta_file))  # get genome path
input_fasta_n = os.path.basename(os.path.abspath(input_fasta_file))

genome_file = "/work/ur376715/TE_analysis/TE_manual_curation/Genome_BLAST_database_for_TE_manual_curation/bgh_dh14_v4.fa"
# check if genome file exist
if not os.path.isfile(genome_file):
    raise FileNotFoundError(f"The fasta file {genome_file} does not exist.")
genome_dir = os.path.dirname(os.path.abspath(genome_file))  # get genome path
genome_n = os.path.basename(os.path.abspath(genome_file))

pfam_dir =  r"/hpcwork/ur376715/Data_base/Pfam_for_scan_pl"
# Check if pfam directory exit
if not os.path.exists(pfam_dir):
    os.mkdir(pfam_dir)

# separate fasta to single files, if fasta header contain "/" or " " convert them to "_"
Separate_fasta = FastaSequenceSeparator(input_fasta_file, Single_files_dir)
Separate_fasta.separate_sequences()  # call this function to separate to single fasta files

# create new object just use to check blast database availability
test_blast_db = SequenceManipulator()

# check if blast database and genome length files are available, otherwise create it at the same directory of genome
test_blast_db.check_database(genome_file)


#####################################################################################################
# Code block: Function do TE annotation
#####################################################################################################

def analyze_sequence(seq_name):

    try:
        # The os.path.abspath function returns the absolute path of a file or directory.
        seq_file = os.path.join(Single_files_dir, seq_name)

        # If sequence length is less than min_length_elongation. Elongation process will be executed.
        seq_record = SeqIO.read(seq_file, "fasta")
        seq_file_len = len(seq_record.seq)
        min_length_elongation = 2000

        # Due to DNA element will be much shorter than LTR and LINE elements, set min_length_elongation to 800.
        if "DNA" in seq_name:
            min_length_elongation = 800

        # The average length of SINE element is around 500 bp, set smaller min_length_elongation number.
        if "SINE" in seq_name:
            min_length_elongation = 100

        if seq_file_len < min_length_elongation:

            seq_elongation_object = SequenceElongation(seq_file, genome_file, output_dir)
            seq_elongation = seq_elongation_object.elongate_query_seq()

            if seq_elongation:
                seq_file = seq_elongation



        # run blast for each single fasta file and return a bed file absolute path
        seq_blast = SequenceManipulator()  # this class contain many functions to do blast and filter blast result
        bed_out_file_dup = seq_blast.blast(seq_file, genome_file, MSA_dir)

    except Exception as e:
        print(f"Error while running blast for sequence: {seq_name}. Error: {str(e)}")
        return

    try:
        # check if blast hit number is equal 0, then skip this sequence
        if seq_blast.blast_hits_count == 0:
            return
        # check if blast hit number is smaller than 10
        elif seq_blast.blast_hits_count != 0 and seq_blast.blast_hits_count <= 10:
            return  # when blast hit number is smaller than 10, code will execute next fasta file
            # with less hits numbers can also refer to young TE, but less blast hit number will hamper to generate
            # consensus sequence by multiple sequence alignment method.
            #TODO check if blast hit lengths occupied 80% of query if so, check if query have LTR and TE domain
            #TODO if it fits critium, keep it to final TE library but label with single copy TE
            #TODO write this fasta name to a file and do 4000 extension

    except Exception as e:
        print(f"Error while checking uniqueness for sequence: {seq_name}. Error: {str(e)}")
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
        bed_out_filter_file = bed_out_filter.process_lines(MSA_dir, threshold=100, top_longest_lines_count=100)

        # extract fast from bed_out_filter_file
        # return fasta_out_flank_file absolute path
        # because have to group MSA the first round extend for left and right side are both 0
        fasta_out_flank_file = seq_blast.extract_fasta(bed_out_filter_file, genome_file, MSA_dir, left_ex=0, right_ex=0)

        # returen False when cluster number is 0, otherwise, return True when divergent column number is smaller than 100
        # otherwise it will return the subset bed and alignment file
        cluster_MSA_result = clean_and_cluster_MSA(fasta_out_flank_file, bed_out_filter_file, MSA_dir,
                                                   gap_threshold=0.8, clean_column_threshold=0.08)
    except Exception as e:
        print(
            f"Error during processing lines, extracting fasta, or clustering MSA for sequence: {seq_name}. Error: {str(e)}")
        traceback.print_exc()
        return

    try:
        if cluster_MSA_result is False:
            return
        elif cluster_MSA_result is True:
            find_boundary_and_crop(bed_out_filter_file, genome_file, MSA_dir, pfam_dir=pfam_dir, left_ex=1000, right_ex=1000)
        else:
            cluster_pattern_alignment_list, cluster_bed_files_list = cluster_MSA_result

            # cluster_pattern_alignment_list have the same index with cluster_bed_files_list
            for i in range(len(cluster_pattern_alignment_list)):
                inner_cluster_MSA_result = clean_and_cluster_MSA(cluster_pattern_alignment_list[i],
                                                                 cluster_bed_files_list[i], MSA_dir,
                                                                 gap_threshold=0.8, clean_column_threshold=0.15)
                if inner_cluster_MSA_result is False:
                    continue
                elif inner_cluster_MSA_result is True:
                    find_boundary_and_crop(cluster_bed_files_list[i], genome_file, MSA_dir, pfam_dir=pfam_dir,
                                           left_ex=1000, right_ex=1000)
                else:
                    inner_cluster_pattern_alignment_list, inner_cluster_bed_files_list = inner_cluster_MSA_result

                    for j in range(len(inner_cluster_bed_files_list)):
                        find_boundary_and_crop(inner_cluster_bed_files_list[j], genome_file, MSA_dir, pfam_dir=pfam_dir,
                                               left_ex=1000, right_ex=1000)
    except Exception as e:
        print(f"Error during boundary finding and cropping for sequence: {seq_name}. Error: {str(e)}")
        traceback.print_exc()

    MSA_dir_representative = os.path.join(output_dir, "Multiple_sequence_alignment_less")
    file_copy_pattern = [
        r".*me_plot.pdf$",
        r".*gap_rm.fa_bou_crop.fa$",
        r".*_proof_anno_me.fa$"
    ]

    file_copier_object = FileCopier(MSA_dir, MSA_dir_representative, file_copy_pattern)
    file_copier_object.move_files()

    # After all processing is done, write the name of the file to the progress file
    progress_file = os.path.join(Single_files_dir, "progress_file.txt")
    with open(progress_file, "a") as f:
        f.write(seq_name + "\n")

def check_progress_file(progress_file_path):
    completed_sequences = []

    # Check if the progress file exists
    if os.path.isfile(progress_file_path):
        # If it does exist, open the file and read the list of completed sequences
        with open(progress_file_path, 'r') as progress_file:
            completed_sequences = progress_file.readlines()

    # Return the list of completed sequences (stripping newline characters)
    return [sequence.strip() for sequence in completed_sequences]

if_continue = False

if if_continue:
    # Check which sequences have already been processed
    completed_sequences = check_progress_file('progress.txt')

    # Filter out already processed sequences from the sequence names to process
    sequence_names = [f for f in os.listdir(Single_files_dir) if f.endswith(".fasta") and f not in completed_sequences]
else:
    sequence_names = [f for f in os.listdir(Single_files_dir) if f.endswith(".fasta")]



num_threads = 30
# Using a ProcessPoolExecutor to run the function in parallel
with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
    executor.map(analyze_sequence, sequence_names)







#TODO deal with single copy intact TE

#TODO linux command line tool development































