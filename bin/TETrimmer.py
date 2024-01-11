# Standard library imports
import os
import shutil
import traceback
from datetime import timedelta, datetime
import multiprocessing as mp
import click
import concurrent.futures
import json
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Local imports
import analyze
from functions import separate_sequences, remove_files_with_start_pattern, \
    change_permissions_recursive, repeatmasker, check_database, cd_hit_est, repeatmasker_output_classify, \
    rename_cons_file, rename_files_based_on_dict, prcyan, prgre, parse_cd_hit_est_result, \
    copy_files_with_start_pattern, fasta_file_to_dict, multi_seq_dotplot

#####################################################################################################
# Code block: Import JSON species_config file and define the default parameters
#####################################################################################################

# Load species-specific default values from the JSON config file
species_config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'species_config.json')
# Load the JSON configuration file
with open(species_config_path, "r") as config_file:
    species_config = json.load(config_file)

#####################################################################################################
# Code block: Main functions of TETrimmer
#####################################################################################################


@click.command(context_settings=dict(max_content_width=120),
               help="""\b
               ##########################################################################################
               \b
               ████████╗███████╗████████╗██████╗ ██╗███╗   ███╗███╗   ███╗███████╗██████╗ 
               ╚══██╔══╝██╔════╝╚══██╔══╝██╔══██╗██║████╗ ████║████╗ ████║██╔════╝██╔══██╗
                  ██║   █████╗     ██║   ██████╔╝██║██╔████╔██║██╔████╔██║█████╗  ██████╔╝
                  ██║   ██╔══╝     ██║   ██╔══██╗██║██║╚██╔╝██║██║╚██╔╝██║██╔══╝  ██╔══██╗
                  ██║   ███████╗   ██║   ██║  ██║██║██║ ╚═╝ ██║██║ ╚═╝ ██║███████╗██║  ██║
                  ╚═╝   ╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝╚═╝     ╚═╝╚═╝     ╚═╝╚══════╝╚═╝  ╚═╝
                  
                Version: v1.0.2 (19/Dec/2023) 

                Github: https://github.com/qjiangzhao/TE-Trimmer

                Developers:                                                                                                       
                Jiangzhao Qian;      RWTH Aachen University;                Email: jqian@bio1.rwth-aachen.de                          
                Hang Xue;            University of California, Berkeley;    Email: hang_xue@berkeley.edu                                
                Stefan Kusch;        Research Center Juelich;               Email: s.kusch@fz-juelich.de

                Funding source:                                                                                                        
                Ralph Panstruga Lab; RWTH Aachen University;                Email: panstruga@bio1.rwth-aachen.de            
                Website: https://www.bio1.rwth-aachen.de/PlantMolCellBiology/index.html                                                                            

                ##########################################################################################              

                python ./path_to_TETrimmer_bin/TETrimmer.py -i <TE_consensus_file> -g <genome_file>

                TETrimmer is designed to replace manual curation of transposable elements (TEs). 

                Two mandatory arguments are required, including 
                <genome file>, the genome FASTA file, and 
                <TE consensus file> from TE annotation software like RepeatModeler, EDTA, or REPET. 
                TETrimmer can do BLAST, sequence extension, multiple sequence alignment, and defining TE boundaries.

""")
@click.option('--input_file', '-i', required=True, type=str,
              help='Path to TE consensus file (FASTA format). Use the output from RepeatModeler, EDTA, REPET, et al.')
@click.option('--genome_file', '-g', required=True, type=str,
              help='Path to genome FASTA file.')
@click.option('--output_dir', '-o', default=os.getcwd(), type=str,
              help='Output directory. Default: current working directory.')
@click.option('--species', '-s', default='fungi', type=click.Choice(species_config.keys()),
              help='Select the type of organism for which to run TETrimmer.')
#@click.option('--engine', '-e', default='blast', type=click.Choice(["blast", "mmseqs"]),
#             help='Select the similar sequence search engine. "blast" or "mmseqs". Default: blast')
@click.option('--continue_analysis', '-ca', default=False, is_flag=True,
              help='Continue analysis after interruption.')
@click.option('--dedup', default=False, is_flag=True,
              help='Remove duplicate sequences in input file.')
@click.option('--genome_anno', '-ga', default=False, is_flag=True,
              help='Perform genome TE annotation using the TETrimmer curated database. Requires RepeatMasker.')
@click.option('--hmm', default=False, is_flag=True,
              help='Generate HMM files for each consensus sequence.')
@click.option('--debug', default=False, is_flag=True,
              help='Open debug mode. This will keep all raw files. WARNING: Many files will be generated.')
@click.option('--fast_mode', default=False, is_flag=True,
              help='Reduce running time at the cost of lower accuracy and specificity.')
#@click.option('--plot_query', default=False, is_flag=True,
#              help='Generate TE_Aid plot for each query sequence before TETrimmer analysis.')
#@click.option('--plot_skip', default=False, is_flag=True,
#              help='Generate TE_Aid plot for skipped elements.')
@click.option('--pfam_dir', '-pd', default=None, type=str,
              help='PFAM database directory. Omit this option if you do not have a local PFAM database '
                   'TETrimmer will download the database automatically in this case.')
@click.option('--cons_thr', type=float,
              help='Threshold used for the final consensus sequence generation. Default: 0.8')
@click.option('--mini_orf', type=int,
              help='Define the minimum ORF length to be predicted by TETrimmer. Default: 200')
@click.option('--max_msa_lines', type=int,
              help='Set the maximum sequence number for multiple sequence alignment. Default: 100')
@click.option('--top_msa_lines', type=int,
              help='When the sequence number of multiple sequence alignment (MSA) is greater than <--max_msa_lines>, '
                   'TETrimmer will sort sequences by length and choose <--top_msa_lines> number '
                   'of sequences. Then, TETrimmer will randomly select sequences from all remaining BLAST hits until '
                   '<--max_msa_lines> sequences are found for the multiple sequence alignment. Default: 70')
@click.option('--min_seq_num', type=int,
              help='The minimum sequence number for each multiple sequence alignment. Note: cannot be smaller than 10. '
                   'Default: 10')
@click.option('--min_blast_len', type=int,
              help='The minimum sequence length for BLAST hits. Default: 150')
@click.option('--max_cluster_num', default=2, type=int,
              help='The maximum cluster number for each multiple sequence alignment. Each multiple sequence alignment '
                   'can be divided into different clusters. TETrimmer will sort clusters by sequence number and choose '
                   'the top <--max_cluster_num> of clusters for further analysis. WARNING: A large number of clusters '
                   'will dramatically increase running time. Default: 2')
@click.option('--ext_thr', type=float,
              help="Sequence similarity threshold used for defining start and end of the consensus sequence, based on "
                   "the multiple sequence alignment. Nucleotides in each position with a similarity proportion smaller "
                   "than <--ext_thr> will be assigned the value N. If no N values are found, the algorithm will extend "
                   "the multiple sequence alignment to determine the limits of the consensus sequence. The lower the value "
                   "of <--ext_thr>, the longer the extensions on both ends. "
                   "Reduce <--ext_thr> if TETrimmer fails to determine the correct ends of repeat elements. Default: 0.7")
@click.option('--ext_check_win', type=str,
              help='Define check window size for extension. Default: 150')
@click.option('--ext_step', type=int,
              help='Number of nucleotides to be added to the left and right ends of the multiple sequence alignment. '
                   'TETrimmer will iteratively add <--ext_step> number of nucleotides until finding the boundary. '
                   'Default: 1000')
@click.option('--max_ext', type=int,
              help='The maximum extension in nucleotides at both ends of the multiple sequence alignment. Default: 7000')
@click.option('--gap_thr', type=float,
              help='If multiple sequence alignment positions (columns) have a gap proportion larger than <--gap_thr> '
                   'and the proportion of the most common nucleotide in this column is less than <--gap_nul_thr>, '
                   'this column will be removed from the consensus. Default: 0.4')
@click.option('--gap_nul_thr', type=float,
              help='Set nucleotide proportion threshold for keeping the column of the multiple sequence alignment. '
                   'Used with the <--gap_thr> option. Default: 0.7')
@click.option('--crop_end_div_thr', type=float,
              help='The crop end by divergence function will convert each nucleotide in the multiple sequence '
                   'alignment into a proportion value. This function will iteratively choose a sliding window from '
                   'each end of each sequence of the MSA and sum up the proportion numbers in this window. '
                   'The cropping will continue until the sum of proportions is larger than <--crop_end_div_thr>. '
                   'Cropped nucleotides will be converted to -. Default: 0.7')
@click.option('--crop_end_div_win', type=int,
              help='Window size used for the end-cropping process. Used with the <--crop_end_div_thr> option. Default: 40')
@click.option('--crop_end_gap_thr', type=float,
              help='The crop end by gap function will iteratively choose a sliding window from each end of each sequence '
                   'of the MSA and calculate the gap proportion in this window. The cropping will continue until the sum '
                   'of gap proportions is smaller than <--crop_end_gap_thr>. Cropped nucleotides will be converted to -. '
                   'Default: 0.1')
@click.option('--crop_end_gap_win', type=int,
              help='Define window size used to crop end by gap, used with the <--crop_end_gap_thr> option. Default: 250')
@click.option('--start_patterns', type=str,
              help='LTR elements always start with a conserved sequence pattern. TETrimmer searches the '
                   'beginning of the consensus sequence for these patterns. If the pattern is not found, '
                   'TETrimmer will extend the search of <--start_patterns> to up to 15 nucleotides from the '
                   'beginning of the consensus sequence and redefine the start of the consensus sequence '
                   'if the pattern is found. Note: The user can provide multiple LTR start patterns in a '
                   'comma-separated list, like: TG,TA,TC (no spaces; the order of patterns determines '
                   'the priority for the search). Default: TG')
@click.option('--end_patterns', type=str,
              help='LTR elements always end with a conserved sequence pattern. TETrimmer searches the '
                   'end of the consensus sequence for these patterns. If the pattern is not found, '
                   'TETrimmer will extend the search of <--end_patterns> to up to 15 nucleotides from the '
                   'end of the consensus sequence and redefine the end of the consensus sequence '
                   'if the pattern is found. Note: The user can provide multiple LTR end patterns in a '
                   'comma-separated list, like: CA,TA,GA (no spaces; the order of patterns determines '
                   'the priority for the search). Default: CA')
@click.option('--num_threads', '-t', default=10, type=int,
              help='Thread number used for TETrimmer. Default: 10')
@click.option('--classify_unknown', default=False, is_flag=True,
              help='Use RepeatClassifier to classify the consensus sequence if the input sequence is not classified or '
                   'is unknown or the processed sequence length by TETrimmer is 2000 bp longer or shorter '
                   'than the query sequence.')
@click.option('--classify_all', default=False, is_flag=True,
              help='Use RepeatClassifier to classify every consensus sequence. WARNING: This may take a long time.')
def main(input_file, genome_file, output_dir, continue_analysis, pfam_dir, min_blast_len, num_threads, max_msa_lines,
         top_msa_lines, min_seq_num, max_cluster_num, cons_thr, ext_thr, ext_step,
         max_ext, gap_thr, gap_nul_thr, crop_end_div_thr, crop_end_div_win, crop_end_gap_thr, crop_end_gap_win,
         start_patterns, end_patterns, mini_orf, species, ext_check_win, dedup, genome_anno, hmm,
         debug, fast_mode, classify_unknown, classify_all):

    # Add this to click options if mmseq2 has been fully tested
    engine = "blast"

    # Set plot_query and plot_skip to true
    plot_query = True
    plot_skip = True
    start_time = datetime.now()
    print(f"\nTETrimmer started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}.\n", flush=True)

    #####################################################################################################
    # Code block: File ends
    #####################################################################################################
    # _fm.bed: final multiple sequence alignment.
    # _u.bed: uniqueness.
    # f.bed: When blast hit more than 100, TETrimmer will only select the top 100 long for MSA
    # _bcl.fa bed clean. Use awk to remove letters that aren't A G C T a g c t in the fasta file
    # _bcln.fa bed clean and use name column as fasta header
    # _n.bed use bed file name column as fasta header
    # _n.fa use bed file name column as fasta header
    # _aln.fa alignment. File after alignment
    # _cl.fa clean. Convert nucleotide with very proportion in MSA into gap
    # _gs.fa gap similarity. Remove gaps in MSA with similarity check
    # _se.fa start end. Selected MSA based on given start and end position
    # _ce.fa crop end for both left and right sides. MSA after crop end by similarity
    # _cel.fa crop end left. Only crop left side of the MSA
    # _cer.fa crop end right. Only crop right side of the MSA
    # _bc.fa boundary crop.
    # _co.fa consensus. Consensus sequence file
    # _orf.txt ORF prediction file
    # _me.fa merged. Merged fasta file
    # _me.pdf merged. Merged pdf file
    # _proof_anno_me.fa proof annotation merged
    # .b blast
    # _orfm.txt ORF modified
    # _orfmt.txt ORF modified table
    # _pf.txt PFAM
    # _pfm.txt PFAM modified
    # _gr.fa gap remove. Remove sequences that contain too many gaps from MSA.
    # ps.pdf pdf file converted from ps format
    # su.pdf scale up pdf file
    # rc.fa reverse complement fasta file
    # cd.fa cd-hit-est output file

    #####################################################################################################
    # Code block: Change permissions of Aliview and TE_Aid
    #####################################################################################################

    # Change TE_Aid permission
    TE_aid_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "TE-Aid-master")
    # Change permissions of the directory and all its content to 755
    # 755 in octal corresponds to rwxr-xr-x
    change_permission = change_permissions_recursive(TE_aid_path, 0o755)
    if not change_permission:
        return

    #####################################################################################################
    # Code block: Define the default options according to the given species
    #####################################################################################################

    default_values = species_config.get(species, {})
    if cons_thr is None:
        cons_thr = default_values.get("cons_thr")

    if max_msa_lines is None:
        max_msa_lines = default_values.get("max_msa_lines")

    if top_msa_lines is None:
        top_msa_lines = default_values.get("top_msa_lines")

    if min_seq_num is None:
        min_seq_num = default_values.get("min_seq_num")
        if min_seq_num < 10:
            min_seq_num = 10

    if min_blast_len is None:
        min_blast_len = default_values.get("min_blast_len")

    if ext_thr is None:
        ext_thr = default_values.get("ext_thr")

    if ext_step is None:
        ext_step = default_values.get("ext_step")

    if max_ext is None:
        max_ext = default_values.get("max_ext")

    if gap_thr is None:
        gap_thr = default_values.get("gap_thr")

    if gap_nul_thr is None:
        gap_nul_thr = default_values.get("gap_nul_thr")

    if crop_end_div_thr is None:
        crop_end_div_thr = default_values.get("crop_end_div_thr")

    if crop_end_div_win is None:
        crop_end_div_win = default_values.get("crop_end_div_win")

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

    if ext_check_win is None:
        ext_check_win = default_values.get("ext_check_win")

    #####################################################################################################
    # Code block: Define input file, output directory, genome, check BLAST database
    #####################################################################################################
    try:
        bin_py_path, output_dir, single_file_dir, MSA_dir, classification_dir, hmm_dir, proof_annotation_dir, \
            low_copy_dir, perfect_proof, good_proof, intermediate_proof, need_check_proof, progress_file, pfam_dir, \
            final_con_file, final_con_file_no_low_copy, final_unknown_con_file, final_classified_con_file, \
            error_files, input_file, genome_file, skipped_dir, cluster_proof_anno_dir\
            = analyze.create_dir(continue_analysis, hmm, pfam_dir, output_dir, input_file, genome_file, plot_skip)
    except Exception:
        return

    #####################################################################################################
    # Code block: Copy TETrimmer_proof_anno_GUI to proof_annotation_dir
    #####################################################################################################
    proof_anno_GUI_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "TETrimmer_proof_anno_GUI")

    # Create the full path for the new directory inside the destination
    proof_anno_GUI_destination_dir = os.path.join(proof_annotation_dir, os.path.basename(proof_anno_GUI_dir))

    # Check if the destination directory exists
    if os.path.exists(proof_anno_GUI_destination_dir):
        # If it exists, remove it first
        shutil.rmtree(proof_anno_GUI_destination_dir)

    # Copy the entire directory
    shutil.copytree(proof_anno_GUI_dir, proof_anno_GUI_destination_dir)

    #####################################################################################################
    # Code block: Remove duplications in input file if required and generate single FASTA file
    #####################################################################################################

    # Generate single files when continue_analysis is false
    if not continue_analysis:
        # Do CD-HIT-EST merge if merge is true and continue_analysis is false
        if dedup:
            click.echo("\nTETrimmer is removing input sequences duplications, this might take some time.\n")
            merge_output = os.path.join(output_dir, f"{input_file}_cd_hit.fa")

            # Remove duplicates
            try:
                cd_hit_est(input_file, merge_output, identity_thr=0.9, aL=0.9, aS=0.9,
                           s=0.9, thread=num_threads)
                # Convert input_file to merged input_file
                input_file = merge_output
                click.echo("Merge finished.\n")
            except Exception as e:
                prcyan("TETrimmer cannot perform the de-duplication step by CD-HIT-EST and will use the "
                       "input sequences directly. This may cause a significantly longer running time but "
                       "will not affect the final result.")
                prgre("You can also run CD-HIT-EST separately to remove redundant sequences:\n"
                      "cd-hit-est -i <input_file> -o <output_file> -T <thread number> -c 0.9 "
                      "-aL 0.9 -aS 0.9 -s 0.9 -l 30")

        # Separate FASTA into single files; if FASTA headers contain "/", " " or ":" convert to "_"
        # Call this function to separate to single FASTA files and create objects from input file
        seq_list = separate_sequences(input_file, single_file_dir, continue_analysis=False)

        # Calculate the total sequence number 
        single_fasta_n = len(seq_list)
        click.echo(f"{single_fasta_n} sequences are detected from the input file")

        # Create new object to check BLAST database availability
        # Check if BLAST database and genome length files are available, otherwise create these in the
        # same directory of genome file
        database_result = check_database(genome_file, search_type=engine)

        # If database returns errors, stop the whole analysis
        if not database_result:
            return
        # Initial call to print 0% progress
        analyze.printProgressBar(0, single_fasta_n, prefix='Progress:', suffix='Complete', length=50)

    else:
        # Check if it can perform continue analysis
        if not os.listdir(single_file_dir):
            prgre("\nWARNING: TETrimmer cannot continue analysis. Please make sure the output directory is "
                  "the same as in the previous interrupted run.\n")
            return

        else:
            click.echo("\nTETrimmer will continue analysis based on previous results.\n")

            # Create seq_list, which contains sequence objects using the single FASTA files.
            seq_list = separate_sequences(input_file, single_file_dir, continue_analysis=True)
            single_fasta_n = len(seq_list)

            # Check which sequences have already been processed
            complete_sequences, skipped_count, low_copy_count, classified_pro = analyze.check_progress_file(
                progress_file)

            # Filter out already complete sequences
            seq_list = [seq for seq in seq_list if seq.name not in complete_sequences]
            click.echo(f"\n{single_fasta_n - len(seq_list)} sequences has been processed previously.\n")

    #####################################################################################################
    # Code block: Enable multiple threads
    #####################################################################################################

    analyze_sequence_params = [
        (seq, genome_file, MSA_dir, min_blast_len, min_seq_num, max_msa_lines,
         top_msa_lines, max_cluster_num, cons_thr, ext_thr, ext_step, classification_dir,
         max_ext, gap_thr, gap_nul_thr, crop_end_div_thr, crop_end_div_win, crop_end_gap_thr, crop_end_gap_win,
         start_patterns, end_patterns, output_dir, pfam_dir, mini_orf, single_fasta_n, hmm, hmm_dir,
         ext_check_win, debug, progress_file, classify_unknown, classify_all,
         final_con_file, final_con_file_no_low_copy, final_unknown_con_file, final_classified_con_file, low_copy_dir,
         fast_mode, error_files, plot_skip, skipped_dir, plot_query, engine, proof_annotation_dir
         ) for seq in seq_list]

    # Using a ProcessPoolExecutor to run the function in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        executor.map(analyze.analyze_sequence_helper, analyze_sequence_params)
    # multiprocessing
    # with mp.Pool(processes=10) as p:
    # p.starmap(analyze_sequence, analyze_sequence_params)

    #####################################################################################################
    # Code block: Check if all sequences are finished
    #####################################################################################################

    # At the end of the program, check if all sequences have been processed
    completed_sequence, skipped_count, low_copy_count, classified_pro = analyze.check_progress_file(progress_file)

    # Calculate the total count
    processed_count = len(completed_sequence)

    if processed_count == single_fasta_n:
        click.echo(f"\n\nAll sequences have been analysed!\n"
                   f"In the analysed sequences {skipped_count} are skipped. Note: not all skipped sequences can have "
                   f"TE Aid plot in the 'TETrimmer_for_proof_annotation' folder.\n"
                   f"In the analysed sequences {low_copy_count} are identified as low copy TE.")

    else:
        remaining = single_fasta_n - processed_count
        click.echo(f"\n\n{remaining} sequences have not been analysed.\n"
                   f"In the analysed sequences {skipped_count} are skipped. Note: not all skipped sequences can have "
                   f"TE Aid plot in the 'TETrimmer_for_proof_annotation' folder.\n"
                   f"In the analysed sequences {low_copy_count} are identified as low copy TE.\n")
        prgre("You might find the reasons why those sequences were not analysed from the 'error_file.txt' in the "
              "'Multiple_sequence_alignment' directory.")

    #####################################################################################################
    # Code block: Finish classifying unknown consensus sequences and write sequences to consensus file
    #####################################################################################################

    # Final RepeatMasker classification is not necessary, skip in case of errors
    try:
        if 0.3 <= classified_pro < 0.99:
            click.echo("\nTETrimmer is doing the final classification. It uses the classified TE to classify "
                       "Unknown elements.")
            temp_repeatmasker_dir = os.path.join(classification_dir, "temp_repeatmasker_classification")

            if os.path.exists(final_unknown_con_file) and os.path.exists(final_classified_con_file):
                os.makedirs(temp_repeatmasker_dir, exist_ok=True)
                classification_out = repeatmasker(final_unknown_con_file, final_classified_con_file,
                                                  temp_repeatmasker_dir,
                                                  thread=num_threads, classify=True)

                if classification_out:
                    repeatmasker_out = os.path.join(temp_repeatmasker_dir,
                                                    "temp_TETrimmer_unknown_consensus.fasta.out")
                    reclassified_dict = repeatmasker_output_classify(repeatmasker_out, progress_file,
                                                                     min_iden=60, min_len=80, min_cov=0.5)
                    if reclassified_dict:
                        click.echo(
                            f"\n{len(reclassified_dict)} TE elements were re-classified by the "
                            f"final classification module.")

                        # Update final consensus file
                        rename_cons_file(final_con_file, reclassified_dict)
                        rename_files_based_on_dict(proof_annotation_dir, reclassified_dict)
                        rename_files_based_on_dict(perfect_proof, reclassified_dict)
                        rename_files_based_on_dict(good_proof, reclassified_dict)
                        rename_files_based_on_dict(intermediate_proof, reclassified_dict)
                        rename_files_based_on_dict(need_check_proof, reclassified_dict)
                        rename_files_based_on_dict(low_copy_dir, reclassified_dict, seq_name=True)
                        if hmm:
                            rename_files_based_on_dict(hmm_dir, reclassified_dict)
                    else:
                        click.echo("0 TE elements were re-classified by the final classification module.")

            else:
                prcyan(f"\nThe final classification module failed.")
                prgre("\nThis does not affect the final TE consensus sequences "
                      "You can choose to ignore this error.\n")
    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"\nFinal RepeatMasker classification is wrong.\n")
            f.write(tb_content + '\n\n')
        prcyan(f"\nThe final classification module failed with error {e}")
        prgre("\nThis does not affect the final TE consensus sequences "
              "You can choose to ignore this error. For traceback content, please refer to 'error_file.txt' "
              "in the 'Multiple_sequence_alignment' directory.\n")

    try:
        # Delete the empty folder inside Classification_dir. For some reason, thr RepeatClassifier folder can't be
        # totally deleted within a short time. Delete all the empty folder here
        for foldername in os.listdir(classification_dir):
            folder_path = os.path.join(classification_dir, foldername)

            # Check if it is a directory
            if os.path.isdir(folder_path):
                # Check if the directory is empty
                if not os.listdir(folder_path):
                    os.rmdir(folder_path)
    except Exception:
        # This won't affect the final result, pass it when any error happens
        pass

    #####################################################################################################
    # Code block: Merge consensus_file to remove output duplications
    #####################################################################################################

    final_merge_success = True
    # Define merged file
    cd_hit_est_final_merged = os.path.join(output_dir, "TETrimmer_consensus_merged.fasta")

    try:
        click.echo("\nTETrimmer is removing sequence duplications.")

        # Do first round of CD-HIT-EST
        cd_hit_merge_output_round1 = os.path.join(classification_dir, "TETrimmer_consensus_merged_round1.fasta")
        cd_hit_merge_output_round1_clstr = f"{cd_hit_merge_output_round1}.clstr"

        # Round 1 merge only requires that the alignment coverage for the shorter sequence is greater than 0.9
        # and the similarity is greater than 0.9
        cd_hit_est(final_con_file, cd_hit_merge_output_round1, identity_thr=0.9, aL=0, aS=0.9, s=0, thread=num_threads)

        # Read progress file
        progress_df = pd.read_csv(progress_file)

        # Create a dictionary with sequence names as keys
        sequence_info = {}
        for index, row in progress_df.iterrows():
            sequence_name = row["consensus_name"]
            evaluation = row["evaluation"] if pd.notna(row["evaluation"]) else "Unknown"  # Default value for NaN
            te_type = row["reclassified_type"] if pd.notna(row["reclassified_type"]) else "Unknown"
            length = row["cons_length"] if pd.notna(row["cons_length"]) else 0  # Default value for NaN
            sequence_info[sequence_name] = {"evaluation": evaluation, "type": te_type, "length": length}

        # Parse cd-hit-est result, clusters is a dictionary, the key is cluster number, value is a list
        # contain all sequence names in this cluster
        clusters, detailed_clusters = parse_cd_hit_est_result(cd_hit_merge_output_round1_clstr)

        # Check if sequences scored "Perfect" and "Good" are included in the cluster and choose the longest sequence
        best_sequences = []  # Define list to store "Perfect" or "Good" sequence names

        # Define list to store sequence in clusters that do not contain "Perfect" or "Good" sequences
        sequence_for_round2 = []
        for cluster_name, sequences in clusters.items():
            perfect_sequences = []
            good_sequences = []
            best_seq = None

            for seq in sequences:
                if seq in sequence_info:
                    evaluation = sequence_info[seq]["evaluation"]
                    length = sequence_info[seq]["length"]
                    if evaluation == "Perfect":
                        perfect_sequences.append((seq, length))
                    elif evaluation == "Good":
                        good_sequences.append((seq, length))
            # Choose the longest "Perfect" sequence, if have Perfect
            if perfect_sequences:
                best_seq = max(perfect_sequences, key=lambda x: x[1])[0]
            # If no "Perfect", choose the longest "Good" sequence
            elif good_sequences:
                best_seq = max(good_sequences, key=lambda x: x[1])[0]

            if best_seq:
                best_sequences.append(best_seq)
            else:
                sequence_for_round2.extend(clusters[cluster_name])  # extend() will create a flat list

        # Read the original consensus file
        consensus_sequences = SeqIO.parse(final_con_file, "fasta")

        # Define temporary file to store "Perfect" and "Good" sequences
        temp_consensus_round1 = os.path.join(classification_dir, "temp_consensus_round1.fasta")

        # Define temporary file to store remaining sequences for second round of CD-HIT-EST
        temp_consensus_round2_input = os.path.join(classification_dir, "temp_consensus_round2_input.fasta")

        # Write sequences to files
        with open(temp_consensus_round1, "w") as high_quality_file, \
                open(temp_consensus_round2_input, 'w') as round2_file:

            for seq_record in consensus_sequences:
                # Sequence names in best_sequences and sequence_for_round2 do not contain classification
                seq_id = seq_record.id.split("#")[0]

                if seq_id in best_sequences:
                    SeqIO.write(seq_record, high_quality_file, "fasta")
                elif seq_id in sequence_for_round2:
                    SeqIO.write(seq_record, round2_file, 'fasta')

        # Do second round of CD-HIT-EST based on temp_consensus_round2_input
        cd_hit_merge_output_round2 = os.path.join(classification_dir, "TETrimmer_consensus_merged_round2.fasta")

        # Round 2 merge requires that the alignment coverage for the long and short sequence are both greater than 0.8
        # and the similarity is greater than 0.85
        cd_hit_est(temp_consensus_round2_input, cd_hit_merge_output_round2, identity_thr=0.85, aL=0.8, aS=0.8, s=0.8,
                   thread=num_threads)

        # Combine the two files into a merged file
        with open(temp_consensus_round1, 'r') as file1, \
                open(cd_hit_merge_output_round2, 'r') as file2, \
                open(cd_hit_est_final_merged, 'w') as combined_file:
            # Write contents of the first file
            for line in file1:
                combined_file.write(line)

            # Write contents of the second file
            for line in file2:
                combined_file.write(line)

        # Find sequence names that are not included inside in cd_hit_est_final_merged file
        # Parse the sequences in the original and merged files
        original_sequences = SeqIO.parse(final_con_file, "fasta")
        merged_sequences = SeqIO.parse(cd_hit_est_final_merged, "fasta")

        # Extract sequence IDs from both files and store to set
        original_ids = {seq_record.id.split("#")[0] for seq_record in original_sequences}
        merged_ids = {seq_record.id.split("#")[0] for seq_record in merged_sequences}

        # Find the difference between the two sets to identify sequence names not included in the merged file
        missing_ids = original_ids - merged_ids

        """
        # Based on missing_ids delete files in proof annotation folder and HMM folder
        for missing_id in missing_ids:

            # if not, set evaluation_level to "Need_check". The "get" method will return the default value
            # when the key does not exist.
            evaluation_level = sequence_info.get(missing_id, {"evaluation": "Need_check"})["evaluation"]

            # Add '#' to the end of missing_id, this can avoid to delete 140 when the id is 14
            missing_id = f"{missing_id}#"
            if evaluation_level == "Perfect":
                remove_files_with_start_pattern(perfect_proof, missing_id, if_seq_name=False)
            elif evaluation_level == "Good":
                remove_files_with_start_pattern(good_proof, missing_id, if_seq_name=False)
            elif evaluation_level == "Reco_check":
                remove_files_with_start_pattern(intermediate_proof, missing_id, if_seq_name=False)
            elif evaluation_level == "Need_check":
                remove_files_with_start_pattern(need_check_proof, missing_id, if_seq_name=False)
            else:
                remove_files_with_start_pattern(low_copy_dir, missing_id, if_seq_name=False)

            if hmm:
                remove_files_with_start_pattern(hmm_dir, missing_ids)
        """
        click.echo(f"\nFinished to remove sequence duplications.\n")

    except Exception as e:
        final_merge_success = False
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"\nFinal CD-HIT-EST deduplication error.\n")
            f.write(tb_content + '\n')
        prcyan("\nThe final CD-HIT-EST merge step cannot be performed. Final TE consensus library redundancy can "
               "be higher but the sensitivity is not affected. You can remove duplicated sequence by yourself.")
        prgre("\nYou can choose to ignore CD-HIT-EST error. For traceback output, please refer to 'error_file.txt' "
              "in the 'Multiple_sequence_alignment' directory.\n")

    #####################################################################################################
    # Code block: Cluster proof annotation files
    #####################################################################################################

    try:
        click.echo("\nTETrimmer is clustering proof annotation files.\n")

        # Create directory to store multiple sequence dotplot
        multi_dotplot_dir = os.path.join(classification_dir, "Multiple_sequence_dotplot")
        os.makedirs(multi_dotplot_dir, exist_ok=True)

        # Load fast file to a dictionary, key is record.id, value is record project
        # When separate_name is true, the key of the dictionary will be the sequence name separated by '#'
        final_con_file_no_low_copy_dict = fasta_file_to_dict(final_con_file_no_low_copy, separate_name=True)

        # Clean cluster_proof_anno_dir when --continue_analysis is on.
        if continue_analysis:
            # When the start pattern isn't given, all files inside the folder will be removed
            remove_files_with_start_pattern(cluster_proof_anno_dir)
            remove_files_with_start_pattern(multi_dotplot_dir)

        # Do CD-HIT-EST for final consensus file without low copy elements
        final_con_file_no_low_copy_cd_out = f"{final_con_file_no_low_copy}_cd.fa"
        final_con_file_no_low_copy_clstr = f"{final_con_file_no_low_copy_cd_out}.clstr"

        # Round 1 merge only requires that the alignment coverage for the shorter sequence is greater than 0.9
        # and the similarity is greater than 0.9
        cd_hit_est(final_con_file_no_low_copy, final_con_file_no_low_copy_cd_out,
                   identity_thr=0.9, aL=0, aS=0.9, s=0, thread=num_threads)
        clusters_proof_anno, detailed_clusters_proof_anno = parse_cd_hit_est_result(
            final_con_file_no_low_copy_clstr)

        for cluster_name_proof_anno, seq_info_proof_anno in detailed_clusters_proof_anno.items():
            # Create cluster folder
            cluster_folder = os.path.join(cluster_proof_anno_dir, cluster_name_proof_anno)
            os.makedirs(cluster_folder, exist_ok=True)

            seq_info_proof_anno_len = len(seq_info_proof_anno)

            cluster_record_list = []
            for i in range(seq_info_proof_anno_len):
                try:
                    seq_length_proof_anno = seq_info_proof_anno[i][0]
                except Exception:
                    seq_length_proof_anno = None

                seq_name_proof_anno = seq_info_proof_anno[i][1]

                try:
                    seq_per_proof_anno = seq_info_proof_anno[i][2]
                except Exception:
                    seq_per_proof_anno = None
                try:
                    seq_direction_proof_anno = seq_info_proof_anno[i][3]
                except Exception:
                    seq_direction_proof_anno = None

                # Copy sequence files into cluster folder
                # if not, set evaluation_level to "Need_check". The "get" method will return the default value
                # when the key does not exist.
                evaluation_level = sequence_info.get(seq_name_proof_anno, {"evaluation": "Need_check"})["evaluation"]

                # Add '#' to the end of seq_name_proof_anno, this can avoid to delete 140 when the id is 14
                seq_name_proof_anno_m = f"{seq_name_proof_anno}#"
                if evaluation_level == "Perfect":
                    copy_files_with_start_pattern(perfect_proof, seq_name_proof_anno_m, cluster_folder,
                                                  seq_length_proof_anno, seq_per_proof_anno, evaluation_level)
                elif evaluation_level == "Good":
                    copy_files_with_start_pattern(good_proof, seq_name_proof_anno_m, cluster_folder,
                                                  seq_length_proof_anno, seq_per_proof_anno, evaluation_level)
                elif evaluation_level == "Reco_check":
                    copy_files_with_start_pattern(intermediate_proof, seq_name_proof_anno_m, cluster_folder,
                                                  seq_length_proof_anno, seq_per_proof_anno, evaluation_level)
                elif evaluation_level == "Need_check":
                    copy_files_with_start_pattern(need_check_proof, seq_name_proof_anno_m, cluster_folder,
                                                  seq_length_proof_anno, seq_per_proof_anno, evaluation_level)

                # Plot multiple sequence dotplot when more than one sequence are included inside one cluster
                if seq_info_proof_anno_len > 1:

                    # Extract record from dictionary
                    cluster_record = final_con_file_no_low_copy_dict.get(seq_name_proof_anno)

                    # When the sequence direction is negative, reverse complement it
                    if cluster_record is not None and seq_direction_proof_anno == '-':
                        rev_comp_cluster_record_seq = cluster_record.seq.reverse_complement()
                        cluster_record = SeqRecord(rev_comp_cluster_record_seq,
                                                   id=cluster_record.id, description="")
                    if cluster_record is not None:
                        cluster_record_list.append(cluster_record)

            if len(cluster_record_list) > 1:
                # Define and write cluster fasta file a
                cluster_fasta = os.path.join(multi_dotplot_dir, f"{cluster_name_proof_anno}.fa")
                SeqIO.write(cluster_record_list, cluster_fasta, "fasta")

                # Do multiple sequence dotplot
                multi_dotplot_pdf = multi_seq_dotplot(cluster_fasta, multi_dotplot_dir, cluster_name_proof_anno)

                # Move muti_dotplot_pdf to proof annotation cluster folder
                if os.path.isfile(multi_dotplot_pdf):
                    shutil.copy(multi_dotplot_pdf, cluster_folder)

        # clear remove_files_with_start_pattern folder
        if not debug and os.path.exists(multi_dotplot_dir):
            shutil.rmtree(multi_dotplot_dir)

    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"\nFinal clustering of proof annotation files failed.\n")
            f.write(tb_content + '\n\n')
        prcyan(f"\nFinal clustering of proof annotation files failed with error {e}")
        prcyan('\n' + tb_content + '')
        prgre("\nThis does not affect the final TE consensus sequences. But this can heavily complicate the "
              "TE proof annotation. If you don't plan to do proof annotation, you can choose to ignore "
              "this error.\n")

    #####################################################################################################
    # Code block: Whole-genome TE annotation
    #####################################################################################################

    try:
        # If 90% of the query sequences have been processed, RepeatMasker is allowed to perform whole genome annotation
        # if processed_count >= single_fasta_n * 0.9:

        # Run RepeatMasker
        if genome_anno:
            click.echo("\nTETrimmer is performing whole-genome TE annotation by RepeatMasker. This could take "
                       "long time. \nThe final TE consensus library has been completed. You can use it now.\n")

            if final_merge_success and os.path.exists(cd_hit_est_final_merged):
                repeatmakser_lib = cd_hit_est_final_merged
            else:
                repeatmakser_lib = final_con_file

            # make a new folder for RepeatMasker output
            repeatmasker_dir = os.path.join(output_dir, "RepeatMasker_result")
            if not os.path.exists(repeatmasker_dir):
                os.mkdir(repeatmasker_dir)
            genome_anno_result = \
                repeatmasker(genome_file, repeatmakser_lib, repeatmasker_dir, thread=num_threads)
            if genome_anno_result:
                click.echo("\nFinished whole genome TE annotation by RepeatMasker\n")

    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"\nGenome TE annotation error.\n")
            f.write(tb_content + '\n\n')

    end_time = datetime.now()
    duration = end_time - start_time

    # Remove microseconds from the duration
    duration_without_microseconds = timedelta(days=duration.days, seconds=duration.seconds)

    # if not debug:
    # Remove all single files after all sequences have been processed
    # shutil.rmtree(single_file_dir)

    analyze.printProgressBar(processed_count, single_fasta_n, prefix='Progress:', suffix='Complete', length=50,
                             final=True)
    print(f"\nTETrimmer analysis finished at {start_time.strftime('%Y-%m-%d %H:%M:%S')}.\n")
    print(f"TETrimmer runtime was {duration_without_microseconds}.")


# The following is necessary to make the script executable, i.e., python myscript.py.
if __name__ == '__main__':
    main()
