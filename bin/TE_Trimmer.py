# Standard library imports
import os
import traceback
from datetime import timedelta, datetime
import multiprocessing as mp
import click
import concurrent.futures
import json
import pandas as pd
from Bio import SeqIO

# Local imports
import analyze
from Function_blast_extension_mafft import separate_sequences, remove_files_with_start_pattern, \
    change_permissions_recursive, repeatmasker, check_database, cd_hit_est, repeatmasker_output_classify, \
    rename_cons_file, rename_files_based_on_dict

#####################################################################################################
# Code block: Import json species_config file and define the default parameters
#####################################################################################################

# Load species-specific default values from the JSON config
species_config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'species_config.json')
# Load the JSON configuration file
with open(species_config_path, "r") as config_file:
    species_config = json.load(config_file)

#####################################################################################################
# Code block: Mian function of TE Trimmer
#####################################################################################################


@click.command(context_settings=dict(max_content_width=120),
               help="""\b
               ##########################################################################################
               \b
               ████████╗███████╗    ████████╗██████╗ ██╗███╗   ███╗███╗   ███╗███████╗██████╗ 
               ╚══██╔══╝██╔════╝    ╚══██╔══╝██╔══██╗██║████╗ ████║████╗ ████║██╔════╝██╔══██╗
                  ██║   █████╗         ██║   ██████╔╝██║██╔████╔██║██╔████╔██║█████╗  ██████╔╝
                  ██║   ██╔══╝         ██║   ██╔══██╗██║██║╚██╔╝██║██║╚██╔╝██║██╔══╝  ██╔══██╗
                  ██║   ███████╗       ██║   ██║  ██║██║██║ ╚═╝ ██║██║ ╚═╝ ██║███████╗██║  ██║
                  ╚═╝   ╚══════╝       ╚═╝   ╚═╝  ╚═╝╚═╝╚═╝     ╚═╝╚═╝     ╚═╝╚══════╝╚═╝  ╚═╝
                  
                Version: v1.1 (27/OCT/2023) 

                Github: https://github.com/qjiangzhao/TE-Trimmer

                Developers:                                                                                                       
                Jiangzhao Qian;  RWTH Aachen University;                Email: jqian@bio1.rwth-aachen.de                          
                Hang Xue;        University of California, Berkeley;    Email: hang_xue@berkeley.edu

                Funding source:                                                                                         
                Panstruga's Lab. Website: https://www.bio1.rwth-aachen.de/PlantMolCellBiology/index.html                 

                Many thanks to Dr. Stefan Kusch                                                           

                ##########################################################################################              

                python ./path_to_TE_Trimmer_bin/TE_Trimmer.py -i <TE_consensus_file> -o <genome_file>

                TE Trimmer is designed to replace transposable element (TE) manual curation. 

                Two mandatory arguments are required including <genome file> and <TE consensus file> from TE 
                annotation software like RepeatModeler or EDTA et al. TE Trimmer can do blast, extension, multiple sequence alignment, and defining TE boundaries.

""")
@click.option('--input_file', '-i', required=True, type=str,
              help='Path to TE consensus file (FASTA format). Use the output from RepeatModeler or EDTA et al.')
@click.option('--genome_file', '-g', required=True, type=str,
              help='Path to genome FASTA file.')
@click.option('--output_dir', '-o', default=os.getcwd(), type=str,
              help='Output directory. Default: current directory.')
@click.option('--species', '-s', default='fungi', type=click.Choice(species_config.keys()),
              help='Select the species for which you want to run TE Trimmer.')
@click.option('--continue_analysis', '-ca', default=False, is_flag=True,
              help='Continue to analysis after interruption.')
@click.option('--dedup', default=False, is_flag=True,
              help='Remove duplicate sequences in input file.')
@click.option('--genome_anno', default=False, is_flag=True,
              help='Perform genome TE annotation using the TE Trimmer curated database. Requires RepeatMasker.')
@click.option('--hmm', default=False, is_flag=True,
              help='Generate HMM files for each consensus sequences.')
@click.option('--debug', default=False, is_flag=True,
              help='Open debug mode. This will keep all raw files. WARNING: Many files will be produced.')
@click.option('--fast_mode', default=False, is_flag=True,
              help='Reduce running time but at the cost of lower accuracy and specificity.')
@click.option('--plot_query', default=False, is_flag=True,
              help='Perform TE_Aid plot for each query sequences before TE Trimmer analysis.')
@click.option('--plot_skip', default=False, is_flag=True,
              help='Perform TE_Aid plot for skipped elements')
@click.option('--pfam_dir', default=None, type=str,
              help='Pfam database directory. Omit this option if you do not have a local PFAM database. '
                   'TE Trimmer will download the database automatically in this case.')
@click.option('--cons_thr', type=float,
              help='Threshold used for the final consensus sequence generation. Default: 0.8')
@click.option('--mini_orf', type=int,
              help='Define the minimum ORF length that will be predicted by TE Trimmer. Default: 200')
@click.option('--max_msa_lines', type=int,
              help='Set the maximum sequences number for multiple sequence alignment. Default: 100')
@click.option('--top_mas_lines', type=int,
              help='When the sequence number of multiple sequence alignment (MSA) is greater than <--max_msa_lines>, '
                   'TE Trimmer will sort sequences by length and choose <--top_msa_lines> number '
                   'of sequences. Then, TE Trimmer will randomly select sequences from all remaining BLAST hits until '
                   '<--max_msa_lines> sequences are found for the multiple sequence alignment. Default: 70')
@click.option('--min_seq_num', type=int,
              help='The minimum sequence number for each multiple sequence alignment. Note: can not smaller than 10. '
                   'Default: 10')
@click.option('--min_blast_len', type=int,
              help='The minimum sequence length for blast hits. Default: 150')
@click.option('--max_cluster_num', type=int,
              help='The maximum cluster number for each multiple sequence alignment. Each multiple '
                   'sequence alignment can be divided into different clusters. TE Trimmer will sort '
                   'clusters by sequence number and choose the top <--max_cluster_num> of clusters for '
                   'further analysis. WARNING: Big number will dramatically increase running time. Default: 2')
@click.option('--ext_thr', type=float,
              help="Threshold used for define the extension extent. The lower the value of <--ext_thr>, the easier the "
                   "extensions on both ends be longer. Reduce <--ext_thr> if TE Trimmer fails to determine the correct "
                   "ends of repeat elements. Default: 0.7")
@click.option('--ext_check_win', type=str,
              help='Define check windows size for extension. Default: 150')
@click.option('--ext_step', type=int,
              help='Number of nucleotides to be added to the left and right ends of the multiple sequence alignment. '
                   'TE_Trimmer will iteratively add <--ext_step> number of nucleotides until finding the boundary. '
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
                   'Cropped nucleotides will be converted to -. Default: 0.8')
@click.option('--crop_end_div_win', type=int,
              help='Window size used for the end-cropping process. Used with --crop_end_div_thr option. Default: 20')
@click.option('--crop_end_gap_thr', type=float,
              help='The crop end by gap function will iteratively choose a sliding window from '
                   'each end of each sequence of the MSA and calculate the gap proportion in this window. '
                   'The cropping will continue until the sum of gap proportions is smaller than <--crop_end_gap_thr>. '
                   'Cropped nucleotides will be converted to -. Default: 0.1')
@click.option('--crop_end_gap_win', type=int,
              help='Define window size used to crop end by gap, used with <--crop_end_gap_thr> option. Default: 250')
@click.option('--start_patterns', type=str,
              help='LTR elements always start with a conserved sequence pattern. TE Trimmer searches the '
                   'beginning of the consensus sequence for these patterns. If the pattern is not found, '
                   'it will extend the search of <--start_patterns> to up to 15 nucleotides from the beginning '
                   'of the consensus sequence and redefine the start of the consensus sequence '
                   'if the pattern is found. Note: The user can provide multiple LTR start patterns in a '
                   'comma-separated list, like: TG,TA,TC (no spaces; the order of patterns determines '
                   'the priority for the search). Default: TG')
@click.option('--end_patterns', type=str,
              help='LTR elements always end with a conserved sequence pattern. TE Trimmer searches the '
                   'end of the consensus sequence for these patterns. If the pattern is not found, '
                   'it will extend the search of <--end_patterns> to up to 15 nucleotides from the end '
                   'of the consensus sequence and redefine the end of the consensus sequence '
                   'if the pattern is found. Note: The user can provide multiple LTR end patterns in a '
                   'comma-separated list, like: CA,TA,GA (no spaces; the order of patterns determines '
                   'the priority for the search). Default: CA')
@click.option('--num_threads', '-t', default=10, type=int,
              help='Threads numbers used for TE Trimmer. Default: 10')
@click.option('--classify_unknown', default=False, is_flag=True,
              help='Use RepeatClassifier to classify the consensus sequence if the input sequence is not classified or '
                   'is unknown.')
@click.option('--classify_all', default=False, is_flag=True,
              help='Use RepeatClassifier to classify every consensus sequence.  WARNING: it will take longer time.')
def main(input_file, genome_file, output_dir, continue_analysis, pfam_dir, min_blast_len, num_threads, max_msa_lines,
         top_mas_lines, min_seq_num, max_cluster_num, cons_thr, ext_thr, ext_step, plot_query, plot_skip,
         max_ext, gap_thr, gap_nul_thr, crop_end_div_thr, crop_end_div_win, crop_end_gap_thr, crop_end_gap_win,
         start_patterns, end_patterns, mini_orf, species, ext_check_win, dedup, genome_anno, hmm,
         debug, fast_mode, classify_unknown, classify_all):
    start_time = datetime.now()
    print(f"\nTE Trimmer started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n", flush=True)

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

    # Change Aliview permission
    aliview_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "aliview")

    # Change permissions of the directory and all its content to 755
    # 755 in octal corresponds to rwxr-xr-x
    change_permissions_recursive(aliview_path, 0o755)

    #####################################################################################################
    # Code block: Define the default options according to the given species
    #####################################################################################################

    default_values = species_config.get(species, {})
    if cons_thr is None:
        cons_thr = default_values.get("cons_thr")

    if max_msa_lines is None:
        max_msa_lines = default_values.get("max_msa_lines")

    if top_mas_lines is None:
        top_mas_lines = default_values.get("top_mas_lines")

    if min_seq_num is None:
        min_seq_num = default_values.get("min_seq_num")
        if min_seq_num < 10:
            min_seq_num = 10

    if min_blast_len is None:
        min_blast_len = default_values.get("min_blast_len")

    if max_cluster_num is None:
        max_cluster_num = default_values.get("max_cluster_num")

        # Convert string "False" to boolean
        if max_cluster_num == "False":
            max_cluster_num = False

        if fast_mode:
            max_cluster_num = 2

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
    # Code block: Define input file, output directory, genome, check blast database
    #####################################################################################################

    bin_py_path, output_dir, single_file_dir, MSA_dir, classification_dir, hmm_dir, proof_annotation_dir, \
        low_copy_dir, perfect_proof, good_proof, intermediate_proof, need_check_proof, progress_file, pfam_dir, \
        final_con_file, final_unknown_con_file, final_classified_con_file, error_files, input_file, genome_file, \
        skipped_dir = analyze.create_dir(continue_analysis, hmm, pfam_dir, output_dir, input_file, genome_file,
                                         plot_skip)

    #####################################################################################################
    # Code block: Remove duplications in input file when it is required and generate single fasta file
    #####################################################################################################

    # Generate single files when continue_analysis is false
    if not continue_analysis:
        # Do cd-hit-est merge when merge is true and continue_analysis is false
        if dedup:
            click.echo("\nTE Trimmer is merging input sequences, this might take some time.\n")
            merge_output = os.path.join(output_dir, f"{input_file}_cd_hit.fa")

            # Set lower identity threshold for the query, this can increase sensitive
            # Merge step will remove single LTR but nested TEs can mask other TEs
            cd_iden_thr = 0.9
            cd_alin_s = 0.85

            if fast_mode:
                cd_iden_thr = 0.8
                cd_alin_s = 0.8

            cd_hit_est(input_file, merge_output, identity_thr=cd_iden_thr, aL=0, aS=cd_alin_s, s=0, thread=num_threads)

            # Convert input_file to merged input_file
            input_file = merge_output
            click.echo("Merge finished.\n")

        # Separate fasta to single files, if fasta header contain "/" or " " or ":" convert them to "_"
        # Call this function to separate to single fasta files and create objects from input file
        seq_list = separate_sequences(input_file, single_file_dir, continue_analysis=False)

        # Calculate the total sequence number 
        single_fasta_n = len(seq_list)
        click.echo(f"{single_fasta_n} sequences are detected from the input file")

        # Create new object to check blast database availability
        # Check if blast database and genome length files are available, otherwise create them at the
        # same directory of genome file
        check_database(genome_file)
        # Initial call to print 0% progress
        analyze.printProgressBar(0, single_fasta_n, prefix='Progress:', suffix='Complete', length=50)

    else:
        # Check if the can perform continue analysis
        if not os.listdir(single_file_dir):
            click.echo("\nWARNING: TE Trimmer can't do continue analysis, please make sure the output directory is same"
                       " with your previous analysis.\n")
            return

        else:
            click.echo("\nTE Trimmer will continue to analyze based on previous results.\n")

            # Create seq_list, which contain sequence objects using the single fasta files.
            seq_list = separate_sequences(input_file, single_file_dir, continue_analysis=True)
            single_fasta_n = len(seq_list)

            # Check which sequences have already been processed
            complete_sequences, skipped_count, low_copy_count, classified_pro = analyze.check_progress_file(
                progress_file)

            # Filter out already complete sequences from the total sequences
            seq_list = [seq for seq in seq_list if seq.name not in complete_sequences]
            click.echo(f"\n{single_fasta_n - len(seq_list)} sequences has been processed previously.\n")

    #####################################################################################################
    # Code block: Enable multiple threads
    #####################################################################################################

    analyze_sequence_params = [
        (seq, genome_file, MSA_dir, min_blast_len, min_seq_num, max_msa_lines,
         top_mas_lines, max_cluster_num, cons_thr, ext_thr, ext_step, classification_dir,
         max_ext, gap_thr, gap_nul_thr, crop_end_div_thr, crop_end_div_win, crop_end_gap_thr, crop_end_gap_win,
         start_patterns, end_patterns, output_dir, pfam_dir, mini_orf, single_fasta_n, hmm,
         ext_check_win, debug, progress_file, classify_unknown, classify_all,
         final_con_file, final_unknown_con_file, final_classified_con_file, low_copy_dir, fast_mode, error_files,
         plot_skip, skipped_dir, plot_query
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
        click.echo(f"All sequences have been analysed!\n"
                   f"In the analysed sequences {skipped_count} are skipped.\n"
                   f"In the analysed sequences {low_copy_count} are identified as low copy TE\n"
                   f"TE Trimmer is doing the final classification, merging, or whole genome annotation step.")

    else:
        remaining = single_fasta_n - processed_count
        click.echo(f"\n\n{remaining} sequences have not been analysed.\n")

    #####################################################################################################
    # Code block: Finish classifying unknown consensus file and writing sequences back to consensus file
    #####################################################################################################

    # Suppress RepeatMasker final classification under fast_mode
    if fast_mode:
        classified_pro = 0.01

    # Final RepeatMasker classification isn't necessary, skip it when errors are there
    try:
        if 0.3 <= classified_pro < 0.99:
            temp_repeatmasker_dir = os.path.join(classification_dir, "temp_repeatmasker_classification")

            if os.path.exists(final_unknown_con_file) and os.path.exists(final_classified_con_file):
                os.makedirs(temp_repeatmasker_dir, exist_ok=True)
                classification_out = repeatmasker(final_unknown_con_file, final_classified_con_file,
                                                  temp_repeatmasker_dir,
                                                  thread=num_threads, classify=True)

                if classification_out:
                    repeatmasker_out = os.path.join(temp_repeatmasker_dir,
                                                    "temp_TE_Trimmer_unknown_consensus.fasta.out")
                    reclassified_dict = repeatmasker_output_classify(repeatmasker_out, progress_file,
                                                                     min_iden=60, min_len=80, min_cov=0.5)
                    if reclassified_dict:
                        click.echo(
                            f"{len(reclassified_dict)} TE elements were re-classified by final classification module")

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
                click.echo("One of the consensus file does not exist, RepeatMasker reclassify does not run, "
                           "This will not affect the final consensus sequences.")
    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"\nFinal RepeatMasker classification is wrong.\n")
            f.write(tb_content + '\n\n')

        click.echo("Note: The final TE classification is skipped, this won't affect the TE consensus result.")

    #####################################################################################################
    # Code block: merge consensus_file to remove duplications
    #####################################################################################################

    try:
        click.echo("TE Trimmer is removing sequence duplications")

        # Do first round cd-hit-est
        cd_hit_merge_output_round1 = os.path.join(classification_dir, "TE_Trimmer_consensus_merged_round1.fasta")
        cd_hit_merge_output_round1_clstr = f"{cd_hit_merge_output_round1}.clstr"

        # Round 1 merge only require that the alignment coverage for the shorter sequence is greater than 0.8
        # and the similarity is greater than 0.85
        cd_hit_est(final_con_file, cd_hit_merge_output_round1, identity_thr=0.85, aS=0.8, s=0, thread=num_threads)

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

        # Parse cd-hit-est result
        clusters = {}
        current_cluster = []
        with open(cd_hit_merge_output_round1_clstr, "r") as f:
            for line in f:
                if line.startswith(">Cluster"):
                    cluster_name = line.strip().replace(" ", "")  # Remove the empty space in the cluster name
                    if current_cluster:  # If the current_cluster isn't empty
                        clusters[cluster_name] = current_cluster
                        current_cluster = []
                else:
                    seq_info = line.split(">")[1].split("...")[0].split("#")[0]
                    current_cluster.append(seq_info)
            if current_cluster:
                clusters[cluster_name] = current_cluster

        # Check if "Perfect" and "Good" level sequences are included inside cluster and choose the longest one
        best_sequences = []  # Define list to store Perfect or Good sequence names

        # Define list to store sequence in clusters that don't contain Perfect or Good sequences
        sequence_for_round2 = []
        for cluster_name, sequences in clusters.items():
            best_seq = None
            best_length = 0
            for seq in sequences:
                if seq in sequence_info:
                    evaluation = sequence_info[seq]["evaluation"]
                    length = sequence_info[seq]["length"]
                    if evaluation in ["Perfect", "Good"] and length > best_length:
                        best_length = length
                        best_seq = seq
            if best_seq:
                best_sequences.append(best_seq)
            else:
                sequence_for_round2.extend(clusters[cluster_name])  # extend() will make a flat a list

        # Read the original consensus file
        consensus_sequences = SeqIO.parse(final_con_file, "fasta")

        # Define temporary file to store Perfect or Good sequences
        temp_consensus_round1 = os.path.join(classification_dir, "temp_consensus_round1.fasta")

        # Define temporary file to store rest sequence for round2 cd-hit-est
        temp_consensus_round2_input = os.path.join(classification_dir, "temp_consensus_round2_input.fasta")

        # Write sequences to files
        with open(temp_consensus_round1, "w") as high_quality_file, \
                open(temp_consensus_round2_input, 'w') as round2_file:

            for seq_record in consensus_sequences:

                # Sequence names in best_sequences and sequence_for_round2 dont' contain classification
                seq_id = seq_record.id.split("#")[0]

                if seq_id in best_sequences:
                    SeqIO.write(seq_record, high_quality_file, "fasta")
                elif seq_id in sequence_for_round2:
                    SeqIO.write(seq_record, round2_file, 'fasta')

        # Do second round cd-hit-est based on temp_consensus_round2_input
        cd_hit_merge_output_round2 = os.path.join(classification_dir, "TE_Trimmer_consensus_merged_round2.fasta")

        # Round 2 merge require that the alignment coverage for the long and short sequence are both greater than 0.8
        # and the similarity is greater than 0.85
        cd_hit_est(temp_consensus_round2_input, cd_hit_merge_output_round2, identity_thr=0.85, aL=0.8, aS=0.8, s=0.8,
                   thread=num_threads)

        # Define merged file
        cd_hit_est_final_merged = os.path.join(output_dir, "TE_Trimmer_consensus_merged.fasta")

        # Combine the two files into merged file
        with open(temp_consensus_round1, 'r') as file1, \
                open(cd_hit_merge_output_round2, 'r') as file2, \
                open(cd_hit_est_final_merged, 'w') as combined_file:
            # Write contents of the first file
            for line in file1:
                combined_file.write(line)

            # Write contents of the second file
            for line in file2:
                combined_file.write(line)

        # Find sequence names that aren't included inside in cd_hit_est_final_merged file
        # Parse the sequences in the original and merged files
        original_sequences = SeqIO.parse(final_con_file, "fasta")
        merged_sequences = SeqIO.parse(cd_hit_est_final_merged, "fasta")

        # Extract sequence IDs from both files and store to set
        original_ids = {seq_record.id.split("#")[0] for seq_record in original_sequences}
        merged_ids = {seq_record.id.split("#")[0] for seq_record in merged_sequences}

        # Find the difference between the two sets to get sequence names not included in the merged file
        missing_ids = original_ids - merged_ids
        click.echo(missing_ids)

        # Based on missing_ids delete files in proof annotation folder and HMM folder
        for missing_id in missing_ids:
            evaluation_leve = sequence_info[missing_id]["evaluation"]

            # Addd # to the end of missing_id
            missing_id = f"{missing_id}#"
            if evaluation_leve == "Perfect":
                remove_files_with_start_pattern(perfect_proof, missing_id, if_seq_name=False)
            elif evaluation_leve == "Good":
                remove_files_with_start_pattern(good_proof, missing_id, if_seq_name=False)
            elif evaluation_leve == "Reco_check":
                remove_files_with_start_pattern(intermediate_proof, missing_id, if_seq_name=False)
            elif evaluation_leve == "Need_check":
                remove_files_with_start_pattern(need_check_proof, missing_id, if_seq_name=False)
            else:
                remove_files_with_start_pattern(low_copy_dir, missing_id, if_seq_name=False)

        if hmm:
            remove_files_with_start_pattern(hmm_dir, missing_ids)
        click.echo("\nFinished to remove sequence duplications.\n")

    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"Final cd-hit-est deduplication error\n")
            f.write(tb_content + '\n\n')
            click.echo("\nThe final cd-hit-est merge step can't be performed. Please remove duplicated sequence"
                       " by yourself\n")

    #####################################################################################################
    # Code block: Whole genome TE annotation
    #####################################################################################################

    # Delete MSA_dir and Classification if they are empty
    """
    if not os.listdir(MSA_dir):
        os.rmdir(MSA_dir)

    if not os.listdir(classification_dir):
        os.rmdir(classification_dir)
    """
    try:
        # If 90% of the query sequences are processed, RepeatMasker is allowed to be performed whole genome annotation
        # if processed_count >= single_fasta_n * 0.9:

        # Run RepeatMasker
        if genome_anno:
            click.echo("\nTE Trimmer is performing whole genome TE annotation by RepeatMasker\n")

            if os.path.exists(cd_hit_est_final_merged):
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

        # else:
        # click.echo(
        # "Less than 90% of the query sequences processed, TE Trimmer can't perform whole genome TE annotation")
    except Exception as e:
        with open(error_files, "a") as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write(f"Genome TE annotation error.\n")
            f.write(tb_content + '\n\n')
            click.echo("Final genome annotation can't be performed. This won't affect the TE consensus library.")

    end_time = datetime.now()
    duration = end_time - start_time

    # Remove microseconds from the duration
    duration_without_microseconds = timedelta(days=duration.days, seconds=duration.seconds)

    # if not debug:
    # Remove all single files when all the sequences are processed
    # shutil.rmtree(single_file_dir)

    analyze.printProgressBar(processed_count, single_fasta_n, prefix='Progress:', suffix='Complete', length=50,
                             final=True)
    print(f"\nTE Trimmer finished at {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    print(f"TE Trimmer runtime was {duration_without_microseconds}")


# The following is necessary to make the script executable, i.e., python myscript.py.
if __name__ == '__main__':
    main()
