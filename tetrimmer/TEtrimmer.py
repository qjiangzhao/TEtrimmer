# Standard library imports
import concurrent.futures
import json
import logging
import os
import shutil
import traceback
import pandas as pd
import warnings
from datetime import datetime, timedelta

import click
from Bio import BiopythonDeprecationWarning

#from .._version import __version__

import analyze
from functions import (
    cd_hit_est,
    decompress_gzip,
    eliminate_curatedlib_by_repeatmasker,
    repeatmasker,
    check_tools,
    init_logging,
    get_genome_length
)

# Suppress all deprecation warnings
warnings.filterwarnings('ignore', category=BiopythonDeprecationWarning)

#####################################################################################################
# Code block: Import JSON species_config file and define the default parameters
#####################################################################################################

# Load species-specific default values from the JSON config file
config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.json')
# Load the JSON configuration file
with open(config_path, 'r') as config_file:
    preset_config = json.load(config_file)

#####################################################################################################
# Code block: Main functions of TEtrimmer
#####################################################################################################


@click.command(
    context_settings={'max_content_width': 120},
    help=f"""\b
               ##########################################################################################
               \b
                ████████\\ ████████\\ ██\\               ██\\
                \\__██  __|██  _____|██ |              \\__|
                   ██ |   ██ |    ██████\\    ██████\\  ██\\ ██████\\████\\  ██████\\████\\   ██████\\   ██████\\
                   ██ |   █████\\  \\_██  _|  ██  __██\\ ██ |██  _██  _██\\ ██  _██  _██\\ ██  __██\\ ██  __██\\
                   ██ |   ██  __|   ██ |    ██ |  \\__|██ |██ / ██ / ██ |██ / ██ / ██ |████████ |██ |  \\__|
                   ██ |   ██ |      ██ |██\\ ██ |      ██ |██ | ██ | ██ |██ | ██ | ██ |██   ____|██ |
                   ██ |   ████████\\ \\████  |██ |      ██ |██ | ██ | ██ |██ | ██ | ██ |\\███████\\ ██ |
                   \\__|   \\________| \\____/ \\__|      \\__|\\__| \\__| \\__|\\__| \\__| \\__| \\_______|\\__|


                Version: 1.5.1

                Github: https://github.com/qjiangzhao/TEtrimmer

                Developers:

                Jiangzhao Qian;      RWTH Aachen University;                Email: jqian@bio1.rwth-aachen.de

                Hang Xue;            University of California, Berkeley;    Email: hang_xue@berkeley.edu

                Stefan Kusch;        Research Center Juelich;               Email: s.kusch@fz-juelich.de

                Funding source:
                Ralph Panstruga Lab; RWTH Aachen University;                Email: panstruga@bio1.rwth-aachen.de
                Website: https://www.bio1.rwth-aachen.de/PlantMolCellBiology/index.html

                ##########################################################################################

                python ./path_to_TEtrimmer_folder/TEtrimmer.py -i <TE_consensus_file> -g <genome_file>

                TEtrimmer is designed to automate the manual curation of transposable elements (TEs).

                Two mandatory arguments are required, including
                <genome file>, the genome FASTA file, and
                <TE consensus file> from TE discovery software like RepeatModeler, EDTA, or REPET.
                TEtrimmer can do BLAST, sequence extension, multiple sequence alignment (MSA), MSA clustering,
                MSA cleaning, TE boundaries definition, and MSA visualization.
""",
)
@click.option(
    '--input_file',
    '-i',
    required=True,
    type=str,
    help='Path to TE consensus library file (FASTA format). Use the output from RepeatModeler, EDTA, REPET, et al.',
)
@click.option(
    '--genome_file',
    '-g',
    required=True,
    type=str,
    help='Path to genome FASTA file (FASTA format).',
)
@click.option(
    '--output_dir',
    '-o',
    default=os.getcwd(),
    type=str,
    help='Path to output directory. Default: current working directory.',
)
@click.option(
    '--preset',
    '-s',
    default='conserved',
    type=click.Choice(preset_config.keys()),
    help='Choose one preset config (conserved or divergent).',
)
# @click.option('--engine', '-e', default='blast', type=click.Choice(["blast", "mmseqs"]),
#             help='Select the similar sequence search engine. "blast" or "mmseqs". Default: blast')
@click.option(
    '--num_threads',
    '-t',
    default=1,
    type=int,
    help='Thread number used for TEtrimmer. Default: 1',
)
@click.option(
    '--classify_unknown',
    default=False,
    is_flag=True,
    help='Use RepeatClassifier to classify the consensus sequence if the input sequence is not classified or '
    'is unknown or the processed sequence length by TEtrimmer is 2000 bp longer or shorter '
    'than the query sequence.',
)
@click.option(
    '--classify_all',
    default=False,
    is_flag=True,
    help='Use RepeatClassifier to classify every consensus sequence. WARNING: This may take a long time.',
)
@click.option(
    '--continue_analysis',
    '-ca',
    default=False,
    is_flag=True,
    help='Continue from previous unfinished TEtrimmer run and would use the same output directory.',
)
@click.option(
    '--dedup',
    default=False,
    is_flag=True,
    help='Remove duplicate sequences in the input file.',
)
@click.option(
    '--curatedlib',
    default=None,
    type=str,
    help='Path to manually curated high-quality TE consensus library file. TEtrimmer eliminates TE consensus '
    'sequence from "--input_file" if the sequence shares more than 95% identity with sequences from '
    '"--curatedlib".',
)
@click.option(
    '--genome_anno',
    '-ga',
    default=False,
    is_flag=True,
    help='Perform genome TE annotation using RepeatMasker with the TEtrimmer curated TE libraries.',
)
@click.option(
    '--hmm',
    default=False,
    is_flag=True,
    help='Generate HMM files for each processed consensus sequence.',
)
@click.option(
    '--debug',
    default=False,
    is_flag=True,
    help='debug mode. This will keep all raw files. WARNING: Many files will be generated.',
)
@click.option(
    '--fast_mode',
    default=False,
    is_flag=True,
    help='Reduce running time at the cost of lower accuracy and specificity.',
)
@click.option(
    '--pfam_dir',
    '-pd',
    default=None,
    type=str,
    help='Pfam database directory. TEtrimmer checks the existence of Pfam database in the provided path and '
    'downloads it automatically when it is not found. By default, Pfam will be downloaded to TEtrimmer source code folder. '
    'For "singularity" user, please use this option to define a local path, TEtrimmer will download the '
    'database to the provided path if Pfam database is not found. If the automatic download fails, you can'
    'download Pfam database by yourself.',
)
@click.option(
    '--cons_thr',
    type=float,
    help='Threshold used to generate final consensus sequences from MSAs. Default: 0.8',
)
@click.option(
    '--mini_orf',
    type=int,
    help='Define the minimum ORF length to be predicted by TEtrimmer. Default: 200',
)
@click.option(
    '--max_msa_lines',
    type=int,
    help='Set the maximum number of sequences to be included in a multiple sequence alignment. Default: 100',
)
@click.option(
    '--top_msa_lines',
    type=int,
    help='If the sequence number of multiple sequence alignment (MSA) is greater than <max_msa_lines>, '
    'TEtrimmer will first sort sequences by length and choose <top_msa_lines> number of sequences. '
    'Then, TEtrimmer will randomly select sequences from all remaining BLAST hits until <max_msa_lines>'
    'sequences are found for the multiple sequence alignment. Default: 100',
)
@click.option(
    '--min_seq_num',
    type=int,
    help='The minimum blast hit number required for the input sequence. We do not recommend decreasing this number. '
    'Default: 10',
)
@click.option(
    '--min_blast_len',
    type=int,
    help='The minimum sequence length for blast hits to be included for further analysis. Default: 150',
)
@click.option(
    '--max_cluster_num',
    default=5,
    type=int,
    help='The maximum number of clusters assigned in each multiple sequence alignment. '
    'Each multiple sequence alignment can be grouped into different clusters based on alignment patterns '
    'WARNING: using a larger number will potentially result in more accurate consensus results but will '
    'also increase the running time. Default: 5',
)
@click.option(
    '--ext_thr',
    type=float,
    help='The threshold to call “N” at a position. For example, if the most conserved nucleotide in a MSA column'
    'has proportion smaller than <ext_thr>, a “N” will be called at this position. Used with <ext_check_win>. '
    'The lower the value of <ext_thr>, the more likely to get longer the extensions on both ends. '
    'You can try reducing <ext_thr> if TEtrimmer fails to find full-length TEs. Default: 0.7',
)
@click.option(
    '--ext_check_win',
    type=int,
    help='the check windows size during defining start and end of the consensus sequence based on the multiple '
    'sequence alignment. Used with <ext_thr>. If <ext_check_win> bp at the end of multiple sequence alignment '
    'has “N” present (ie. positions have similarity proportion smaller than <ext_thr>), the extension will stop, '
    'which defines the edge of the consensus sequence. Default: 150',
)
@click.option(
    '--ext_step',
    type=int,
    help='the number of nucleotides to be added to the left and right ends of the multiple sequence alignment in each '
    'extension step. TE_Trimmer will iteratively add <ext_step> nucleotides until finding the TE boundary or '
    'reaching <max_ext>. Default: 1000',
)
@click.option(
    '--max_ext',
    type=int,
    help='The maximum extension in nucleotides at both ends of the multiple sequence alignment. Default: 7000',
)
@click.option(
    '--gap_thr',
    type=float,
    help='If a single column in the multiple sequence alignment has a gap proportion larger than <gap_thr> '
    'and the proportion of the most common nucleotide in this column is less than <gap_nul_thr>, '
    'this column will be removed from the consensus. Default: 0.4',
)
@click.option(
    '--gap_nul_thr',
    type=float,
    help='The nucleotide proportion threshold for keeping the column of the multiple sequence alignment. '
    'Used with the <gap_thr> option. i.e. if this column has <40% gap and the portion of T (or any other) nucleotide '
    'is >70% in this particular column, this column will be kept. Default: 0.7',
)
@click.option(
    '--crop_end_div_thr',
    type=float,
    help='The crop end by divergence function will convert each nucleotide in the multiple sequence '
    'alignment into a proportion value. This function will iteratively choose a sliding window from '
    'each end of each sequence of the MSA and sum up the proportion numbers in this window. '
    'The cropping will continue until the sum of proportions is larger than <--crop_end_div_thr>. '
    'Cropped nucleotides will be converted to -. Default: 0.7',
)
@click.option(
    '--crop_end_div_win',
    type=int,
    help='Window size used for the end-cropping process. Used with the <--crop_end_div_thr> option. Default: 40',
)
@click.option(
    '--crop_end_gap_thr',
    type=float,
    help='The crop end by gap function will iteratively choose a sliding window from each end of each sequence '
    'of the MSA and calculate the gap proportion in this window. The cropping will continue until the sum '
    'of gap proportions is smaller than <--crop_end_gap_thr>. Cropped nucleotides will be converted to -. '
    'Default: 0.1',
)
@click.option(
    '--crop_end_gap_win',
    type=int,
    help='Define window size used to crop end by gap. Used with the <--crop_end_gap_thr> option. Default: 250',
)
@click.option(
    '--start_patterns',
    type=str,
    default = 'TG',
    help='LTR elements always start with a conserved sequence pattern. TEtrimmer searches the '
    'beginning of the consensus sequence for these patterns. If the pattern is not found, '
    'TEtrimmer will extend the search of <--start_patterns> to up to 15 nucleotides from the '
    'beginning of the consensus sequence and redefine the start of the consensus sequence '
    'if the pattern is found. Note: The user can provide multiple LTR start patterns in a '
    'comma-separated list, like: TG,TA,TC (no spaces; the order of patterns determines '
    'the priority for the search). Default: TG',
)
@click.option(
    '--end_patterns',
    type=str,
    default = 'CA',
    help='LTR elements always end with a conserved sequence pattern. TEtrimmer searches the '
    'end of the consensus sequence for these patterns. If the pattern is not found, '
    'TEtrimmer will extend the search of <--end_patterns> to up to 15 nucleotides from the '
    'end of the consensus sequence and redefine the end of the consensus sequence '
    'if the pattern is found. Note: The user can provide multiple LTR end patterns in a '
    'comma-separated list, like: CA,TA,GA (no spaces; the order of patterns determines '
    'the priority for the search). Default: CA',
)
@click.option(
    '--poly_patterns',
    type=str,
    default = 'A',
    help="The 3' end of LINE and SINE elements often contains characteristic sequences such as poly(A), "
         "poly(T), or short tandem repeats. TEtrimmer identifies the presence of those feature sequences "
         "to help to define the 3' end boundary of LINE or SINE elements. "
         "You can provide multiple end patterns in a comma-separate list, like: A,T,TA (No space; the order of "
         "patterns determines the priority for the search). Default: A"
)
@click.option(
    '--poly_len',
    type=int,
    help='Define the minimum length requirement of the poly pattern from the parameter --poly_patterns. Default: 10'
)
@click.option(
    '--define_perfect',
    type=int,
    default = 30,
    help='Define the minimum copy number that the output TE consensus sequence can be evaluated as "Perfect". Default: 30'
)
@click.option(
    '--logfile',
    '-l',
    default=None,
    type=str,
    help='Path to log file.'
)
@click.option(
    '--loglevel',
    '-ll',
    default='INFO',
    type=str,
    help='Log level. [DEBUG, INFO, WARNING, ERROR, CRITICAL]',
)
@click.version_option("1.5.1", prog_name='TEtrimmer')
def main(
    input_file,
    genome_file,
    output_dir,
    continue_analysis,
    pfam_dir,
    min_blast_len,
    num_threads,
    max_msa_lines,
    top_msa_lines,
    min_seq_num,
    max_cluster_num,
    cons_thr,
    ext_thr,
    ext_step,
    max_ext,
    gap_thr,
    gap_nul_thr,
    crop_end_div_thr,
    crop_end_div_win,
    crop_end_gap_thr,
    crop_end_gap_win,
    start_patterns,
    end_patterns,
    mini_orf,
    preset,
    ext_check_win,
    dedup,
    genome_anno,
    hmm,
    debug,
    fast_mode,
    classify_unknown,
    classify_all,
    curatedlib,
    poly_patterns,
    poly_len,
    define_perfect,
    logfile,
    loglevel,
):
    perfect_seq_num = define_perfect

    # Add this to click options if mmseq2 has been fully tested
    engine = 'blast'

    # Check for required programs.
    required_tools = [
        'bedtools',
        'samtools',
        'cd-hit-est',
        'ps2pdf',
        'polydot',
        'mafft',
        'iqtree',
        'getorf',
        'hmmsearch',
        'pfam_scan.pl',
        'repeatmasker',
        'repeatmodeler'
    ]
    optional_tools = ['blastn', 'makeblastdb', 'mmseqs']

    if engine == 'blast':
        required_tools.append('blastn')
        required_tools.append('makeblastdb')
        optional_tools = ['mmseqs']
    elif engine == 'mmseqs':
        required_tools.append('mmseqs')
        optional_tools = ['blastn','makeblastdb']

    check_tools(required_tools=required_tools)

    # Set the log file path
    if not logfile:
        # Create output directory if it does not exist
        os.makedirs(output_dir, exist_ok=True)
        # Set the log file path
        logfile = os.path.join(output_dir, 'TEtrimmer.log')

    # Initialize the logging system
    init_logging(loglevel=loglevel, logfile=logfile)

    # Set plot_query, plot_skip, and fast_mode to true
    plot_query = True
    plot_skip = True
    fast_mode = True
    start_time = datetime.now()
    logging.info(f'TEtrimmer started at {start_time.strftime("%Y-%m-%d %H:%M:%S")}.\n')

    #####################################################################################################
    # Code block: File ends
    #####################################################################################################
    # _fm.bed: final multiple sequence alignment.
    # _u.bed: uniqueness.
    # f.bed: When blast hit more than 100, TEtrimmer will only select the top 100 long for MSA
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
    # _proof_anno_me.fa proof curation merged
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
    # fasta_rc reverse complementary fasta file
    # _nm.fa fasta header name modified

    #####################################################################################################
    # Code block: Change permissions of Aliview and TE_Aid
    #####################################################################################################

    # Change TE_Aid permission
    TE_aid_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'TE-Aid-master'
    )
    # Change permissions of the directory and all its content to 755
    # 755 in octal corresponds to rwxr-xr-x
    change_permission = analyze.change_permissions_recursive(TE_aid_path, 0o775)
    if not change_permission:
        pass

    #####################################################################################################
    # Code block: Define the default options according to the given species
    #####################################################################################################

    # Define the default values for the parameters
    default_values = preset_config.get(preset, {})

    # Update the given parameters with default values where user input is None
    if cons_thr is None:
        cons_thr = default_values.get('cons_thr')

    if max_msa_lines is None:
        max_msa_lines = default_values.get('max_msa_lines')

    if top_msa_lines is None:
        top_msa_lines = default_values.get('top_msa_lines')

    if min_seq_num is None:
        min_seq_num = default_values.get('min_seq_num')
        if min_seq_num < 10:
            min_seq_num = 10

    if min_blast_len is None:
        min_blast_len = default_values.get('min_blast_len')

    if ext_thr is None:
        ext_thr = default_values.get('ext_thr')

    if ext_step is None:
        ext_step = default_values.get('ext_step')

    if max_ext is None:
        max_ext = default_values.get('max_ext')

    if gap_thr is None:
        gap_thr = default_values.get('gap_thr')

    if gap_nul_thr is None:
        gap_nul_thr = default_values.get('gap_nul_thr')

    if crop_end_div_thr is None:
        crop_end_div_thr = default_values.get('crop_end_div_thr')

    if crop_end_div_win is None:
        crop_end_div_win = default_values.get('crop_end_div_win')

    if crop_end_gap_thr is None:
        crop_end_gap_thr = default_values.get('crop_end_gap_thr')

    if crop_end_gap_win is None:
        crop_end_gap_win = default_values.get('crop_end_gap_win')

    if poly_len is None:
        poly_len = default_values.get("poly_len")

    if start_patterns is None:
        start_patterns = default_values.get('start_patterns')

    if end_patterns is None:
        end_patterns = default_values.get('end_patterns')

    if mini_orf is None:
        mini_orf = default_values.get('mini_orf')

    if ext_check_win is None:
        ext_check_win = default_values.get('ext_check_win')

    # Log all used parameters
    parameters = {
        "input_file": input_file,
        "genome_file": genome_file,
        "output_dir": output_dir,
        "continue_analysis": continue_analysis,
        "pfam_dir": pfam_dir,
        "min_blast_len": min_blast_len,
        "num_threads": num_threads,
        "max_msa_lines": max_msa_lines,
        "top_msa_lines": top_msa_lines,
        "min_seq_num": min_seq_num,
        "max_cluster_num": max_cluster_num,
        "cons_thr": cons_thr,
        "ext_thr": ext_thr,
        "ext_step": ext_step,
        "max_ext": max_ext,
        "gap_thr": gap_thr,
        "gap_nul_thr": gap_nul_thr,
        "crop_end_div_thr": crop_end_div_thr,
        "crop_end_div_win": crop_end_div_win,
        "crop_end_gap_thr": crop_end_gap_thr,
        "crop_end_gap_win": crop_end_gap_win,
        "start_patterns": start_patterns,
        "end_patterns": end_patterns,
        "mini_orf": mini_orf,
        "preset": preset,
        "ext_check_win": ext_check_win,
        "dedup": dedup,
        "genome_anno": genome_anno,
        "hmm": hmm,
        "debug": debug,
        "fast_mode": fast_mode,
        "classify_unknown": classify_unknown,
        "classify_all": classify_all,
        "curatedlib": curatedlib,
        "poly_patterns": poly_patterns,
        "poly_len": poly_len,
        "define_perfect": define_perfect,
        "logfile": logfile,
        "loglevel": loglevel,
    }

    logging.info("TEtrimmer parameters:")
    for key, value in parameters.items():
        logging.info(f"  {key}: {value}")

    #####################################################################################################
    # Code block: Define input file, output directory, genome
    #####################################################################################################
    # If genome_file is gzipped make a copy of the genome file and unzip it
    # Check if the genome file is gzipped
    is_gzipped = genome_file.endswith('.gz')

    if is_gzipped:
        decompressed_genome_file = decompress_gzip(genome_file)
    else:
        decompressed_genome_file = genome_file

    try:
        (
            bin_py_path,
            output_dir,
            single_file_dir,
            MSA_dir,
            classification_dir,
            hmm_dir,
            proof_curation_dir,
            low_copy_dir,
            perfect_proof,
            good_proof,
            intermediate_proof,
            need_check_proof,
            progress_file,
            pfam_dir,
            final_con_file,
            final_con_file_no_low_copy,
            final_unknown_con_file,
            final_classified_con_file,
            error_files,
            input_file,
            decompressed_genome_file,
            skipped_dir,
            cluster_proof_anno_dir,
        ) = analyze.create_dir(
            continue_analysis,
            hmm,
            pfam_dir,
            output_dir,
            input_file,
            decompressed_genome_file,
            plot_skip,
        )
    except Exception:
        return

    #####################################################################################################
    # Code block: Remove duplications in input file if required, generate single FASTA file and check BLAST database
    #####################################################################################################

    # Generate single files when continue_analysis is false
    if not continue_analysis:
        # When --curatedlib is not None, check if the provided file exist.
        if curatedlib is not None:
            if os.path.isfile(curatedlib):
                logging.info(
                    f'User specified curated TE consensus library found: {curatedlib}'
                )

                # Define path to store curatedlib analysis. Store it to classification_dir
                curatedlib_dir = os.path.join(
                    classification_dir, 'Filter_input_file_based_on_curatedlib'
                )
                os.makedirs(curatedlib_dir, exist_ok=True)
                curatedlib_check = eliminate_curatedlib_by_repeatmasker(
                    curatedlib, input_file, curatedlib_dir
                )

                if curatedlib_check:
                    input_file = curatedlib_check

            else:
                logging.error(
                    f'User specified curated TE consensus library not found: {curatedlib}'
                )
                raise FileNotFoundError(
                    f'User specified curated TE consensus library not found: {curatedlib}'
                )

        # Do CD-HIT-EST merge if merge is true and continue_analysis is false
        if dedup:
            logging.info(
                'TEtrimmer is removing input sequences duplications with CD-HIT-EST, this might take some time.'
            )

            # Define the output file
            merge_output = os.path.join(output_dir, f'{input_file}_cd_hit.fa')

            # Remove duplicates
            try:

                cd_hit_est(
                    input_file,
                    merge_output,
                    identity_thr=0.9,
                    aL=0.95,
                    aS=0.95,
                    s=0.9,
                    thread=num_threads,
                )

                # Convert input_file to merged input_file
                input_file = merge_output
                logging.info(
                    f'Finished merging input sequences with CD-HIT-EST.\nOutput file: {merge_output}'
                )

            except Exception as e:
                # Log the error
                logging.error(f'An error occurred: {e}')
                logging.warning(
                    'TEtrimmer cannot perform the de-duplication step by CD-HIT-EST and will use the '
                    'input sequences directly. This may cause a significantly longer running time but '
                    'will not affect the final result.'
                )

                logging.warning(
                    'You can also run CD-HIT-EST separately to remove redundant sequences:\n'
                    'cd-hit-est -i <input_file> -o <output_file> -T <thread number> -c 0.9 '
                    '-aL 0.9 -aS 0.9 -s 0.9 -l 30'
                )

        # Separate FASTA into single files; if FASTA headers contain "/", " " or ":" convert to "_"
        # Call this function to separate to single FASTA files and create objects from input file
        seq_list, single_fasta_n = analyze.separate_sequences(
            input_file, single_file_dir, continue_analysis=False
        )

        logging.info(f'{single_fasta_n} TE sequences are detected from the input file')

        # Check if BLAST database and genome length files are available, otherwise create these in the
        # same dir as the genome or the output directory specified by the user
        # Set idx_dir to None, the genome blast database files will be stored in the same directory with the
        # genome file. This is required because the TE-Aid package needs those files
        (
            mmseqs_database_dir,
            database_dir,
            database_name,
            length_file,
            fai_file,
        ) = analyze.check_database(
            decompressed_genome_file, idx_dir=None, search_type=engine
        )

        blast_database_path = os.path.join(database_dir, database_name)

        # Initial call to print 0% progress
        analyze.printProgressBar(
            0,
            single_fasta_n,
            prefix='Progress:',
            suffix='Complete',
            length=50
        )

    else:
        # Check if it can perform continue analysis
        if not os.listdir(single_file_dir):
            logging.error(
                'TEtrimmer cannot continue analysis. Please make sure the output directory is '
                'the same as in the previous interrupted run.'
            )
            exit(1)

        else:
            logging.info('TEtrimmer will continue analysis based on previous results.')

            (
                mmseqs_database_dir,
                database_dir,
                database_name,
                length_file,
                fai_file,
            ) = analyze.check_database(
                decompressed_genome_file, idx_dir=None, search_type=engine
            )

            blast_database_path = os.path.join(database_dir, database_name)

            # Create seq_list, which contains sequence objects using the single FASTA files.
            seq_list, single_fasta_n = analyze.separate_sequences(
                input_file, single_file_dir, continue_analysis=True
            )

            # Check which sequences have already been processed
            complete_sequences, skipped_count, low_copy_count, classified_pro = (
                analyze.check_progress_file(progress_file)
            )

            # Filter out already complete sequences
            seq_list = [seq for seq in seq_list if seq.name not in complete_sequences]
            logging.warning(
                f'{single_fasta_n - len(seq_list)} sequences have been processed previously.'
            )

    #####################################################################################################
    # Code block: Calculate genome length
    #####################################################################################################
    try:
        genome_length = get_genome_length(decompressed_genome_file)
    
    except Exception as e:
        with open(error_files, 'a') as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write('\nFailed to calculate genome length file.\n')
            f.write(tb_content + '\n\n')
        logging.warning(
            'This does not affect the final TE consensus sequences '
            'Failed to calculate genome length. This will only affect one column of the Summary.txt'
        )
        genome_length = None

    #####################################################################################################
    # Code block: Enable multiple threads
    #####################################################################################################

    analyze_sequence_params = [
        (
            seq,
            decompressed_genome_file,
            MSA_dir,
            min_blast_len,
            min_seq_num,
            max_msa_lines,
            top_msa_lines,
            max_cluster_num,
            cons_thr,
            ext_thr,
            ext_step,
            classification_dir,
            max_ext,
            gap_thr,
            gap_nul_thr,
            crop_end_div_thr,
            crop_end_div_win,
            crop_end_gap_thr,
            crop_end_gap_win,
            start_patterns,
            end_patterns,
            output_dir,
            pfam_dir,
            mini_orf,
            single_fasta_n,
            hmm,
            hmm_dir,
            ext_check_win,
            debug,
            progress_file,
            classify_unknown,
            classify_all,
            final_con_file,
            final_con_file_no_low_copy,
            final_unknown_con_file,
            final_classified_con_file,
            low_copy_dir,
            fast_mode,
            error_files,
            plot_skip,
            skipped_dir,
            plot_query,
            engine,
            proof_curation_dir,
            poly_patterns, 
            poly_len,
            perfect_seq_num,
            database_dir,
            blast_database_path,
            mmseqs_database_dir,
            loglevel,
            logfile
        )
        for seq in seq_list
    ]

    # Using a ProcessPoolExecutor to run the function in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
        executor.map(analyze.analyze_sequence_helper, analyze_sequence_params)

    #####################################################################################################
    # Code block: Check if all sequences are finished
    #####################################################################################################

    # At the end of the program, check if all sequences have been processed
    completed_sequence, skipped_count, low_copy_count, classified_pro = (
        analyze.check_progress_file(progress_file)
    )

    # Calculate the total count
    processed_count = len(completed_sequence)

    if processed_count == single_fasta_n:
        logging.info(
            f'All sequences have been analysed!\n'
            f'In the analysed sequences {skipped_count} are skipped. Note: not all skipped sequences can have '
            f"TE Aid plot in the 'TEtrimmer_for_proof_curation' folder.\n"
            f'In the analysed sequences {low_copy_count} are identified as low copy TE.\n'
        )

    else:
        remaining = single_fasta_n - processed_count
        logging.warning(
            f'\n\n{remaining} sequences have not been analysed.\n'
            f'In the analysed sequences {skipped_count} are skipped. Note: not all skipped sequences can have '
            f"TE Aid plot in the 'TEtrimmer_for_proof_curation' folder.\n"
            f'In the analysed sequences {low_copy_count} are identified as low copy TE.\n'
        )
        logging.warning(
            "You might find the reasons why some sequences were not analysed from the 'error_file.txt' in the "
            "'Multiple_sequence_alignment' directory."
        )

    #####################################################################################################
    # Code block: Finish classifying unknown consensus sequences and write sequences to consensus file
    #####################################################################################################

    # Final RepeatMasker classification is not necessary, skip in case of errors
    try:
        if 0.3 <= classified_pro < 0.99:
            logging.info(
                'TEtrimmer is doing the final classification. It uses the classified TE to classify '
                'Unknown elements.'
            )
            analyze.repeatmasker_classification(
                final_unknown_con_file,
                final_classified_con_file,
                classification_dir,
                num_threads,
                progress_file,
                final_con_file,
                proof_curation_dir,
                perfect_proof,
                good_proof,
                intermediate_proof,
                need_check_proof,
                low_copy_dir,
                hmm,
                hmm_dir,
            )
        elif classified_pro >= 0.99:
            logging.warning(
                "More than 99% TE are classified, TEtrimmer won't classify 'Unknown' TE by classified TE.\n"
            )
        elif classified_pro < 0.3:
            logging.warning(
                "Less than 30% TE are classified, TEtrimmer won't classify 'Unknown' TE by classified TE.\n"
            )

    except Exception as e:
        with open(error_files, 'a') as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write('\nFinal RepeatMasker classification is wrong.\n')
            f.write(tb_content + '\n\n')
        logging.error(f'The final classification module failed with error {e}')
        logging.warning(
            'This does not affect the final TE consensus sequences '
            "You can choose to ignore this error. For traceback content, please refer to 'error_file.txt' "
            "in the 'Multiple_sequence_alignment' directory.\n"
        )

    #####################################################################################################
    # Code block: Delete the empty folder inside Classification_dir
    #####################################################################################################

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
    # Code block: Read sequence information from the Summary.txt file
    #####################################################################################################
    try:
        # Read progress file
        progress_df = pd.read_csv(progress_file)

        # Create a dictionary with sequence names as keys
        summary_sequence_info = {}

        # Sequence_info directory will contain detailed information for each sequence
        for _index, row in progress_df.iterrows():
            sequence_name = row['output_name']
            evaluation = (
                row['evaluation'] if pd.notna(row['evaluation']) else 'Unknown'
            )  # Default value for NaN
            te_type = (
                row['output_TE_type'] if pd.notna(row['output_TE_type']) else 'Unknown'
            )
            length = (
                row['output_length'] if pd.notna(row['output_length']) else 0
            )  # Default value for NaN

            output_genome_cov_len = (
                row['output_genome_cov_len'] if pd.notna(row['output_genome_cov_len']) else 0
            )

            output_te_type = (
                row['output_TE_type'] if pd.notna(row['output_TE_type']) else 0
            )

            summary_sequence_info[sequence_name] = {
                'evaluation': evaluation,
                'te_type': te_type,
                'length': length,
                'output_genome_cov_len': output_genome_cov_len,
                'output_TE_type' : output_te_type
            }

    except Exception as e:

        with open(error_files, 'a') as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write('\nFinal read Summary.txt file error.\n')
            f.write(tb_content + '\n')
        logging.error(
            f'The final reading Summary.txt file failed.\n Error: {e} \n {tb_content}\n'
        )
        exit(1)
    #####################################################################################################
    # Code block: Merge consensus_file to remove output duplications
    #####################################################################################################

    final_merge_success = True
    # Define merged file
    cd_hit_est_final_merged = os.path.join(
        output_dir, 'TEtrimmer_consensus_merged.fasta'
    )

    try:
        logging.info(
            'TEtrimmer is removing sequence duplications. This might take long time when many sequences'
            'are included into the final consensus library. Please be patient!'
        )
        analyze.merge_cons(
            classification_dir,
            final_con_file,
            summary_sequence_info,
            cd_hit_est_final_merged,
            num_threads,
        )  # Do first round of CD-HIT-EST

    except Exception as e:
        final_merge_success = False
        with open(error_files, 'a') as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write('\nFinal CD-HIT-EST deduplication error.\n')
            f.write(tb_content + '\n')
        logging.error(
            'The final CD-HIT-EST merge step cannot be performed. Final TE consensus library redundancy can '
            f'be higher but the sensitivity is not affected. You can remove duplicated sequence by yourself.\n Error: {e}'
        )
        logging.warning(
            "You can choose to ignore CD-HIT-EST error. For traceback output, please refer to 'error_file.txt' "
            "in the 'Multiple_sequence_alignment' directory.\n"
        )
        exit(1)

    #####################################################################################################
    # Code block: Cluster proof curation files
    #####################################################################################################

    try:
        logging.info(
            'TEtrimmer is clustering TE consensus library. This can potentially take long time when many '
            'sequences exist in the consensus library. Please be patient!'
        )
        multi_dotplot_dir = os.path.join(
            classification_dir, 'Multiple_sequence_dotplot'
        )
        os.makedirs(multi_dotplot_dir, exist_ok=True)
        analyze.cluster_proof_anno_file(
            multi_dotplot_dir,
            final_con_file_no_low_copy,
            continue_analysis,
            cluster_proof_anno_dir,
            num_threads,
            summary_sequence_info,
            perfect_proof,
            good_proof,
            intermediate_proof,
            need_check_proof,
            genome_length
        )

        # clear remove_files_with_start_pattern folder
        if not debug and os.path.exists(multi_dotplot_dir):
            shutil.rmtree(multi_dotplot_dir)

    except Exception as e:
        with open(error_files, 'a') as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write('\nFinal clustering of proof curation files failed.\n')
            f.write(tb_content + '\n\n')
        logging.error(f'Final clustering of proof curation files failed with error: \n{e}')
        logging.error('\n' + tb_content + '')
        logging.warning(
            'This does not affect the final TE consensus sequences. But this can heavily complicate the '
            "TE proof curation. If you don't plan to do proof curation, you can choose to ignore "
            'this error.\n'
        )

    #####################################################################################################
    # Code block: Whole-genome TE annotation
    #####################################################################################################

    try:
        # If 90% of the query sequences have been processed, RepeatMasker is allowed to perform whole genome annotation
        # if processed_count >= single_fasta_n * 0.9:

        # Run RepeatMasker
        if genome_anno:
            logging.info(
                'TEtrimmer is performing whole-genome TE annotation by RepeatMasker. This could take a '
                'long time. \nThe final TE consensus library has been completed. You can use it now.'
            )

            if final_merge_success and os.path.exists(cd_hit_est_final_merged):
                repeatmakser_lib = cd_hit_est_final_merged
            else:
                repeatmakser_lib = final_con_file

            # make a new folder for RepeatMasker output
            repeatmasker_dir = os.path.join(output_dir, 'RepeatMasker_result')
            if not os.path.exists(repeatmasker_dir):
                os.mkdir(repeatmasker_dir)
            genome_anno_result = repeatmasker(
                decompressed_genome_file,
                repeatmakser_lib,
                repeatmasker_dir,
                thread=num_threads,
            )
            if genome_anno_result:
                logging.info('\nFinished whole genome TE annotation by RepeatMasker\n')

    except Exception as e:
        logging.error(
            f'Whole-genome TE annotation by RepeatMasker failed with error: {e}'
        )
        with open(error_files, 'a') as f:
            # Get the traceback content as a string
            tb_content = traceback.format_exc()
            f.write('\nGenome TE annotation error.\n')
            f.write(tb_content + '\n\n')

    #####################################################################################################
    # Code block: End
    #####################################################################################################

    # Remove the decompressed genome file if it was created
    if is_gzipped and os.path.isfile(decompressed_genome_file):
        logging.info(
            f'Removing the decompressed genome file: {decompressed_genome_file}'
        )
        os.remove(decompressed_genome_file)

    end_time = datetime.now()
    duration = end_time - start_time

    # Remove microseconds from the duration
    duration_without_microseconds = timedelta(
        days=duration.days, seconds=duration.seconds
    )

    # if not debug:
    # Remove all single files after all sequences have been processed
    # shutil.rmtree(single_file_dir)

    analyze.printProgressBar(
        processed_count,
        single_fasta_n,
        prefix='Progress:',
        suffix='Complete',
        length=50,
        final=True,
    )
    logging.info(
        f'TEtrimmer analysis finished at {start_time.strftime("%Y-%m-%d %H:%M:%S")}.'
    )
    logging.info(f'TEtrimmer runtime was {duration_without_microseconds}.')


# The following is necessary to make the script executable, i.e., python myscript.py.
if __name__ == '__main__':
    main()
