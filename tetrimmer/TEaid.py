import logging
import os
import subprocess
import traceback

from functions import (
    blast,
    check_terminal_repeat,
    check_blast_full_length,
    sum_non_overlapping_lengths,
    file_exists_and_not_empty
)


def low_copy_full_blast_and_terminal_check_plus_teaid_plotting(
    seq_obj,
    seq_file,
    output_dir,
    genome_file,
    blast_database_path,
    mmseqs_database_dir,
    blast_hits_count,
    blast_out_file,
):
    # Check full length blast for the input sequence
    all_blast_hit_n, blast_full_length_n = check_blast_full_length(
        seq_obj,
        blast_out_file,
        identity=75,
        coverage=0.8,
        min_hit_length=50,
        te_aid_blast=False,
        check_query=True
    )

    try:
        # Calculate input TE genome coverage length
        te_genome_coverage_len = sum_non_overlapping_lengths(blast_out_file, te_aid_blast=False)

    except Exception as e:
        te_genome_coverage_len = 'NaN'

    TE_aid_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'TE-Aid-master'
    )
    # At least 2 lines need to meet the requirement
    if blast_full_length_n >= 1:
        check_blast = True
    else:
        check_blast = False

    # Check self-alignment of terminal repeats
    TE_aid_object = TEAid(seq_file, output_dir, genome_file, blast_database_path,
        mmseqs_database_dir, TE_aid_dir=TE_aid_path)
    TE_aid_plot = TE_aid_object.run(title="before")
    found_match = TE_aid_object.teaid_check_termina_repeat()

    # Convert found_match to True when LTR or TIR is found
    # Require at least 2 full length blast hits for LTR and 1 for TIR
    if (found_match == 'LTR' and blast_full_length_n >= 2) or found_match == 'TIR':
        found_match_boolean = True
    else:
        found_match_boolean = False

    # Update_low_copy returns 'True' when check_80 and found_match are both true
    check_low_copy = seq_obj.update_low_copy(check_blast, found_match_boolean)

    # Update input sequence terminal repeat information
    seq_obj.set_old_terminal_repeat(found_match)

    seq_obj.set_input_genome_cov_len(te_genome_coverage_len)

    return check_low_copy, blast_full_length_n, found_match, TE_aid_plot


class TEAid:
    def __init__(
        self,
        input_file,
        output_dir,
        genome_file,
        blast_database_path,
        mmseqs_database_dir,
        TE_aid_dir,
        error_file=None,
        min_orf=200,
        full_length_threshold=0.9,
    ):
        """
        TE-Aid class offer functions to use TE-Aid sequence blast or self-blast result to check full length blast
        and terminal repeat
        :param input_file: str, absolute path of input file
        :param output_dir: str, absolute directory of output file
        :param genome_file: str, absolute path of genome file
        :param TE_aid_dir: str, absolute path of executable of TE-Aid software
        :param min_orf: num default 200, minimum ORF size
        :param full_length_threshold: num (0-1), default 0.9, threshold to classify as intact TE against consensus sequences
        """
        self.input_file = input_file
        self.output_dir = output_dir
        self.genome_file = genome_file
        self.blast_database_path = blast_database_path
        self.mmseqs_database_dir = mmseqs_database_dir
        self.TE_aid_dir = TE_aid_dir
        self.error_file = error_file
        self.min_orf = min_orf
        self.full_length_threshold = full_length_threshold
        # Make a folder to store TE-Aid result.
        self.TE_aid_output_dir = os.path.join(self.output_dir, f'{os.path.basename(self.input_file)}_TEaid')
        os.makedirs(self.TE_aid_output_dir, exist_ok=True)

    def run(self, title="after", v_x_line_1=0, v_x_line_2=0):
        # Define TE-Aid software executable path
        TE_aid = os.path.join(self.TE_aid_dir, 'TE-Aid')

        # Check if TE-Aid exists
        if not os.path.exists(TE_aid):
            raise FileNotFoundError(f'The TE-Aid executable at {TE_aid} does not exist.')

        # TE-Aid will create plot with end name c2g.pdf
        final_pdf_file = os.path.join(
            self.TE_aid_output_dir, f'{os.path.basename(self.input_file)}.c2g.pdf'
        )

        # Check if the PDF c2g.pdf exists. If so, skip the downstream analysis
        # One query file can contain multiple clusters,
        # TE-Aid will test if this has been created before.
        if os.path.exists(final_pdf_file):
            return final_pdf_file

        command = [
            TE_aid,
            '-q', self.input_file,
            '-g', self.genome_file,
            '-o', self.TE_aid_output_dir,
            '-m', str(self.min_orf),
            '-f', str(self.full_length_threshold),
            '-v_x_line_1', str(v_x_line_1),
            '-v_x_line_2', str(v_x_line_2),
            '-title', str(title),
            '-T'  # keep self-BLAST file from TE-Aid file,

        ]

        try:
            te_aid_result = subprocess.run(
                command,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )

            if te_aid_result.stderr:
                if self.error_file is not None:
                    with open(self.error_file, 'a') as f:
                        f.write(
                            f'\nTE Aid error for {os.path.basename(self.input_file)} with error \n '
                            f'{te_aid_result.stderr}'
                        )

            return final_pdf_file

        except subprocess.CalledProcessError as e:
            if self.error_file is not None:
                with open(self.error_file, 'a') as f:
                    f.write(
                        f'\nTE Aid error for {os.path.basename(self.input_file)} with error code {e.returncode}'
                    )
                    f.write(f'\n{e.stdout}')
                    f.write(f'\n{e.stderr}\n')

    #####################################################################################################
    # Code block: Check terminal repeat
    #####################################################################################################
    def teaid_check_termina_repeat(self):

        # Check the presence of self-alignment of terminal repeats in self-blast.pairs.txt
        # TE-Aid generated self-blast result, use it when available
        self_blast_txt = os.path.join(
            self.TE_aid_output_dir,
            f'{os.path.basename(self.input_file)}.self-blast.pairs.txt',
        )

        # When blast_out is None means TEAid didn't generate the self-blast file successfully.
        # the function check_terminal_repeat will do sequence self-blast
        if not os.path.exists(self_blast_txt):
            blast_out = None
        else:
            blast_out = self_blast_txt

        LTR_boundary, TIR_boundary, found_match = check_terminal_repeat(
            self.input_file,
            self.TE_aid_output_dir,
            teaid_blast_out=blast_out,
        )

        return found_match

    #####################################################################################################
    # Code block: Check number of BLAST full-length hits
    #####################################################################################################
    def check_blast_full_n(self, seq_obj, check_query=True, engine='blast'):

        # query_check determine to use the query input sequence or use TEtrimmer processed sequence for the
        # full length blast number calculation

        te_aid_blast_file = os.path.join(self.TE_aid_output_dir, 'blastn.txt')

        # If blast.txt file was found, use the TE-Aid output directly. Otherwise, do BLAST search.
        if os.path.exists(te_aid_blast_file):
            all_blast_hit_n, full_length_n = check_blast_full_length(
                seq_obj,
                te_aid_blast_file,
                identity=85,
                coverage=0.9,
                min_hit_length=50,
                te_aid_blast=True,
                check_query=check_query
            )
        else:
            bed_out_file, blast_hits_count, blast_out_file = blast(
                self.input_file, self.genome_file, self.blast_database_path,
                self.mmseqs_database_dir, self.output_dir, search_type=engine
            )
            all_blast_hit_n, full_length_n = check_blast_full_length(
                seq_obj,
                blast_out_file,
                identity=85,
                coverage=0.9,
                min_hit_length=50,
                te_aid_blast=False,
                check_query=check_query
            )

        return all_blast_hit_n, full_length_n

    #####################################################################################################
    # Code block: Calculate genome coverage length for each TE
    #####################################################################################################
    def te_genome_coverage(self):

        # query_check determine to use the query input sequence or use TEtrimmer processed sequence for the
        # full length blast number calculation

        te_aid_blast_file = os.path.join(self.TE_aid_output_dir, 'blastn.txt')

        # If blast.txt file was found, use the TE-Aid output directly. Otherwise, do BLAST search.
        if os.path.exists(te_aid_blast_file):
            try:
                te_genome_coverage_len = sum_non_overlapping_lengths(te_aid_blast_file, te_aid_blast=True)
            except Exception:
                te_genome_coverage_len = 'NaN'

        else:
            try:
                bed_out_file, blast_hits_count, blast_out_file = blast(
                    self.input_file, self.genome_file, self.blast_database_path,
                    self.mmseqs_database_dir, self.output_dir
                )

                te_genome_coverage_len = sum_non_overlapping_lengths(blast_out_file, te_aid_blast=False)

            except Exception:
                te_genome_coverage_len = 'NaN'

        return te_genome_coverage_len