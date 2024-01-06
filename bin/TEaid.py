import subprocess
import os
from Bio import SeqIO
import click
import pandas as pd
from functions import blast, check_terminal_repeat, file_exists_and_not_empty


def check_self_alignment(seq_obj, seq_file, output_dir, genome_file, blast_hits_count, blast_out_file, plot_skip=True):
    """
    "plot_skip" uses TE-Aid to plot the skipped query sequences.
    """
    blast_full_length_n = check_blast_full_length(seq_obj, blast_out_file, identity=85, coverage=0.8,
                                                  min_hit_length=100, te_aid_blast=False, if_low_copy=True)
    TE_aid_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "TE-Aid-master")

    # At least 2 lines need to meet the requirement
    if blast_full_length_n >= 1:
        check_blast = True

        # Check self-alignment of terminal repeats
        TE_aid_object = TEAid(seq_file, output_dir, genome_file, TE_aid_dir=TE_aid_path)
        TE_aid_plot, found_match = TE_aid_object.run(low_copy=True, label=False)

        seq_obj.update_blast_hit_n(blast_hits_count)

        # Convert found_match to True when LTR or TIR is found
        # Require at least 2 full length blast hits for LTR and 1 for TIR
        if (found_match == "LTR" and blast_full_length_n >= 2) or found_match == "TIR":
            found_match_boolean = True
        else:
            found_match_boolean = False

        # Update_low_copy returns 'True' when check_80 and found_match are both true
        check_low_copy = seq_obj.update_low_copy(check_blast, found_match_boolean)

    else:
        if plot_skip:
            # Plot skipped elements if required
            TE_aid_object = TEAid(seq_file, output_dir, genome_file, TE_aid_dir=TE_aid_path)
            TE_aid_plot, found_match = TE_aid_object.run(low_copy=True, label=False)
        else:
            TE_aid_plot = None
        check_low_copy = False
        found_match = False

    return check_low_copy, blast_full_length_n, found_match, TE_aid_plot


def check_blast_full_length(seq_obj, blast_out_file, identity=80, coverage=0.9, min_hit_length=100, te_aid_blast=False,
                            if_low_copy=False):

    if not file_exists_and_not_empty(blast_out_file):
        return 0

    # The TE-Aid BLAST output file has a header
    if te_aid_blast:
        df = pd.read_csv(blast_out_file, sep=r"\s+", skiprows=1, header=None)
    else:
        df = pd.read_csv(blast_out_file, sep=r"\s+", header=None)

    # Extract sequence length
    if if_low_copy:
        seq_length = seq_obj.old_length
    else:
        # If "if_low_copy" is 'False', the sequence is not a low copy element. seq_obj is consi_obj (see "seqclass.py").
        seq_length = seq_obj.new_length

    identity_condition = df[2] > identity
    coverage_condition = df[3] / seq_length >= coverage
    length_condition = df[3] > min_hit_length

    # Filter the DataFrame
    filtered_df = df[identity_condition & coverage_condition & length_condition]
    blast_full_length_n = filtered_df.shape[0]

    return blast_full_length_n


class TEAid:

    def __init__(self, input_file, output_dir, genome_file, TE_aid_dir, error_file=None,
                 min_orf=200, full_length_threshold=0.9):
        """
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
        self.TE_aid_dir = TE_aid_dir
        self.error_file = error_file
        self.min_orf = min_orf
        self.full_length_threshold = full_length_threshold

    def run(self, low_copy=False, label=True):
        # Define TE-Aid software executable path
        TE_aid = os.path.join(self.TE_aid_dir, "TE-Aid")

        # Check if TE-Aid exists
        if not os.path.exists(TE_aid):
            raise FileNotFoundError(f"The TE-Aid executable at {TE_aid} does not exist.")

        # Make a folder to store TE-Aid result.
        TE_aid_output_dir = os.path.join(self.output_dir, f"{os.path.basename(self.input_file)}_TEaid")
        os.makedirs(TE_aid_output_dir, exist_ok=True)

        # TE-Aid will create plot with end name c2g.pdf
        final_pdf_file = os.path.join(TE_aid_output_dir, f"{os.path.basename(self.input_file)}.c2g.pdf")
        found_match = False

        # Check if the PDF c2g.pdf exists. If so, skip the downstream analysis, if low_copy is 'False'.
        # TETrimmer will plot the query sequence using TE Aid. One query file can contain multiple clusters,
        # TE-Aid will test if this has been created before.
        if not low_copy and os.path.exists(final_pdf_file):
            return final_pdf_file, found_match

        command = [
            TE_aid,
            "-q", self.input_file,
            "-g", self.genome_file,
            "-o", TE_aid_output_dir,
            "-m", str(self.min_orf),
            "-f", str(self.full_length_threshold)
            ]

        # If it is low copy element, add '-t' option to enable self-BLAST from TE-Aid file
        if low_copy:
            command.extend(["-T"])
        if label:
            command.extend(["-TM"])

        try:
            subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        except subprocess.CalledProcessError as e:
            if self.error_file is not None:
                with open(self.error_file, 'a') as f:
                    f.write(f"\nTE Aid error for {os.path.basename(self.input_file)} with error code {e.returncode}")
                    f.write('\n' + e.stderr + '\n')
            pass

        if low_copy:

            #####################################################################################################
            # Code block: Check terminal repeat
            #####################################################################################################
            # Check the presence of self-alignment of terminal repeats in self-blast.pairs.txt
            self_blast_txt = os.path.join(TE_aid_output_dir, f"{os.path.basename(self.input_file)}.self-blast.pairs.txt")

            if not os.path.exists(self_blast_txt):
                blast_out = None
            else:
                blast_out = self_blast_txt

            LTR_boundary, TIR_boundary, found_match = check_terminal_repeat(self.input_file, TE_aid_output_dir,
                                                                            teaid_blast_out=blast_out, TIR_adj=30,
                                                                            LTR_adj=50)

        return final_pdf_file, found_match

    #####################################################################################################
    # Code block: Check number of BLAST full-length hits
    #####################################################################################################
    def check_blast_full_n(self, seq_obj, engine="blast"):

        # Make a folder to store TE-Aid result.
        TE_aid_output_dir = os.path.join(self.output_dir, f"{os.path.basename(self.input_file)}_TEaid")
        if not os.path.isdir(TE_aid_output_dir):
            os.makedirs(TE_aid_output_dir)

        te_aid_blast_file = os.path.join(TE_aid_output_dir, "blastn.txt")

        # If blast.txt file was found, use the TE-Aid output directly. Otherwise, do BLAST search. This may save resources.
        if os.path.exists(te_aid_blast_file):
            full_length_n = check_blast_full_length(seq_obj, te_aid_blast_file, identity=85, coverage=0.9,
                                                    min_hit_length=100, te_aid_blast=True)
        else:
            bed_out_file, blast_hits_count, blast_out_file = blast(self.input_file, self.genome_file,
                                                                   self.output_dir, search_type=engine)
            full_length_n = check_blast_full_length(seq_obj, blast_out_file, identity=85, coverage=0.9,
                                                    min_hit_length=100, te_aid_blast=False)

        return full_length_n











