import subprocess
import os
from Bio import SeqIO
import click
import pandas as pd
from Class_blast_extension_mafft import SequenceManipulator
from Function_blast_extension_mafft import blast


def check_self_alignment(seq_obj, seq_file, output_dir, genome_file, blast_hits_count, blast_out_file):

    blast_full_length_n = check_blast_low_copy(seq_obj, blast_out_file, identity=90, coverage=0.9, min_hit_length=100,
                                               te_aid_blast=False, if_low_copy=True)

    # At least 2 lines need to meet the requirement
    if blast_full_length_n >= 2:

        check_blast = True

        # Check self-alignment of terminal repeats
        TE_aid_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "TE-Aid-master")

        # Run TE_aid
        TE_aid_object = TEAid(seq_file, output_dir, genome_file, TE_aid_dir=TE_aid_path)
        TE_aid_plot, found_match = TE_aid_object.run(low_copy=True)

        seq_obj.update_blast_hit_n(blast_hits_count)

        # Convert found_match to True when LTR or TIR is found
        if found_match == "LTR" or found_match == "TIR":
            found_match_boolean = True
        else:
            found_match_boolean = False

        # Update_low_copy return true when check_80 and found_match are both true
        check_low_copy = seq_obj.update_low_copy(check_blast, found_match_boolean)

    else:
        check_low_copy = False
        found_match = False

    return check_low_copy, blast_full_length_n, found_match


def check_blast_low_copy(seq_obj, blast_out_file, identity=90, coverage=0.9, min_hit_length=100, te_aid_blast=False,
                         if_low_copy=False):

    # The TE Aid blast output have a header
    if te_aid_blast:
        df = pd.read_csv(blast_out_file, sep="\s+", skiprows=1, header=None)
    else:
        df = pd.read_csv(blast_out_file, sep="\s+", header=None)

    if if_low_copy:
        seq_length = seq_obj.old_length
    else:
        seq_length = seq_obj.new_length

    identity_condition = df[2] > identity
    coverage_condition = df[3] / seq_length > coverage
    length_condition = df[3] > min_hit_length

    # Filter the DataFrame
    filtered_df = df[identity_condition & coverage_condition & length_condition]

    blast_full_length_n = filtered_df.shape[0]

    return blast_full_length_n


class TEAid:

    def __init__(self, input_file, output_dir, genome_file, TE_aid_dir, min_orf=200, full_length_threshold=0.9):
        """
        :param input_file: str, absolute path of input file
        :param output_dir: str, absolute directory of output file
        :param genome_file: str, absolute path of genome file
        :param TE_aid_dir: str, absolute path of executable TE-Aid software
        :param min_orf: num default 200, minimum orf size
        :param full_length_threshold: num (0-1) default 0.9, threshold to classify as intact TE against consensus sequence
        """
        self.input_file = input_file
        self.output_dir = output_dir
        self.genome_file = genome_file
        self.TE_aid_dir = TE_aid_dir
        self.min_orf = min_orf
        self.full_length_threshold = full_length_threshold

    def run(self, low_copy=False):

        # Define TE_Aid software executable path
        TE_aid = os.path.join(self.TE_aid_dir, "TE-Aid")

        # Check if TE_aid exists
        if not os.path.exists(TE_aid):
            raise FileNotFoundError(f"The TE-Aid executable at {TE_aid} does not exist.")

        # Make a folder to store TE_aid result.
        TE_aid_output_dir = os.path.join(self.output_dir, f"{os.path.basename(self.input_file)}_TEaid")
        os.makedirs(TE_aid_output_dir, exist_ok=True)

        found_match = False

        command = [
            TE_aid,
            "-q", self.input_file,
            "-g", self.genome_file,
            "-o", TE_aid_output_dir,
            "-m", str(self.min_orf),
            "-f", str(self.full_length_threshold)
            ]

        # If it is low copy element add -t option to enable to keep self blast file from TE_Aid
        if low_copy:

            command.extend(["-T"])

        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        if result.stderr:
            pass
            #click.echo(f"Error encountered: {self.input_file}\n{result.stderr.decode('utf-8')}")

        final_pdf_file = os.path.join(TE_aid_output_dir, f"{os.path.basename(self.input_file)}.c2g.pdf")
    
        if low_copy:

            #####################################################################################################
            # Code block: Check terminal repeat
            #####################################################################################################

            # Read input file and get sequence length
            record = SeqIO.read(self.input_file, "fasta")
            record_len = len(record.seq)

            # Check the presence of self-alignment of terminal repeats in self-blast.pairs.txt
            self_blast_txt = os.path.join(TE_aid_output_dir, f"{os.path.basename(self.input_file)}.self-blast.pairs.txt")

            # Sometime TE_Aid won't give self blast result
            if not os.path.exists(self_blast_txt):

                database_file = os.path.join(TE_aid_output_dir, "Tem_blast_database")
                makeblastdb_cmd = f"makeblastdb -in {self.input_file} -dbtype nucl -out {database_file}"
                subprocess.run(makeblastdb_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

                blast_cmd = f"blastn -query {self.input_file} -db {database_file} " \
                            f"-outfmt \"6 qseqid qstart qend sstart send \" " \
                            f"-evalue 0.05"

                # Execute the command
                result = subprocess.run(blast_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

                # If there is an error, print it
                if result.returncode != 0:
                    click.echo(
                        f"An error occurred self blast: {os.path.basename(self.input_file)}\n{result.stderr.decode('utf-8')}")
                    return
                else:
                    blast_out = result.stdout.decode('utf-8')

                    with open(self_blast_txt, 'w') as f:

                        # Because TE_Aid result have header, give a header to blast result
                        f.write("qseqid\tqstart\tqend\tsstart\tsend\n")
                        f.write(blast_out)

            data_list = []
            # Read the input file
            with open(self_blast_txt, "r") as file:
                # Skip the header line
                next(file)
                # Iterate through the lines in the file
                for line in file:
                    # Split each line into a list of strings using whitespace as the separator
                    parts = line.strip().split()[1:]
                    # Convert the string elements to integers
                    data = [int(part) for part in parts]
                    # Append the data list to the data_list
                    data_list.append(data)
            
            # Iterate through the lists
            for i, lst1 in enumerate(data_list):
                for j, lst2 in enumerate(data_list):
                    # Check if the beginning and end aligned to itself.
                    # The length of LTR have to be more than 100 bp. The start point of LTR has to be smaller than 15
                    # The length of LTR can't be longer than 1/5 of query sequence
                    if i != j and lst1[0] <= 15 and 100 <= lst1[1] - lst1[0] <= record_len / 5 \
                            and lst1[:2] == lst2[2:] and lst1[2:] == lst2[:2]:
                        found_match = "LTR"
                        break

                    # Find the reverse repeat
                    # Make sure reverse repeat has to be longer than 20 and can't longer than 1/2 of query sequence
                    if i != j and lst1[0] <= 15 and lst1[0] == lst2[3] and lst1[1] == lst2[2] and lst1[2] == lst2[1] \
                            and lst1[3] == lst2[0] and 20 <= lst1[1] - lst1[0] <= record_len / 2:
                        found_match = "TIR"
                        break
                if found_match == "LTR" or found_match == "TIR":
                    break

        return final_pdf_file, found_match

    #####################################################################################################
    # Code block: Check blast full length number
    #####################################################################################################
    def check_blast_full_n(self, seq_obj):

        # Make a folder to store TE_aid result.
        TE_aid_output_dir = os.path.join(self.output_dir, f"{os.path.basename(self.input_file)}_TEaid")
        if not os.path.isdir(TE_aid_output_dir):
            os.makedirs(TE_aid_output_dir)

        te_aid_blast_file = os.path.join(TE_aid_output_dir, "blastn.txt")

        # If blast.txt file is found use the TE Aid output directly. Otherwise, do blast
        if os.path.exists(te_aid_blast_file):

            full_length_n = check_blast_low_copy(seq_obj, te_aid_blast_file, identity=90, coverage=0.9,
                                                 min_hit_length=100, te_aid_blast=True)
        else:
            blast_obj = SequenceManipulator()
            bed_out_file, blast_hits_count, blast_out_file = blast_obj.blast(self.input_file, self.genome_file,
                                                                             self.output_dir)
            full_length_n = check_blast_low_copy(seq_obj, blast_out_file, identity=90, coverage=0.9,
                                                 min_hit_length=100, te_aid_blast=False)

        return full_length_n











