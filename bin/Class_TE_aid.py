import subprocess
import os

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


    def run(self):

        TE_aid = os.path.join(self.TE_aid_dir, "TE-Aid")

        # Make a folder to store TE_aid result.
        TE_aid_output_dir = os.path.join(self.output_dir, f"{os.path.basename(self.input_file)}_TEaid_folder")
        if not os.path.isdir(TE_aid_output_dir):
            os.makedirs(TE_aid_output_dir)

        command = [
            TE_aid,
            "-q", self.input_file,
            "-g", self.genome_file,
            "-o", TE_aid_output_dir,
            "-m", str(self.min_orf),
            "-f", str(self.full_length_threshold)
            ]

        subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        final_pdf_file = os.path.join(TE_aid_output_dir, f"{os.path.basename(self.input_file)}.c2g.pdf")

        return final_pdf_file