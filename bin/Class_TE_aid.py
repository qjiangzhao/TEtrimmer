import subprocess
import os
from Class_blast_extension_mafft import SequenceManipulator
import click

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

    def run(self, low_copy = False):

        change_permission_object = SequenceManipulator()
        # Change permissions of the directory and all its content to 755
        # 755 in octal corresponds to rwxr-xr-x
        change_permission_object.change_permissions_recursive(self.TE_aid_dir, 0o755)

        TE_aid = os.path.join(self.TE_aid_dir, "TE-Aid")

        # Check if TE_aid exists
        if not os.path.exists(TE_aid):
            raise FileNotFoundError(f"The TE-Aid executable at {TE_aid} does not exist.")

        # Change TE_aid permission
        # Because I added bash before TE_aid, this can be skipped
        # os.chmod(TE_aid, 0o755)

        # Make a folder to store TE_aid result.
        TE_aid_output_dir = os.path.join(self.output_dir, f"{os.path.basename(self.input_file)}_TEaid_folder")
        if not os.path.isdir(TE_aid_output_dir):
            os.makedirs(TE_aid_output_dir)

        found_match = False

        command = [
            TE_aid,
            "-q", self.input_file,
            "-g", self.genome_file,
            "-o", TE_aid_output_dir,
            "-m", str(self.min_orf),
            "-f", str(self.full_length_threshold)
            ]
        
        if low_copy:
            print (f"{self.input_file} run TE aid low copy")
            command.extend(["-t"])

        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        if result.stderr:
            pass
            #click.echo(f"Error encountered: {self.input_file}{result.stderr.decode('utf-8')}")

        final_pdf_file = os.path.join(TE_aid_output_dir, f"{os.path.basename(self.input_file)}.c2g.pdf")
    
        if low_copy:
            # check the presence of self-alignment of terminal repeats in self-blast.pairs.txt 
            self_blast_txt = os.path.join(TE_aid_output_dir, f"{os.path.basename(self.input_file)}.self-blast.pairs.txt")
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
                    data = [part for part in parts]
                    # Append the data list to the data_list
                    data_list.append(data)
            
            # Iterate through the lists
            for i, lst1 in enumerate(data_list):
                for j, lst2 in enumerate(data_list):
                    # check if the beginning and end aligned to itself. 
                    if i != j and lst1[:2] == lst2[2:] and lst1[2:] == lst2[:2] :
                        found_match = True
                        break
                    if i != j and lst1[0] == lst2[3] and lst1[1] == lst2[2] and lst1[2] == lst2[1] and lst1[3] == lst2[0]:
                        found_match = True
                        break
                if found_match:
                    print("found self-alignment of terminal repeats")
                    break                

        return final_pdf_file, found_match
