import subprocess
import os

class FastaFromBlast:
    """
    A class do local blast and get fasta file.
    """
    def __init__(self, genome, fasta_in, out_dir, script_dir,
                 min_length=100, left_side_extension=1000, right_side_extension=1000):
        """
        :param genome: str, the absolute path of genome file.
        :param fasta_in: str, the absolute path of
        :param out_dir:
        :param script_dir:
        :param min_length:
        :param left_side_extension:
        :param right_side_extension:
        """
        self.genome = genome
        self.fasta_in = fasta_in
        self.min_length = min_length
        self.left_side_extension = left_side_extension
        self.right_side_extension = right_side_extension
        self.out_dir = out_dir
        self.script_dir = script_dir
        self.run()

    def run(self):
        """
        Run function to run bash script "make_fasta_from_blast"
        :return: a file called {file_name}_blast_bed_{left_side_extension}_{right_side_extension}
        """
        # Build the command as a list of strings
        make_fasta_from_blast = os.path.join(self.script_dir, "make_fasta_from_blast.sh")
        command = [make_fasta_from_blast,
                   self.genome,
                   self.fasta_in,
                   str(self.min_length),
                   str(self.left_side_extension),
                   str(self.right_side_extension),
                   self.out_dir,
                   self.script_dir]

        # Run the command
        process = subprocess.run(command, capture_output=True, text=True)

        # Check if the command was successful
        if process.returncode != 0:
            print(f"Error: {process.stderr}")
            return None

        # Return the output of the command
        return process.stdout
