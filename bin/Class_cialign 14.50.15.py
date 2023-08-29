import subprocess

class CIAlignRunner:
    """
    A class to run CIAlign using subprocess.

    Methods
    -------
    run():
        Executes the CIAlign command with specified input and output files.
    """
    def __init__(self, input_file, output_file):
        """
        Initialize the CIAlign with input and output files, and optional parameters.

        :param input_file: str, path to the input file in FASTA format
        :param output_dir: str, path to the output directory
        :param crop_ends_mingap_perc Minimum proportion of the sequence length (excluding gaps) that is the threshold for change in gap numbers.
        """
        self.input_file = input_file
        self.output_file = output_file

    def run(self):
        """
        Executes the CIAlign command with specified input and output files.

        Returns
        -------
        process
            a CompletedProcess instance which holds information about the process run.
        """
        command = ["cialign", self.input_file, self.output_file]
        process = subprocess.run(command, capture_output=True)
        return process





