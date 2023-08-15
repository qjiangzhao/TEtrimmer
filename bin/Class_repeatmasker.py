import subprocess

class RepeatMaskerRunner:

    #TODO set -s option into __init__
    #TODO set -div to __int__

    """
    A class to run RepeatMasker on an input file.

      Methods
    -------
    run():
        Executes the RepeatMasker command with specified input and output files.
    """

    def __init__(self, input_file, output_dir, genome_file, lib_file, num_threads=1):
        """
        Initialize the RepeatMaskerRunner with input and output files, and optional parameters.

        :param input_file: str, path to the input file in FASTA format
        :param output_dir: str, path to the output directory
        :param genome_file: str, path to genome file in FASTA format
        :param lib_file: str, path to TE consensus library in FASTA format
            NOTE: consensus file have to follow RepeatMasker format like #LTR/Copia
        :param num_threads: int, number of threads to use for RepeatMasker (default: 1)
        """
        self.input_file = input_file
        self.output_dir = output_dir
        self.genome_file = genome_file
        self.lib_file = lib_file
        self.num_threads = num_threads

    def run(self):
        """
        Run RepeatMasker with the provided parameters.
        """
        # Construct the RepeatMasker command
        command = ["RepeatMasker",
                   self.genome_file,
                   "-lib", self.lib_file,
                   "-pa", str(self.num_threads),
                   "-dir", self.output_dir,
                   "-s",    # Slow search; 0-5% more sensitive, 2-3 times slower than default
                   "-gff",  # Creates an additional Gene Feature Finding format output
                   "-xm",   # Creates an additional output file in cross_match format (for parsing)
                   "-a",    # Writes alignments in .align output file
                   ]



        # Run RepeatMasker using subprocess
        subprocess.run(command, check=True)
