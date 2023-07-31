import subprocess

#TODO I can't simply use cd-hit-est to cluster sequences. Because nested TEs are still exist

class CDHitEst:
    """
    A wrapper for the CD-HIT-EST program, a sequence clustering tool.

    Usage:
        input_file = "input_sequences.fasta"
        output_file = "clustered_sequences.fasta"
        cd_hit_est = CDHitEst(input_file, output_file)
        cd_hit_est.run()
    """

    def __init__(self, input_file, output_file, threshold=0.95, word_length=10,
                 global_identity=1, band_width=200, memory_limit=4000,
                 num_threads=1, min_length=80, max_description=0, alignment_coverage=0.8, accurate=False):
        """
        Initialize the CDHitEst class with input and output files, and optional parameters.

        :param input_file: str, path to the input file in FASTA format
        :param output_file: str, path to the output file
        :param threshold: float, sequence identity threshold (default: 0.95)
        :param word_length: int, The word length refers to the size of the short subsequences (k-mers)
            that are used for the initial sequence comparison to speed up the clustering process.
        :param global_identity: use global sequence identity, default 1. If set to 0, then use local sequence identity,
            calculated as: number of identical amino acids or bases in alignment divided by length of the alignment
            NOTE!!! don't use -G 0 unless you use alignment coverage controls; see options -aL, -AL, -aS, -AS
        :param band_width: Using a smaller bandwidth value can make the program run faster, but it might also increase
            the chances of missing some true sequence matches. On the other hand, using a larger bandwidth value can
            increase sensitivity, but it will also result in slower processing times.
        :param memory_limit: memory limit (in MB) for the program, default 4000; 0 for unlimitted
        :param num_threads: number of threads, default 1; with 0, all CPUs will be used
        :param min_length: length of throw_away_sequences, default 80
        :param max_description: length of description in .clstr file, default 0.
            if set to 0, it takes the fasta defline and stops at first space
        :param alignment_coverage: alignment coverage for the shorter sequence, default 0.8.
            if set to 0.9, the alignment must covers 90% of the sequence
        :param accurate: 1 or 0, default 0
            by cd-hit's default algorithm, a sequence is clustered to the first cluster that meet the
            threshold (fast cluster). If set to 1, the program will cluster it into the most similar cluster
            that meet the threshold (accurate but slow mode)
            but either 1 or 0 won't change the representatives of final clusters
        """
        self.input_file = input_file
        self.output_file = output_file
        self.threshold = threshold
        self.word_length = word_length
        self.global_identity = global_identity
        self.band_width = band_width
        self.memory_limit = memory_limit
        self.num_threads = num_threads
        self.min_length = min_length
        self.max_description = max_description
        self.alignment_coverage = alignment_coverage
        self.accurate = accurate

    def run(self):
        """
        Run the CD-HIT-EST program with the provided parameters.
        """
        command = [
            "cd-hit-est",
            "-i", self.input_file,
            "-o", self.output_file,
            "-c", str(self.threshold),
            "-n", str(self.word_length),
            "-G", str(self.global_identity),
            "-b", str(self.band_width),
            "-M", str(self.memory_limit),
            "-T", str(self.num_threads),
            "-l", str(self.min_length),
            "-d", str(self.max_description),
            "-aS", str(self.alignment_coverage),
            "-g", str(self.accurate)
        ]

        subprocess.run(command, check=True)
