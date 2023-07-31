import os.path
from Bio.Align import AlignInfo
from Bio import AlignIO


class StartEndChecker:
    def __init__(self, input_file, start, end, start_pattern="TGT", end_pattern="ACA", threshold=0.7, ambiguous="X"):
        self.input_file = input_file
        self.start = start
        self.end = end
        self.start_pattern = start_pattern
        self.end_pattern = end_pattern
        self.threshold = threshold
        self.ambiguous = ambiguous
        self.alignment = None
        self.consensus_seq = None

    # Check if file name contain LTR
    def check_LTR(self):
        input_file_name = os.path.basename(self.input_file)
        return 'LTR' in input_file_name

    # Generate consensus sequence
    def con_generater(self):
        self.alignment = AlignIO.read(self.input_file, "fasta")
        summary = AlignInfo.SummaryInfo(self.alignment)
        self.consensus_seq = list(summary.dumb_consensus(threshold=self.threshold, ambiguous=self.ambiguous).upper())

    # Check if the start and end are matchable with the give pattern
    def check_start_end(self):
        MSA_len = self.alignment.get_alignment_length()
        if MSA_len > 1000:  # Only check when sequence length is longer than 1000 bp
            if self.consensus_seq[self.start: self.start + len(self.start_pattern)] == list(self.start_pattern) and \
                    self.consensus_seq[self.end - len(self.end_pattern): self.end] == list(self.end_pattern):
                return True
            else:
                return False

    def check_and_update(self, other_end_pattern=None):
        """
        Check if the start and end of MSA equal to the given patterns
        :param other_end_pattern: list, contains pattern elements
        :return: boolean (if find start and end patterns), start position, end position
        """
        MSA_len = self.alignment.get_alignment_length()

        if MSA_len > 1000:

            # Define other_end_pattern, for example some LTR element is end by AGA
            if other_end_pattern is not None:
                for other_end_patterns in other_end_pattern:
                    other_end_patterns = other_end_patterns.upper()  # Because consensus is upper case.
                    if self.consensus_seq[self.start: self.start + len(self.start_pattern)] == list(self.start_pattern) and \
                            self.consensus_seq[self.end - len(other_end_patterns): self.end] == list(other_end_patterns):
                        return True, self.start, self.end

            if self.consensus_seq[self.start: self.start + len(self.start_pattern)] == list(self.start_pattern) and \
                    self.consensus_seq[self.end - len(self.end_pattern): self.end] == list(self.end_pattern):
                return True, self.start, self.end

            if self.consensus_seq[self.start: self.start + len(self.start_pattern)] != list(self.start_pattern):
                start_window = self.consensus_seq[max(0, self.start - 15): min(self.start + 15, MSA_len)]
                for i in range(len(start_window) - len(self.start_pattern) + 1):
                    if start_window[i: i + len(self.start_pattern)] == list(self.start_pattern):
                        self.start = max(0, self.start - 15) + i
                        break
            if self.consensus_seq[self.end - len(self.end_pattern): self.end] != list(self.end_pattern):
                end_window = self.consensus_seq[max(0, self.end - 15): min(self.end + 15, MSA_len)]
                for i in reversed(range(len(end_window) - len(self.end_pattern) + 1)):
                    if end_window[i: i + len(self.end_pattern)] == list(self.end_pattern):
                        self.end = max(0, self.end - 15) + i + len(self.end_pattern)
                        break
            if self.consensus_seq[self.start: self.start + len(self.start_pattern)] == list(self.start_pattern) and \
                    self.consensus_seq[self.end - len(self.end_pattern): self.end] == list(self.end_pattern):
                return True, self.start, self.end
        return False, self.start, self.end





