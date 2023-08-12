import os.path
from Bio.Align import AlignInfo
from Bio import AlignIO


class StartEndChecker:
    """
    Check if the MSA start and end with the given pattern
    """

    def __init__(self, input_file, start, end, start_patterns=None, end_patterns=None, threshold=0.7, ambiguous="X"):
        """
        :param input_file: str, input file directory
        :param start: num, start point
        :param end: num, end point
        :param start_patterns: str default: None, give start patterns separated by comma like ATA, AGA. Note: the order
        of pattern is important because TE Trimmer will start check from the first pattern
        :param end_patterns: str default: None, give end patterns separated by comma like ACA, ATA, AGA
        :param threshold: num 0-1 default: 0.7, threshold to generate consensus sequence
        :param ambiguous: str default: X, this letter will be used for ambiguous consensus letter
        """
        self.input_file = input_file
        self.start = start
        self.end = end
        self.start_patterns = [pattern.upper() for pattern in start_patterns.split(
            ',')] if start_patterns else None  # Ensure all patterns are upper case
        self.end_patterns = [pattern.upper() for pattern in end_patterns.split(',')] if end_patterns else None
        self.threshold = threshold
        self.ambiguous = ambiguous
        self.alignment = None
        self.consensus_seq = None

    # Check if file name contains LTR
    def check_LTR(self):
        input_file_name = os.path.basename(self.input_file)
        return 'LTR' in input_file_name

    # Generate consensus sequence
    def con_generater(self):
        self.alignment = AlignIO.read(self.input_file, "fasta")
        summary = AlignInfo.SummaryInfo(self.alignment)
        self.consensus_seq = list(summary.dumb_consensus(threshold=self.threshold, ambiguous=self.ambiguous).upper())

    # Check if the start and end are matchable with the given pattern
    def check_start_end(self):
        start_matched = end_matched = True

        if self.start_patterns:  # Check if start_patterns is None
            for start_pattern in self.start_patterns:
                if self.consensus_seq[self.start: self.start + len(start_pattern)] == list(start_pattern):
                    break
            else:  # No break, means no match found
                start_matched = False
        if self.end_patterns:  # Check if end_patterns is None
            for end_pattern in self.end_patterns:
                if self.consensus_seq[self.end - len(end_pattern): self.end] == list(end_pattern):
                    break
            else:  # No break, means no match found
                end_matched = False

        return start_matched, end_matched

    def check_and_update(self):
        """
        Check if the start and end of MSA equal to the given patterns
        """
        start_matched, end_matched = self.check_start_end()

        if not start_matched or not end_matched:

            # If exact position matching fails, try sliding window approach
            if not start_matched and self.start_patterns:
                for start_pattern in self.start_patterns:
                    start_window = self.consensus_seq[max(0, self.start - 15): self.start + 15]
                    for i in range(len(start_window) - len(start_pattern) + 1):
                        if start_window[i: i + len(start_pattern)] == list(start_pattern):
                            self.start = max(0, self.start - 15) + i
                            break

            if not end_matched and self.end_patterns:
                for end_pattern in self.end_patterns:
                    end_window = self.consensus_seq[max(0, self.end - 15): self.end + 15]
                    for i in reversed(range(len(end_window) - len(end_pattern) + 1)):
                        if end_window[i: i + len(end_pattern)] == list(end_pattern):
                            self.end = max(0, self.end - 15) + i + len(end_pattern)
                            break

        # Check again after updating positions
        start_matched, end_matched = self.check_start_end()
        return start_matched, end_matched, self.start, self.end
