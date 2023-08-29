from Bio.Align import AlignInfo, MultipleSeqAlignment
from Bio import AlignIO
import os

class DefineBoundary():

    def __init__(self, input_file, threshold=0.8, check_window=200, max_X=50, if_con_generater=True):
        self.threshold = threshold
        self.input_file = input_file
        self.alignment = None
        self.check_window = check_window
        self.max_X=max_X
        self.if_con_generater=True
        self.consensus_seq = []
        self.nucl = ["A", "G", "C", "T", "a", "g", "c", "t"]
        self.ambiguous = "X"
        self.start_post = None
        self.end_post = None
        self.right_ext = False
        self.left_ext = False
        self.cut_seqs = []
        if if_con_generater:  # Default to use normal consensus sequence generation method
            self.con_generater()
            self.boundary_position()
            self.extension_check()
        else:  # Otherwise, ues consensus generation method that have minimum sequence requirement for each column
            self.con_generater_select_column()
            self.boundary_position()
            self.extension_check()

    def con_generater(self):
        # Read input file
        self.alignment = AlignIO.read(self.input_file, "fasta")
        summary = AlignInfo.SummaryInfo(self.alignment)
        # Get consensus sequence and convert to a list element to enable single element mutation
        self.consensus_seq = list(summary.dumb_consensus(threshold=self.threshold, ambiguous=self.ambiguous).upper())

    # Generate consensus sequences
    def con_generater_select_column(self):
        # Read input file
        self.alignment = AlignIO.read(self.input_file, "fasta")
        summary = AlignInfo.SummaryInfo(self.alignment)
        # Get consensus sequence and convert to a list element to enable single element mutation
        self.consensus_seq = list(summary.dumb_consensus(threshold=self.threshold, ambiguous=self.ambiguous).upper())

        for i in range(len(self.consensus_seq)): # iterate over columns
            column = self.alignment[:, i]
            # set(column) get an unordered collection of unique nucleotide for that column
            nucleotide_counts = {nucleotide: column.count(nucleotide) for nucleotide in set(column) if nucleotide in self.nucl}

            # Check if the column has at least 5 nucleotides and is greater than one fourth of alignemnt sequence number
            if sum(nucleotide_counts.values()) < 5 or sum(nucleotide_counts.values()) <= round(len(self.alignment) / 10):
                self.consensus_seq[i] = self.ambiguous  # if column has fewer than 5 nucleotides, mark as ambiguous
        # Convert list to a string by join function
        self.consensus_seq = "".join(self.consensus_seq)
        return self.consensus_seq

    # Check start and end position, consensus_seq from summary_align.dumb_consensus()
    def boundary_position(self):
        for i, letter in enumerate(self.consensus_seq):
            if letter in self.nucl:
                if i + self.check_window <= len(self.consensus_seq):
                    Xnum = self.consensus_seq[i: i + self.check_window].count(self.ambiguous)
                    # Check if the number of ambiguous smaller than 30% of check_window
                    if Xnum <= self.max_X:
                        self.start_post = i
                        break

        # Check the end position
        for i, letter in reversed(list(enumerate(self.consensus_seq))):
            if letter in self.nucl:
                i = i + 1
                if i - self.check_window >= 0:
                    Xnum = self.consensus_seq[i - self.check_window: i].count(self.ambiguous)
                    # Check if the number of ambiguous smaller than 10% of check_window
                    if Xnum <= self.max_X:
                        self.end_post = i
                        break

    # Check if need more extension for Multiple sequence alignment
    def extension_check(self):
        if self.start_post <= 150:
            # print("Need more extension for left side")
            self.left_ext = True

        if self.end_post >= len(self.consensus_seq) - 150:
            # print("Need more extension for right side")
            self.right_ext = True


    # Iterate through each sequence in the MSA
    def crop_MSA(self, output_dir, crop_extension=0):

        MSA_len = self.alignment.get_alignment_length()

        # Ensure the start position is within the sequence range
        start_col = max(self.start_post - crop_extension, 0)

        # Ensure the end position is within the sequence range
        end_col = min(self.end_post + crop_extension, MSA_len)

        # Select the window columns from the alignment
        selected_alignment = self.alignment[:, start_col:end_col]

        # Create a new MultipleSeqAlignment object with the selected alignment
        selected_alignment = MultipleSeqAlignment(selected_alignment)

        # Create a new MultipleSeqAlignment object with the cut sequences
        selected_alignment = MultipleSeqAlignment(selected_alignment)

        # Write the cut MSA to a file
        output_file = os.path.join(output_dir, f"{os.path.basename(self.input_file)}_bou_crop.fa")

        with open(output_file, "w") as f:
            AlignIO.write(selected_alignment, f, "fasta")
        return output_file


"""
test = DefineBoundary(r"G:\TE_manual_curation\Software_develop\test_crop_end\Test_crop_end.fasta", if_con_generater=False)
print("".join(test.consensus_seq))
"""