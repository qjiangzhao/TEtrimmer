import os.path
import numpy as np
from Bio import AlignIO
from Bio.Seq import Seq
import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


class CropEnd:
    """
    Crop each single sequence end of MSA by the nucleotide divergence.
    """

    def __init__(self, input_file, threshold=16, window_size=20):
        """
        :param input_file: str, path to the multiple sequence alignment
        :param threshold: default 16, nucleotides number inside the check window whose proportion greater than 80%
        :param window_size: default 20, check window size to define start and end position
        """
        self.input_file = input_file
        self.alignment = AlignIO.read(self.input_file, "fasta")
        self.threshold = threshold
        self.window_size = window_size
        # Define a dictionary the key are sequence names, the values are a list contains nucleotides proportion
        self.proportions_dict = {record.id: [] for record in self.alignment}
        # Define a dictionary to hold start and end positions
        self.position_dict = {record.id: [0, 0] for record in self.alignment}
        # Define an empty dataframe to store proportion information
        self.df = None
        self.cropped_alignment = []

    def pro_calculation(self):
        """
        :function pro_calculation: calculate nucleotides proporation at each column per sequence
        :return: a data frame contains all sequence names and nucleotide proporation information
        """
        # Loop through each column of the alignment
        for i in range(self.alignment.get_alignment_length()):
            # Count the number of each nucleotide in this column
            counts = {"a": 0, "c": 0, "g": 0, "t": 0}
            for record in self.alignment:
                nucleotide = record.seq[i]
                if nucleotide in counts:
                    counts[nucleotide] += 1

            # Calculate the proportion of each nucleotide
            total = sum(counts.values())

            # Generate a dictionary named proportions contains nucleotide proportion for this column
            if total < 5:  # Ignore columns less 5 nucleotides
                proportions = {nucleotide: 0 for nucleotide in counts}
            else:
                proportions = {nucleotide: count / total for nucleotide, count in counts.items()}

            # Add the proportion of the nucleotide at this position to each sequence
            for record in self.alignment:  # This will loop each sequences in alignment
                nucleotide = record.seq[i]  # Refer to that column
                if nucleotide in proportions:
                    # Write proportion information into proportions_dict
                    self.proportions_dict[record.id].append(proportions[nucleotide])
                else:
                    # When there is a gap, use number 0 replace proportion
                    self.proportions_dict[record.id].append(np.nan)

            # Convert the dictionary to a DataFrame
        self.df = pd.DataFrame(self.proportions_dict)
        self.df = self.df.transpose()  # transpose the DataFrame so that each row represents a sequence
        self.df.columns = range(1,
                                self.alignment.get_alignment_length() + 1)  # rename the columns to represent positions
        # Convert to two decimal numbers
        self.df = self.df.round(2)

    def average_proportion_per_column(self):
        """
        Calculate the average proportion for each column across all sequences.
        """
        # Ensure the DataFrame is available
        if not hasattr(self, 'df'):
            raise ValueError("The DataFrame has not been created yet. Please run pro_calculation() first.")

        # Calculate the average for each column
        average_proportions = self.df.mean()

        # Calculate the overall mean of the column averages
        overall_average = average_proportions.mean()

        return average_proportions, overall_average

    def find_positions(self):
        """
         this function will define the start and end position for each sequence
            by nucleotide proportions.

        :return: a dictionary that contain sequence name, start, and end positions
        """
        # Loop over the DataFrame's rows
        for index, row in self.df.iterrows():
            # Find start position
            for i in range(len(row) - self.window_size + 1):
                window = row[i:i+self.window_size]
                if window.sum() > self.threshold:
                    self.position_dict[index][0] = i
                    break
            # Find end position
            for i in range(len(row) - 1, self.window_size - 2, -1):
                window = row[i-self.window_size+1:i+1]
                if window.sum() > self.threshold:
                    self.position_dict[index][1] = i+1  # add 1 to make the position 1-indexed
                    break

    def crop_alignment(self):
        # Create a new list to hold the cropped sequences

        # Loop through each sequence in the alignment
        for record in self.alignment:
            # Create a new string with the cropped sequence
            cropped_seq = "-" * self.position_dict[record.id][0] + \
                          str(record.seq[self.position_dict[record.id][0]:self.position_dict[record.id][1]]) + \
                          "-" * (len(record.seq) - self.position_dict[record.id][1])
            # Create a new SeqRecord with the cropped sequence and add it to the list
            self.cropped_alignment.append(SeqRecord(Seq(cropped_seq), id=record.id, description=""))
        # Convert the list of SeqRecords into a MultipleSeqAlignment
        self.cropped_alignment = MultipleSeqAlignment(self.cropped_alignment)

        return self.cropped_alignment

    def write_to_file(self, output_dir):
        output_file = os.path.join(output_dir, f"{os.path.basename(self.input_file)}_ce.fa")
        with open(output_file, "w") as f:
            AlignIO.write(self.cropped_alignment, f, "fasta")
        return output_file


