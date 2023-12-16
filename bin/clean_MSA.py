import numpy as np
import pandas as pd
import os
from Bio import AlignIO
from Bio.Seq import Seq
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
        self.pro_calculation()
        self.find_positions()

    def pro_calculation(self):
        """
        :function pro_calculation: calculate nucleotides proportion at each column per sequence
        :return: a data frame contains all sequence names and nucleotide proportion information
        """
        # Loop through each column of the alignment
        for i in range(self.alignment.get_alignment_length()):
            # Count the number of each nucleotide in this column
            counts = {"a": 0, "c": 0, "g": 0, "t": 0}
            for record in self.alignment:
                nucleotide = record.seq[i].lower()
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
                window = row[i:i + self.window_size]
                if window.sum() > self.threshold:
                    self.position_dict[index][0] = i
                    break
            # Find end position
            for i in range(len(row) - 1, self.window_size - 2, -1):
                window = row[i - self.window_size + 1:i + 1]
                if window.sum() > self.threshold:
                    self.position_dict[index][1] = i + 1  # add 1 to make the position 1-indexed
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

    def write_to_file(self, output_dir):
        output_file = os.path.join(output_dir, f"{os.path.basename(self.input_file)}_ce.fa")
        with open(output_file, "w") as f:
            AlignIO.write(self.cropped_alignment, f, "fasta")
        return output_file


class CropEndByGap:

    def __init__(self, input_file, gap_threshold=0.05, window_size=300):
        self.input_file = input_file  # Path to the multiple sequence alignment file
        self.alignment = AlignIO.read(self.input_file, "fasta")  # Read the alignment file in FASTA format
        self.gap_threshold = gap_threshold  # Threshold for gap proportion. Default is 0.05
        self.window_size = window_size  # The size of the window to check for gaps. Default is 50

        # Initialize a dictionary to hold the start and end positions of each sequence
        self.position_dict = {record.id: [0, 0] for record in self.alignment}
        self.cropped_alignment = []  # Initialize an empty list to hold the cropped alignment
        self.find_positions()  # Call the method to find the start and end positions of each sequence

    # This method finds the start and end positions of each sequence in the alignment based on the gap proportion
    def find_positions(self):
        for record in self.alignment:  # Loop through each sequence in the alignment
            seq_str = str(record.seq)  # Convert the sequence to a string
            len_seq = len(seq_str)  # Get the length of the sequence

            # Find the start position
            for i in range(len_seq - self.window_size + 1):  # Loop through the sequence with a sliding window
                window = seq_str[i:i + self.window_size]  # Get the subsequence in the window
                gap_proportion = window.count('-') / self.window_size  # Calculate the gap proportion in the window

                # If the gap proportion is less than or equal to the threshold, set this position as the start position
                # and break the loop
                if gap_proportion <= self.gap_threshold:
                    self.position_dict[record.id][0] = i
                    break

            # Find the end position (similar to finding the start position, but from the end of the sequence)
            for i in range(len_seq - 1, self.window_size - 2, -1):
                window = seq_str[i - self.window_size + 1:i + 1]
                gap_proportion = window.count('-') / self.window_size
                if gap_proportion <= self.gap_threshold:
                    end_position = i + 1  # Note: add 1 to make the position 1-indexed

                    # Ensure end position is not smaller than the start position
                    if end_position <= self.position_dict[record.id][0]:
                        self.position_dict[record.id][1] = self.position_dict[record.id][0]
                    else:
                        self.position_dict[record.id][1] = end_position

                    break

    # New method to find sequences with cropped region > 50% of the alignment length
    def find_large_crops(self):
        large_crop_ids = []  # List to store the IDs of sequences with large cropped regions
        remaining_sequence_ids = []  # List to store the IDs of the remaining sequences
        total_length = self.alignment.get_alignment_length()  # Total length of the alignment

        for seq_id, positions in self.position_dict.items():
            # Remove '-' and '+' from the sequence ID
            seq_id = seq_id.replace('(-)', '').replace('(+)', '')
            start, end = positions
            cropped_length = start + (total_length - end)

            # Check if the cropped region is greater than 50% of the total alignment length
            if cropped_length > total_length * 0.9:
                large_crop_ids.append(seq_id)  # Add the sequence ID to the list
            else:
                remaining_sequence_ids.append(seq_id)  # Add the sequence ID to the remaining sequences list

        # The number of sequences remaining can be calculated as the length of remaining_sequence_ids
        remaining_sequences = len(remaining_sequence_ids)
        large_crop_ids = [int(seq_id) for seq_id in large_crop_ids]
        remaining_sequence_ids = [int(seq_id) for seq_id in remaining_sequence_ids]

        return large_crop_ids, remaining_sequences, remaining_sequence_ids

    # This method crops the alignment based on the start and end positions found
    def crop_alignment(self):
        for record in self.alignment:  # Loop through each sequence in the alignment

            # Create a new sequence by cropping the original sequence based on the start and end positions
            cropped_seq = "-" * self.position_dict[record.id][0] + \
                          str(record.seq[self.position_dict[record.id][0]:self.position_dict[record.id][1]]) + \
                          "-" * (len(record.seq) - self.position_dict[record.id][1])

            # Create a new SeqRecord with the cropped sequence and add it to the list
            self.cropped_alignment.append(SeqRecord(Seq(cropped_seq), id=record.id, description=""))
        # Convert the list of SeqRecords into a MultipleSeqAlignment
        self.cropped_alignment = MultipleSeqAlignment(self.cropped_alignment)
        return self.cropped_alignment

    # Write the cropped alignment to a new file
    def write_to_file(self, output_dir, cropped_alignment):
        output_file = os.path.join(output_dir, f"{os.path.basename(self.input_file)}_ceg.fa")
        with open(output_file, "w") as f:
            AlignIO.write(cropped_alignment, f, "fasta")
        return output_file
