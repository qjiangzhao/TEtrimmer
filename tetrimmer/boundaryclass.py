import os
import warnings

import numpy as np
import pandas as pd
from Bio import AlignIO, BiopythonDeprecationWarning
from Bio.Align import AlignInfo, MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Suppress all deprecation warnings
warnings.filterwarnings('ignore', category=BiopythonDeprecationWarning)


class CropEnd:
    """
    Crop each single sequence end of the MSA by nucleotide divergence.
    """

    def __init__(
        self, input_file, threshold=16, window_size=20, crop_l=True, crop_r=True
    ):
        """
        :param input_file: str, path to the multiple sequence alignment
        :param threshold: int, number of nucleotide sites inside the checking window with a proportion greater than 80%. Default: 16
        :param window_size: int, check window size to define start and end position. Default: 20
        :param crop_l and crop_r: decide if cropping one or both ends
        """
        self.input_file = input_file
        self.alignment = AlignIO.read(self.input_file, 'fasta')
        self.threshold = threshold
        self.window_size = window_size
        self.crop_l = crop_l
        self.crop_r = crop_r
        self.alignment_len = self.alignment.get_alignment_length()
        # Define a dictionary; keys are sequence names and values are a list containing nucleotide proportions
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
        :function pro_calculation: calculate nucleotide proportions for each column and sequence
        :return: a dataframe containing all sequence names and nucleotide proportion information
        """
        # Loop through each column of the alignment
        for i in range(self.alignment_len):
            # Count the number of each nucleotide in this column
            counts = {'a': 0, 'c': 0, 'g': 0, 't': 0}
            for record in self.alignment:
                nucleotide = record.seq[i].lower()
                if nucleotide in counts:
                    counts[nucleotide] += 1

            # Calculate the proportion of each nucleotide for this column
            total = sum(counts.values())

            # Generate a dictionary named proportions containing nucleotide proportions for this column
            if total < 5:  # Ignore columns with fewer than 5 nucleotides
                proportions = {nucleotide: 0 for nucleotide in counts}
            else:
                proportions = {
                    nucleotide: count / total for nucleotide, count in counts.items()
                }

            # Add the proportion of each nucleotide at this position to each sequence
            for record in (
                self.alignment
            ):  # This will loop through all sequences in the alignment
                nucleotide = record.seq[i]  # Refer to that column
                if nucleotide in proportions:
                    # Write proportion information into proportions_dict
                    self.proportions_dict[record.id].append(proportions[nucleotide])
                else:
                    # If there is a gap, use number 0 to replace proportion
                    self.proportions_dict[record.id].append(np.nan)

        # Convert the dictionary into a DataFrame
        self.df = pd.DataFrame(self.proportions_dict)
        self.df = (
            self.df.transpose()
        )  # transpose the DataFrame so that each row represents a sequence
        self.df.columns = range(
            1, self.alignment_len + 1
        )  # rename the columns to represent positions
        # Convert to two decimal numbers
        self.df = self.df.round(2)

    def find_positions(self):
        """
         This function will define the start and end position for each sequence
            by nucleotide proportions.

        :return: a dictionary that contains sequence name, start, and end positions
        """
        # Loop over the DataFrame rows
        for index, row in self.df.iterrows():
            if self.crop_l:
                # Find start position
                for i in range(len(row) - self.window_size + 1):
                    window = row[i : i + self.window_size]
                    if window.sum() > self.threshold:
                        self.position_dict[index][0] = i
                        break
            else:
                # Set the start position to 0 if crop_l is not used
                self.position_dict[index][0] = 0

            if self.crop_r:
                # Find end position
                for i in range(len(row) - 1, self.window_size - 1, -1):
                    window = row[i - self.window_size : i + 1]
                    if window.sum() > self.threshold:
                        self.position_dict[index][1] = (
                            i + 1
                        )  # add 1 to make the position 1-indexed
                        break
            else:
                # Set the end position to the end of the alignment if crop_r is not used
                self.position_dict[index][1] = self.alignment_len

    def crop_alignment(self):
        # Create a new list to hold the cropped sequences
        # Loop through each sequence in the alignment
        for record in self.alignment:
            # Create a new string with the cropped sequence
            cropped_seq = (
                '-' * self.position_dict[record.id][0]
                + str(
                    record.seq[
                        self.position_dict[record.id][0] : self.position_dict[
                            record.id
                        ][1]
                    ]
                )
                + '-' * (len(record.seq) - self.position_dict[record.id][1])
            )
            # Create a new SeqRecord with the cropped sequence and add it to the list
            self.cropped_alignment.append(
                SeqRecord(Seq(cropped_seq), id=record.id, description='')
            )
        # Convert the list of SeqRecords into a MultipleSeqAlignment
        self.cropped_alignment = MultipleSeqAlignment(self.cropped_alignment)

        return self.cropped_alignment

    def average_proportion_per_column(self):
        """
        Calculate the average proportion for each column across all sequences.
        """
        # Ensure the DataFrame is available
        if not hasattr(self, 'df'):
            raise ValueError(
                'The DataFrame has not been created yet. Please run pro_calculation() first.'
            )

        # Calculate the average for each column
        average_proportions = self.df.mean()

        # Calculate the overall mean of the column averages
        overall_average = average_proportions.mean()

        return average_proportions, overall_average

    def write_to_file(self, output_dir):
        # Define different names for the different directions of cropping.
        if self.crop_l and self.crop_r:
            output_file = os.path.join(
                output_dir, f'{os.path.basename(self.input_file)}_ce.fa'
            )
        elif self.crop_l and not self.crop_r:
            output_file = os.path.join(
                output_dir, f'{os.path.basename(self.input_file)}_cel.fa'
            )
        elif not self.crop_l and self.crop_r:
            output_file = os.path.join(
                output_dir, f'{os.path.basename(self.input_file)}_cer.fa'
            )

        with open(output_file, 'w') as f:
            AlignIO.write(self.cropped_alignment, f, 'fasta')
        return output_file


class CropEndByGap:
    def __init__(self, input_file, gap_threshold=0.05, window_size=300):
        self.input_file = input_file  # Path to the multiple sequence alignment file
        self.alignment = AlignIO.read(
            self.input_file, 'fasta'
        )  # Read the alignment file in FASTA format
        self.gap_threshold = (
            gap_threshold  # Threshold for gap proportion. Default is 0.05
        )
        self.window_size = (
            window_size  # The size of the window to check for gaps. Default is 50
        )

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
            for i in range(
                len_seq - self.window_size + 1
            ):  # Loop through the sequence using a sliding window
                window = seq_str[
                    i : i + self.window_size
                ]  # Get the subsequence in the window
                gap_proportion = (
                    window.count('-') / self.window_size
                )  # Calculate the gap proportion in the window

                # If the gap proportion is less than or equal to the threshold, set this position as the start position
                # and break the loop
                if gap_proportion <= self.gap_threshold:
                    self.position_dict[record.id][0] = i
                    break

            # Find the end position (similar to finding the start position, but from the end of the sequence)
            for i in range(len_seq - 1, self.window_size - 2, -1):
                window = seq_str[i - self.window_size + 1 : i + 1]
                gap_proportion = window.count('-') / self.window_size
                if gap_proportion <= self.gap_threshold:
                    end_position = i + 1  # Note: add 1 to make the position 1-indexed

                    # Ensure end position is not smaller than the start position
                    if end_position <= self.position_dict[record.id][0]:
                        self.position_dict[record.id][1] = self.position_dict[
                            record.id
                        ][0]
                    else:
                        self.position_dict[record.id][1] = end_position

                    break

    # New method to find sequences with cropped region > 90% of the alignment length
    def find_large_crops(self):
        large_crop_ids = []  # List to store the IDs of sequences with large cropped regions
        remaining_sequence_ids = []  # List to store the IDs of the remaining sequences
        total_length = (
            self.alignment.get_alignment_length()
        )  # Total length of the alignment

        for seq_id, positions in self.position_dict.items():
            # Remove '-' and '+' from the sequence ID
            seq_id = seq_id.replace('(-)', '').replace('(+)', '')
            start, end = positions
            cropped_length = start + (total_length - end)

            # Check if the cropped region is greater than 90% of the total alignment length
            if cropped_length > total_length * 0.9:
                large_crop_ids.append(seq_id)  # Add the sequence ID to the list
            else:
                remaining_sequence_ids.append(
                    seq_id
                )  # Add the sequence ID to the remaining sequences list

        # The number of sequences that remain can be calculated as the length of remaining_sequence_ids
        remaining_sequences = len(remaining_sequence_ids)
        large_crop_ids = [int(seq_id) for seq_id in large_crop_ids]
        remaining_sequence_ids = [int(seq_id) for seq_id in remaining_sequence_ids]

        return large_crop_ids, remaining_sequences, remaining_sequence_ids

    # This method crops the alignment based on the start and end positions found
    def crop_alignment(self):
        for record in self.alignment:  # Loop through each sequence in the alignment
            # Create a new sequence by cropping the original sequence based on the start and end positions
            cropped_seq = (
                '-' * self.position_dict[record.id][0]
                + str(
                    record.seq[
                        self.position_dict[record.id][0] : self.position_dict[
                            record.id
                        ][1]
                    ]
                )
                + '-' * (len(record.seq) - self.position_dict[record.id][1])
            )

            # Create a new SeqRecord with the cropped sequence and add it to the list
            self.cropped_alignment.append(
                SeqRecord(Seq(cropped_seq), id=record.id, description='')
            )
        # Convert the list of SeqRecords into a MultipleSeqAlignment
        self.cropped_alignment = MultipleSeqAlignment(self.cropped_alignment)
        return self.cropped_alignment

    # Write the cropped alignment into a new file
    def write_to_file(self, output_dir, cropped_alignment):
        output_file = os.path.join(
            output_dir, f'{os.path.basename(self.input_file)}_ceg.fa'
        )
        with open(output_file, 'w') as f:
            AlignIO.write(cropped_alignment, f, 'fasta')
        return output_file


class DefineBoundary:
    def __init__(
        self,
        input_file,
        threshold=0.8,
        check_window=200,
        max_X=0.25,
        if_con_generater=True,
        extension_buffer=150,
        end_position=None,
    ):
        self.threshold = threshold
        self.input_file = input_file
        self.alignment = None
        self.check_window = check_window
        self.max_X = max_X
        self.if_con_generater = True
        self.consensus_seq = []
        self.nucl = ['A', 'G', 'C', 'T', 'a', 'g', 'c', 't']
        self.ambiguous = 'N'
        self.start_post = None
        self.end_post = None
        self.right_ext = False
        self.left_ext = False
        self.if_continue = True
        self.cut_seqs = []
        self.extension_buffer = extension_buffer
        self.end_position = end_position
        if (
            if_con_generater
        ):  # Default to use standard consensus sequence generation method
            self.con_generator()
            self.boundary_position()
            self.extension_check()
        else:  # Otherwise, use consensus generation method with minimum sequence requirements for each column
            self.con_generator_select_column()
            self.boundary_position()
            self.extension_check()

    def con_generator(self):
        # Read input file
        self.alignment = AlignIO.read(self.input_file, 'fasta')
        summary = AlignInfo.SummaryInfo(self.alignment)
        # Get consensus sequence and convert to a list element to enable single element mutation
        self.consensus_seq = list(
            summary.dumb_consensus(
                threshold=self.threshold, ambiguous=self.ambiguous
            ).upper()
        )

    # Generate consensus sequences
    def con_generator_select_column(self):
        """
        Generate consensus sequence if the column has more than 5 nucleotides. For columns with less than 5 nucleotides,
        write letter indicating ambiguity.
        """
        # Read input file
        self.alignment = AlignIO.read(self.input_file, 'fasta')
        summary = AlignInfo.SummaryInfo(self.alignment)
        # Get consensus sequence and convert to a list element to enable single element mutation
        self.consensus_seq = list(
            summary.dumb_consensus(
                threshold=self.threshold, ambiguous=self.ambiguous
            ).upper()
        )

        for i in range(len(self.consensus_seq)):  # iterate over columns
            column = self.alignment[:, i]
            # set(column) returns an unordered collection of unique nucleotides in that column
            nucleotide_counts = {
                nucleotide: column.count(nucleotide)
                for nucleotide in set(column)
                if nucleotide in self.nucl
            }

            # Check if the column has at least 5 nucleotides and is greater than one-tenth of alignment sequence number
            # if sum(nucleotide_counts.values()) < 5 or
            # sum(nucleotide_counts.values()) <= round(len(self.alignment) / 10):
            if sum(nucleotide_counts.values()) < 5:
                self.consensus_seq[i] = (
                    self.ambiguous
                )  # if column has fewer than 5 nucleotides, mark as ambiguous
        # Convert list to a string using the join function
        self.consensus_seq = ''.join(self.consensus_seq)
        return self.consensus_seq

    # Check start and end position, consensus_seq from summary_align.dumb_consensus()
    def boundary_position(self):
        # Check the start position
        for i, letter in enumerate(self.consensus_seq):
            if letter in self.nucl:
                if i + self.check_window <= len(self.consensus_seq):
                    Xnum = self.consensus_seq[i : i + self.check_window].count(
                        self.ambiguous
                    )
                    # Check if the number of ambiguous sites is smaller than 30% of check_window
                    if Xnum <= (self.max_X * self.check_window):
                        self.start_post = i
                        break

        # Set start_post to the sequence length if it is not defined
        if self.start_post is None:
            self.start_post = len(self.consensus_seq)

        if self.end_position is not None:
            self.end_post = self.end_position
        else:
            # Check the end position
            for i, letter in reversed(list(enumerate(self.consensus_seq))):
                if letter in self.nucl:
                    i = i + 1
                    if i - self.check_window >= 0:
                        Xnum = self.consensus_seq[i - self.check_window : i].count(
                            self.ambiguous
                        )
                        # Check if the number of ambiguous sites is smaller than 10% of check_window
                        if Xnum <= (self.max_X * self.check_window):
                            self.end_post = i
                            break

        # Set end_post to 1 if it is not defined
        if self.end_post is None:
            self.end_post = 1

    # Check if MSA needs further extension
    def extension_check(self):
        # Check if the start position is smaller than the end position.
        if self.start_post < self.end_post:
            if self.start_post <= self.extension_buffer:
                # print("Needs additional extension on the left side of the MSA")
                self.left_ext = True

            if self.end_post >= len(self.consensus_seq) - self.extension_buffer:
                # print("Needs additional extension on the right side of the MSA")
                self.right_ext = True
        else:
            # If if_continue is 'False', the MSA is too short
            self.if_continue = False

    # Iterate through each sequence in the MSA
    def crop_MSA(self, output_dir, crop_extension=0):
        MSA_len = self.alignment.get_alignment_length()

        # Ensure the start position is within the sequence range
        start_col = max(self.start_post - crop_extension, 0)
        # Ensure the end position is within the sequence range
        end_col = min(self.end_post + crop_extension, MSA_len)
        # Select the window columns from the alignment
        selected_alignment = self.alignment[:, start_col:end_col]
        # Create a new MultipleSeqAlignment object using the selected alignment
        selected_alignment = MultipleSeqAlignment(selected_alignment)
        # Write the cut MSA to a file
        output_file = os.path.join(
            output_dir, f'{os.path.basename(self.input_file)}_bc.fa'
        )

        with open(output_file, 'w') as f:
            AlignIO.write(selected_alignment, f, 'fasta')
        return output_file
