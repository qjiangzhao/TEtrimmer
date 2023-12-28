import subprocess
import sys


def install_and_import(required_packages_dic):

    for package in required_packages:
        try:
            __import__(package)
        except ImportError:
            try:
                print(f"{package} was not found. Installing it automatically.")
                subprocess.check_call([sys.executable, "-m", "pip", "install", required_packages_dic[package]])
                print(f"{package} was successfully installed.")
            except subprocess.CalledProcessError as e:
                print(f"\nRequired Python packages are missing and cannot be installed automatically. Installation failed with error {e.stderr}"
                      "\nPlease install 'click', 'numpy', 'pandas', and 'biopython' using 'pip install'.\n")
                return


required_packages = {'click': 'click', 'Bio': 'biopython', 'numpy': 'numpy', 'pandas': 'pandas'}
install_and_import(required_packages)


import click
import os
import numpy as np
from Bio import AlignIO
from Bio.Seq import Seq
import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


class CropEnd:
    """
    Crop each single sequence end of the MSA by nucleotide divergence.
    """

    def __init__(self, input_file, threshold=16, window_size=20):
        """
        :param input_file: str, path to the multiple sequence alignment
        :param threshold: default 16, nucleotide number inside the checking window whose proportion should be greater than 80%
        :param window_size: default 20, checking window size to define start and end position
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
        :function pro_calculation: calculate nucleotide proportion in each column of the MSA
        :return: a data frame containing all sequence names and the respective nucleotide proportion information for each position
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
        This function will define the start and end position for each sequence
            by nucleotide proportions.

        :return: a dictionary that contains sequence name, start, and end positions
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
    


@click.command()
@click.option("--input_file", "-i", required="True", type=str,
              help="Multiple sequence alignment FASTA file path")
@click.option("--output_file", "-o", required="True", type=str,
              help="Output file")
@click.option("--threshold", "-thr", default=0.8, type=float,
              help="Higher number means more stringent cropping. Range: 0-1")
@click.option("--window_size", "-ws", default=20, type=int,
              help="Window size used for cropping end")
def crop_end_div(input_file, output_file, threshold, window_size):

    if os.path.isfile(input_file):

        crop_end_thr = threshold * window_size
        crop_end_div_object = CropEnd(input_file, threshold=crop_end_thr, window_size=window_size)
        crop_alignment = crop_end_div_object.crop_alignment()
        with open(output_file, "w") as f:
            AlignIO.write(crop_alignment, f, "fasta")
    else:
        raise FileNotFoundError(f"The file '{input_file}' does not exist.")


if __name__ == '__main__':
    crop_end_div()
