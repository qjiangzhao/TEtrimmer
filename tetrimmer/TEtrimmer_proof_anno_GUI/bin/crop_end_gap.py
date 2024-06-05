import subprocess
import sys


def install_and_import(required_packages_dict):
    for package in required_packages_dict:
        try:
            __import__(package)
        except ImportError:
            try:
                print(f"{package} was not installed. Installing it automatically.")
                subprocess.check_call([sys.executable, "-m", "pip", "install", required_packages_dict[package]])
                print(f"{package} was successfully installed.")
            except subprocess.CalledProcessError as e:
                print(
                    f"\nRequired Python packages are missing and cannot be installed automatically. Installation failed with error {e.stderr}"
                    "\nPlease install 'click' and 'biopython' using 'pip install'.\n")
                return


required_packages = {'click': 'click', 'Bio': 'biopython'}
install_and_import(required_packages)


import click
import os
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


class CropEndByGap:

    def __init__(self, input_file, gap_threshold=0.05, window_size=300):
        self.input_file = input_file  # Path to the multiple sequence alignment file
        self.alignment = AlignIO.read(self.input_file, "fasta")  # Read the alignment file in FASTA format
        self.gap_threshold = gap_threshold  # Threshold for gap proportion. Default is 0.05
        self.window_size = window_size  # The size of the window to check for gaps. Default is 50

        # Initialize a dictionary to keep the start and end positions of each sequence
        self.position_dict = {record.id: [0, 0] for record in self.alignment}
        self.cropped_alignment = []  # Initialize an empty list to keep the cropped alignment
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
            if cropped_length > total_length * 0.7:
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
        # Convert the list of SeqRecords into a MultipleSeqAlignment (MSA)
        self.cropped_alignment = MultipleSeqAlignment(self.cropped_alignment)
        return self.cropped_alignment
    
    # Write the cropped alignment to a new file
    def write_to_file(self, output_dir, cropped_alignment):
        output_file = os.path.join(output_dir, f"{os.path.basename(self.input_file)}_ceg.fa")
        with open(output_file, "w") as f:
            AlignIO.write(cropped_alignment, f, "fasta")
        return output_file


def crop_end_gap(input_file, output_file, gap_threshold, window_size):

    if os.path.isfile(input_file):

        crop_end_gap_object = CropEndByGap(input_file, gap_threshold=gap_threshold, window_size=window_size)
        crop_alignment = crop_end_gap_object.crop_alignment()
        with open(output_file, "w") as f:
            AlignIO.write(crop_alignment, f, "fasta")
    else:
        raise FileNotFoundError(f"The file '{input_file}' does not exist.")



@click.command()
@click.option("--input_file", "-i", required="True", type=str,
              help="Multiple sequence alignment FASTA file path")
@click.option("--output_file", "-o", required="True", type=str,
              help="Output file")
@click.option("--gap_threshold", "-thr", default=0.05, type=float,
              help="If the gap proportion in the selected window is greater than --gap_threshold, convert all "
                   "nucleotides to '-' in this window")
@click.option("--window_size", "-ws", default=300, type=int,
              help="Window size used for cropping end")
def crop_end_gap_click(input_file, output_file, gap_threshold, window_size):

    crop_end_gap(input_file, output_file, gap_threshold, window_size)


if __name__ == '__main__':
    crop_end_gap_click()
