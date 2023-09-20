import os.path
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

        # Initialize a dictionary to hold the start and end positions of each sequence
        self.position_dict = {record.id: [0, 0] for record in self.alignment}
        self.cropped_alignment = []  # Initialize an empty list to hold the cropped alignment
        self.find_positions()  # Call the method to find the start and end positions of each sequence
        self.crop_alignment()  # Call the method to crop the alignment based on the positions found

    # This method finds the start and end positions of each sequence in the alignment based on the gap proportion
    def find_positions(self):
        for record in self.alignment:  # Loop through each sequence in the alignment
            seq_str = str(record.seq)  # Convert the sequence to a string
            len_seq = len(seq_str)  # Get the length of the sequence

            # Find the start position
            for i in range(len_seq - self.window_size + 1):  # Loop through the sequence with a sliding window
                window = seq_str[i:i+self.window_size]  # Get the subsequence in the window
                gap_proportion = window.count('-') / self.window_size  # Calculate the gap proportion in the window
                # If the gap proportion is less than or equal to the threshold, set this position as the start position
                # and break the loop
                if gap_proportion <= self.gap_threshold:
                    self.position_dict[record.id][0] = i
                    break

            # Find the end position (similar to finding the start position, but from the end of the sequence)
            for i in range(len_seq - 1, self.window_size - 2, -1):
                window = seq_str[i-self.window_size+1:i+1]
                gap_proportion = window.count('-') / self.window_size
                if gap_proportion <= self.gap_threshold:
                    self.position_dict[record.id][1] = i+1  # Note: add 1 to make the position 1-indexed
                    break

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

    # Write the cropped alignment to a new file
    def write_to_file(self, output_dir):
        output_file = os.path.join(output_dir, f"{os.path.basename(self.input_file)}_ceg.fa")
        with open(output_file, "w") as f:
            AlignIO.write(self.cropped_alignment, f, "fasta")
        return output_file


