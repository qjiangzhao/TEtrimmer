import random
import os
import shutil

class BEDFile:
    """
    The bed generated based on blast might contain duplicate elements. This class can remove duplications.
    When bed file is too long, choose desired line number, which can faster multiple sequence alignment.
    """
    def __init__(self, input_file):
        self.input_file = input_file
        self.lines = self.read_bed_file()
        self.remove_duplicates()

    def read_bed_file(self):
        """ Open bed file"""
        with open(self.input_file, 'r') as file:
            lines = file.readlines()
        return [line.strip().split('\t') for line in lines]

    def remove_duplicates(self):
        unique_lines = []
        line_set = set()

        for line in self.lines:
            line_key = tuple(line[:3])
            if line_key not in line_set:
                line_set.add(line_key)
                unique_lines.append(line)

        self.lines = unique_lines

    def get_line_length(self, line):
        return int(line[2]) - int(line[1])

    def select_top_longest_lines(self, n):
        return sorted(self.lines, key=self.get_line_length, reverse=True)[:n]

    def select_random_lines(self, n, remaining_lines):
        random.shuffle(remaining_lines)
        return remaining_lines[:n]

    def process_lines(self, output_dir, threshold=100, top_longest_lines_count=100):
        """
        Select desired line numbers from bed file.

        :param output_dir: str, the absolute path of output folder directory.
        :param threshold: num default=100, the maximum line number for bed file to keep.
        :param top_longest_lines_count: num default=100 smaller or equal to threshold.
        When the bed file line number is excess than the threshold, sort bed file by sequence length and choose
        "top_longest_lines_count" lines. When this number is smaller than the threshold, randomly choose the rest
        number of lines from the bed file.

        :return: the absolute selected bed file path.
        """
        # Test if top_longest_lines_count is equal or smaller than threshold.
        if top_longest_lines_count > threshold:
            raise ValueError("top_mas_lines must be equal to or smaller than max_msa_lines.")

        random_lines_count = threshold - top_longest_lines_count
        bed_out_filter_file = os.path.join(output_dir, f"{os.path.basename(self.input_file)}.fil.bed" )

        if len(self.lines) > threshold:
            top_longest_lines = self.select_top_longest_lines(top_longest_lines_count)
            # Eliminate selected lines
            remaining_lines = [line for line in self.lines if line not in top_longest_lines]
            # Randomly choose lines
            random_lines = self.select_random_lines(random_lines_count, remaining_lines)
            selected_lines = top_longest_lines + random_lines

            # Write the selected lines to the output file.
            with open(bed_out_filter_file, 'w') as file:
                for line in selected_lines:
                    file.write("\t".join(line) + "\n")
        else:
            # Copy the original file to the output directory with the new name.
            # Because I have to keep the file name consistency.
            shutil.copy(self.input_file, bed_out_filter_file)

        return bed_out_filter_file
