from Bio import AlignIO
import os


class CleanAndSelectColumn:
    """
    Class to eliminate noise nucleotide and select the diverse columns, which can be used for MSA clustering
    """
    def __init__(self, input_file, threshold=0.08):
        """
        :param input_file: str, absolute path of input file
        :param threshold: nucleotide percentage (gap not count) lower than threshold will be converted to "-"
        """
        self.input_file = input_file
        self.alignment = AlignIO.read(input_file, "fasta")
        self.proportions = {}
        self.threshold = threshold
        self.alignment_length = self.alignment.get_alignment_length()
        self.alignment_filtered = None
        self.if_need_cluster = True
        self.calculation_proportion()

    def calculation_proportion(self):
        """
        Method: calculate nucleotide proportion at each column
        """
        sequences = [record.seq for record in self.alignment]

        # Loop through each column of the alignment
        for i in range(self.alignment_length):
            # Count the number of each nucleotide in this column
            counts = {"a": 0, "c": 0, "g": 0, "t": 0}
            for sequence in sequences:
                nucleotide = sequence[i]
                if nucleotide in counts:
                    counts[nucleotide] += 1

            # Calculate the proportion of each nucleotide
            total = sum(counts.values())
            self.proportions[i] = {nucleotide: count / total for nucleotide, count in counts.items()}

    def clean_column(self, output_dir):
        """
        Replace nucleotides that are less than threshold with a gap character
        :return: the absolute path of column clean MSA
        """
        for i in range(self.alignment_length):
            for j, record in enumerate(self.alignment):
                nucleotide = record.seq[i]
                if nucleotide in self.proportions[i] and self.proportions[i][nucleotide] < self.threshold:
                    self.alignment[j].seq = record.seq[:i] + "-" + record.seq[i + 1:]

        output_file = os.path.join(output_dir, f"{os.path.basename(self.input_file)}_cl.fa")
        with open(output_file, 'w') as f:
            AlignIO.write(self.alignment, f, 'fasta')
        return output_file

    def select_divergent_column(self, threshold=0.8):
        """
        Select distinct columns from multiple sequence alignment file, which will be used for clustering
        :param threshold: num (0-1) default 0.8, Check if any nucleotide proportion is greater than threshold,
        then delete that column
        :return: if_need_cluster: boolean, decide if need clustering
        """
        columns_to_delete = []
        for i in range(self.alignment_length):
            if any(proportion > threshold for proportion in self.proportions[i].values()):
                columns_to_delete.append(i)

        columns_to_keep = []
        for i in range(self.alignment_length):
            if i not in columns_to_delete:
                columns_to_keep.append(i)

        """
        When alignment length is greater than 1000, the distinct column number have to be more than 100
        otherwise, it will be ten percent of the MSA length but not less than 50
        """
        min_length = max(50, int(0.1 * self.alignment_length)) if self.alignment_length <= 1000 else 100

        if len(columns_to_keep) > min_length:  # Set the threshold number to decide if perform cluster
            self.alignment_filtered = self.alignment[:, columns_to_keep[0]:columns_to_keep[0] + 1]
            for i in columns_to_keep[1:]:
                self.alignment_filtered += self.alignment[:, i:i + 1]
        else:
            self.if_need_cluster = False  # use this values to decide if execute "Class_group_MSA"

        return self.if_need_cluster

    def write_alignment_filtered(self, output_dir):
        output_file = os.path.join(output_dir, f"{os.path.basename(self.input_file)}_pat_MSA.fa")
        with open(output_file, 'w') as f:
            AlignIO.write(self.alignment_filtered, f, 'fasta')
        return output_file


"""
os.chdir("G:\\TE_manual_curation\\Software_develop\\test_column_cleaning")
test = CleanAndSelectColumn("G:\\TE_manual_curation\\Software_develop\\test_column_cleaning\\LTRR8.fasta")
test2 = test.clean_column()
print(test2)
with open("G:\\TE_manual_curation\\Software_develop\\test_column_cleaning\\LTRR8_column_clean.fasta", "w") as f:
    AlignIO.write(test.alignment, f, "fasta")

test.select_divergent_column()
test.write_alignment_filtered("G:\\TE_manual_curation\\Software_develop\\test_column_cleaning")
"""