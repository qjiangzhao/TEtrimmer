from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import os
from functions import select_gaps_block_with_similarity_check


class CleanAndSelectColumn:
    """
    Class to eliminate noise nucleotide and select the diverse columns, which can be used for MSA clustering
    """
    def __init__(self, input_file, threshold=0.05):
        """
        :param input_file: str, absolute path of input file
        :param threshold: nucleotide percentage (gap not count) lower than threshold will be converted to "-"
        """
        self.input_file = input_file
        self.alignment = AlignIO.read(input_file, "fasta")
        self.alignment_seq_num = len(self.alignment)
        self.proportions = {}
        self.threshold = threshold
        self.alignment_length = self.alignment.get_alignment_length()
        self.alignment_filtered = None
        self.if_need_cluster = True
        self.divergent_column_len = None
        self.gap_alignment_filter = None
        self.alignment_filtered_len = None
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

    def select_divergent_column(self, cluster_col_thr=500, dis_col_threshold=0.8):
        """
        Select distinct columns from multiple sequence alignment file, which will be used for clustering
        :param dis_col_threshold: num (0-1) default 0.8, Check if any nucleotide proportion is greater than threshold,
        then delete that column
        :param cluster_col_thr: int default 500, if sequence length longer than cluster_col_thr * 20, the minimum
        column number that required for clustering will be sequence_length * 0.05
        :return: if_need_cluster: boolean, decide if need clustering
        """
        if self.alignment_seq_num >= 90:
            dis_col_threshold = 0.85

        # Delete highly conserved columns
        columns_to_delete = []
        for i in range(self.alignment_length):
            if any(proportion > dis_col_threshold for proportion in self.proportions[i].values()):
                columns_to_delete.append(i)

        columns_to_keep = []
        for i in range(self.alignment_length):
            if i not in columns_to_delete:
                columns_to_keep.append(i)

        # It is only meaningful to calculate the gap block when there are enough distinct columns.
        #min_length_to_include_gap_block = max(25, int(0.1 * self.alignment_length)) \
            #if self.alignment_length <= 1000 else 100

        #if len(columns_to_keep) > min_length_to_include_gap_block:

        gap_block_to_keep = select_gaps_block_with_similarity_check(self.input_file)

        # Concatenate columns with high divergence and gap block columns
        if gap_block_to_keep and (len(columns_to_keep) >= 30):
            divergence_len = len((set(columns_to_keep + gap_block_to_keep)))
        else:
            divergence_len = len(columns_to_keep)
        self.alignment_filtered_len = len(columns_to_keep)
        columns_to_keep = sorted(columns_to_keep)

        """
        When alignment length is greater than 1000, the distinct column number have to be more than 100
        otherwise, it will be ten percent of the MSA length but not less than 50
        """
        min_length = max(50, int(0.05 * self.alignment_length)) if self.alignment_length <= cluster_col_thr / 0.05 \
            else cluster_col_thr

        if divergence_len > min_length:  # Set the threshold number to decide if perform cluster
            self.alignment_filtered = self.alignment[:, columns_to_keep[0]:columns_to_keep[0] + 1]
            for i in columns_to_keep[1:]:
                self.alignment_filtered += self.alignment[:, i:i + 1]
        else:
            self.if_need_cluster = False  # use this values to decide if execute "Class_group_MSA"

        return self.if_need_cluster

    # Remove sequences that contain many gaps. This can cause problem for iqtree clustring.
    def filter_out_big_gap_seq(self, gap_threshold=1):
        gap_alignment_filter_list = []
        for record in self.alignment_filtered:
            gap_count = record.seq.count("-")
            gap_fraction = gap_count / self.alignment_filtered_len

            if gap_fraction < gap_threshold:
                gap_alignment_filter_list.append(record)
        self.gap_alignment_filter = MultipleSeqAlignment(gap_alignment_filter_list)

    def write_alignment_filtered(self, output_dir):
        output_file = os.path.join(output_dir, f"{os.path.basename(self.input_file)}_pat_MSA.fa")
        self.filter_out_big_gap_seq()
        with open(output_file, 'w') as f:
            AlignIO.write(self.gap_alignment_filter, f, 'fasta')
        return output_file