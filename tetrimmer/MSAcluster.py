import logging
import os
import os.path
import re
import subprocess
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import AlignIO, Phylo
from matplotlib.colors import ListedColormap
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA

from functions import (
    align_sequences,
    filter_out_big_gap_seq,
    muscle_align,
    famsa_align,
    remove_gaps_with_similarity_check,
    select_gaps_block_with_similarity_check,
)

try:
    pd.set_option('future.no_silent_downcasting', True)
except Exception:
    pass


class CleanAndSelectColumn:
    """
    Class to eliminate noise nucleotides and select diverse columns, which can be used for MSA clustering
    """

    def __init__(self, input_file, threshold=0.05):
        """
        :param input_file: str, absolute path of input file
        :param threshold: nucleotide percentage (excluding gaps) lower than threshold will be converted to "-"
        """
        self.input_file = input_file
        self.alignment = AlignIO.read(input_file, 'fasta')
        self.alignment_seq_num = len(self.alignment)
        self.proportions = {}
        self.threshold = threshold
        self.alignment_length = self.alignment.get_alignment_length()
        self.alignment_filtered = None
        self.if_need_cluster = True
        self.divergent_column_len = None
        self.alignment_filtered_len = None
        self.calculation_proportion()

    def calculation_proportion(self):
        """
        Method: calculate nucleotide proportion in each column
        """
        sequences = [record.seq for record in self.alignment]

        # Loop through each column of the alignment
        for i in range(self.alignment_length):
            # Count the number of each nucleotide in this column
            counts = {'a': 0, 'c': 0, 'g': 0, 't': 0}
            for sequence in sequences:
                nucleotide = sequence[i]
                if nucleotide in counts:
                    counts[nucleotide] += 1

            # Calculate the proportion of each nucleotide
            total = sum(counts.values())
            self.proportions[i] = {
                nucleotide: count / total for nucleotide, count in counts.items()
            }

    def clean_column(self, output_dir):
        """
        Replace nucleotides with proportions below threshold with a gap character
        :return: the absolute path of column clean MSA
        """
        for i in range(self.alignment_length):
            for j, record in enumerate(self.alignment):
                nucleotide = record.seq[i]
                if (
                    nucleotide in self.proportions[i]
                    and self.proportions[i][nucleotide] < self.threshold
                ):
                    self.alignment[j].seq = record.seq[:i] + '-' + record.seq[i + 1 :]

        output_file = os.path.join(
            output_dir, f'{os.path.basename(self.input_file)}_cl.fa'
        )
        with open(output_file, 'w') as f:
            AlignIO.write(self.alignment, f, 'fasta')
        return output_file

    def select_divergent_column(self, cluster_col_thr=500, dis_col_threshold=0.8):
        """
        Select distinct columns from the multiple sequence alignment file, which will be used for clustering
        :param dis_col_threshold: float, check if any nucleotide proportion is greater than threshold, if true,
        delete that column. Default: 0.8
        :param cluster_col_thr: int, if sequence length is longer than cluster_col_thr * 20, the minimum
        column number required for clustering will be sequence_length * 0.05. Default: 500
        :return: if_need_cluster: boolean, decide if clustering is required

        """
        if self.alignment_seq_num >= 90:
            dis_col_threshold = 0.85

        # Delete highly conserved columns
        columns_to_delete = []
        for i in range(self.alignment_length):
            if any(
                proportion > dis_col_threshold
                for proportion in self.proportions[i].values()
            ):
                columns_to_delete.append(i)

        columns_to_keep = []
        for i in range(self.alignment_length):
            if i not in columns_to_delete:
                columns_to_keep.append(i)

        # It is only meaningful to calculate the gap block if there are enough distinct columns.
        # min_length_to_include_gap_block = max(25, int(0.1 * self.alignment_length)) \
        # if self.alignment_length <= 1000 else 100

        # if len(columns_to_keep) > min_length_to_include_gap_block:

        gap_block_to_keep = select_gaps_block_with_similarity_check(self.input_file)

        # Concatenate columns with high divergence and gap block columns
        if gap_block_to_keep and (len(columns_to_keep) >= 30):
            divergence_len = len((set(columns_to_keep + gap_block_to_keep)))
        else:
            divergence_len = len(columns_to_keep)
        self.alignment_filtered_len = len(columns_to_keep)
        columns_to_keep = sorted(columns_to_keep)

        """
        If alignment length is greater than 2000, the number of distinct columns has to be greater than 100.
        Otherwise, the number will be ten percent of the MSA length but no less than 50 columns.
        """
        min_length = (
            max(50, int(0.05 * self.alignment_length))
            if self.alignment_length <= cluster_col_thr / 0.05
            else cluster_col_thr
        )

        if (
            divergence_len > min_length
        ):  # Set the threshold number to decide if performing clustering
            self.alignment_filtered = self.alignment[
                :, columns_to_keep[0] : columns_to_keep[0] + 1
            ]
            for i in columns_to_keep[1:]:
                self.alignment_filtered += self.alignment[:, i : i + 1]
        else:
            self.if_need_cluster = (
                False  # Use this value to decide if executing "Class_group_MSA"
            )

        return self.if_need_cluster

    def write_alignment_filtered(self, output_dir):
        output_file = os.path.join(
            output_dir, f'{os.path.basename(self.input_file)}_pat_MSA.fa'
        )

        # Remove sequences that contain many gaps. Gaps can cause problems for IQ-TREE clustering.
        gap_alignment_filter = filter_out_big_gap_seq(
            self.alignment_filtered, gap_threshold=0.9
        )
        with open(output_file, 'w') as f:
            AlignIO.write(gap_alignment_filter, f, 'fasta')
        return output_file


def calculate_tree_dis(input_file):
    # Read the tree from a Newick file
    tree = Phylo.read(input_file, 'newick')
    output_file = f'{input_file}.tdistance'
    # Get all terminals (leaves) in the tree
    terminals = tree.get_terminals()
    num_terminals = len(terminals)

    # Initialize an empty distance matrix with sequence names
    distance_matrix = [[''] * (num_terminals + 1) for _ in range(num_terminals + 1)]

    # Set the first row and column as sequence names
    for i in range(1, num_terminals + 1):
        distance_matrix[i][0] = terminals[i - 1].name
        distance_matrix[0][i] = terminals[i - 1].name

    # Calculate pairwise distances
    for i in range(1, num_terminals + 1):
        for j in range(1, num_terminals + 1):
            # If i and j are equal, set distance to 1
            if i == j:
                distance_matrix[i][j] = 0
            else:
                # Calculate the total branch length between terminal nodes i and j
                total_branch_length = tree.distance(terminals[i - 1], terminals[j - 1])

                # Assign the calculated distance to the matrix
                distance_matrix[i][j] = total_branch_length

    # Write the distance matrix to a file
    with open(output_file, 'w') as f:
        for row in distance_matrix:
            f.write('\t'.join(map(str, row)) + '\n')
    return output_file


def read_msa(input_file):
    """Read the MSA file"""
    return AlignIO.read(input_file, 'fasta')


def cluster_msa_iqtree_DBSCAN(alignment, min_cluster_size=10, max_cluster=2):
    # this function uses iqtree to separate alignment into branches and then do dbscan_cluster on sequences
    # based on maximum likelihood tree distance
    # dependencies: calculate_tree_dis, dbscan_cluster
    # temp using K2P+I
    iqtree_command = [
        'iqtree',
        '-m',
        'K2P+I',
        '--redo-tree',
        '--seqtype',
        'DNA',
        '-s',
        alignment,
    ]
    try:
        subprocess.run(
            iqtree_command,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

    except FileNotFoundError as e:
        logging.error(
            "'iqtree' command not found. Please ensure 'iqtree' is correctly installed."
        )
        raise Exception from e

    except subprocess.CalledProcessError as e:
        logging.error(f'\niqtree failed with error code {e.returncode}\n{e.stdout}\n{e.stderr}\n')
        raise Exception('iqtree error') from e

    treefile = f'{alignment}.treefile'
    distance_file = calculate_tree_dis(treefile)
    cluster, sequence_names = dbscan_cluster(distance_file, pca=False)

    # map sequence_names name to cluster
    sequence_cluster_mapping = dict(zip(sequence_names, cluster))

    """
    # if all cluster size < 10, there is no meaningful cluster, if_cluster = False, skip this sequence
    # When no -1 cluster number is greater or equal to max_cluster, use no -1 cluster for further analysis
    # otherwise, check if -1 cluster have more than or equal 15 sequence and the sequence number more than 60% of the
    # total sequence, if so count -1 cluster into the further analysis process.
    # Use Counter to count occurrences
    -1 cluster has to be treated carefully.
    """
    counter = Counter(cluster)

    # Find cluster with size > min_cluster_size
    filter_cluster = [
        element for element, count in counter.items() if count >= min_cluster_size
    ]
    seq_cluster_list = []

    # pre-define if_cluster to true
    if_cluster = True

    # Return empty list and False directly when filter_cluster is empty. This  means no cluster have sequence number
    # greater than min_cluster_size
    if not filter_cluster:
        return [], False

    negative_n = counter.get(-1, 0)

    # Filter and sort clusters excluding -1
    top_cluster_obj = Counter(
        {
            element: count
            for element, count in counter.items()
            if element in filter_cluster and element != -1
        }
    ).most_common
    if negative_n >= 15 and negative_n >= 0.6 * len(sequence_names):
        # top_cluster_obj() returns all elements sorted by frequency
        if len(top_cluster_obj()) >= max_cluster:
            # top_cluster is assigned to be a list
            top_cluster = top_cluster_obj(max_cluster)
        else:
            # Add -1 cluster when it has more than 15 sequences and the total sequence inside this cluster occupy
            # more than 60% of the total sequence
            top_cluster = top_cluster_obj(max_cluster)
            top_cluster.append((-1, negative_n))
    else:
        # Pick the top max_cluster
        top_cluster = top_cluster_obj(max_cluster)
        # Find the seq name in selected cluster

    if len(top_cluster) > 0:
        for i in top_cluster:
            # Extract sequence id for each cluster
            # i[0] represents the cluster number
            seq_records = [
                sequence
                for sequence, cluster_label in sequence_cluster_mapping.items()
                if cluster_label == i[0]
            ]
            seq_cluster_list.append(seq_records)
    else:
        # len(top_cluster) == 0 means no cluster fit the requirement, set if_cluster to False
        if_cluster = False

    return seq_cluster_list, if_cluster


def dbscan_cluster(input_file, pca=True):
    # Read the distance matrix from the file
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # Extract distance values
    distance_values = []
    sequence_names = []  # Add this line to store sequence names
    for line in lines[1:]:
        values = line.strip().split('\t')
        sequence_names.append(
            values[0]
        )  # Assuming the sequence name is in the first column
        float_values = [
            float(value) for value in values[1:]
        ]  # Adjust the index based on your data
        distance_values.append(float_values)

    # Convert distance values to a NumPy array
    distance_np = np.array(distance_values, dtype=np.float64)

    # Replace NaN values with a large number, NaN value can hamper DBSCAN clustering
    distance_np[np.isnan(distance_np)] = np.nanmax(distance_np) + 1

    # Apply DBSCAN clustering with the recommended eps value
    dbscan = DBSCAN(eps=0.1, min_samples=3, metric='precomputed')
    labels = dbscan.fit_predict(distance_np)
    if pca:
        title = os.path.basename(input_file)
        plot_pca(distance_np, labels, input_file, title=title)

    return labels, sequence_names


def plot_pca(distance_matrix, labels, output_file, title='PCA Plot'):
    # Calculate PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(distance_matrix)

    # Create a DataFrame with the results
    pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])
    pca_df['Cluster'] = labels

    # Plot the PCA with colored clusters
    plt.figure(figsize=(10, 6))
    cluster_counts = {}  # Dictionary to hold the count of points in each cluster

    for label in set(labels):
        indices = np.where(labels == label)[0]
        cluster_points = pca_df.loc[indices, ['PC1', 'PC2']]
        plt.scatter(
            cluster_points['PC1'], cluster_points['PC2'], label=f'Cluster {label}'
        )

        # Calculate the count of points in the cluster
        cluster_counts[label] = len(indices)

        # Annotate the cluster with its point count
        # The annotation is placed at the mean position of the cluster's points
        plt.annotate(
            f'Count: {len(indices)}',
            (cluster_points['PC1'].mean(), cluster_points['PC2'].mean()),
            textcoords='offset points',
            xytext=(0, 10),
            ha='center',
        )

    plt.title(f'{title} - PCA Plot')
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.legend()

    # Save the plot to a file
    output_file_n = f'{output_file}_PCA.pdf'
    plt.savefig(output_file_n)
    plt.close()


def process_labels(input_file, filtered_cluster_records):
    """
    During the IQtree process, the sequence name is modified to remove special characters.
    This step extracts BED info based on the index of sequence name in the BEDfile created by nameOnly
    :param input_file: str, the absolute path of the BED file
    :return: a list containing the clustered BED files
    """
    bed_dfs = []
    bed_df = pd.read_csv(input_file, sep='\t', header=None)
    for cluster in filtered_cluster_records:
        # Use regular expression to find all numbers
        ids = [re.findall(r'\d+', seq_id)[0] for seq_id in cluster]
        # Convert the numbers to integers
        ids = [int(num) for num in ids]
        # select rows where the fourth column is in filter_cluster_records
        cluster_df = bed_df[bed_df[3].isin(ids)]
        bed_dfs.append(cluster_df)

    return bed_dfs


def subset_alignment_dis(filtered_cluster_records, output_dir, input_align_dis_path):
    """
    This function will subset the input file means the pattern MSA.
    """
    output_file_list = []
    alignment_dis = AlignIO.read(input_align_dis_path, 'fasta')
    for i, cluster in enumerate(filtered_cluster_records):
        output_file = os.path.join(
            output_dir, f'{os.path.basename(input_align_dis_path)}_g_{i + 1}.fa'
        )
        with open(output_file, 'w') as file:
            for seq_id in cluster:
                record = [
                    seq
                    for seq in alignment_dis
                    if seq.id.split('(')[0] == seq_id.split('_')[0]
                ]
                if record:
                    record = record[0]
                    file.write(f'>{record.id}\n')
                    file.write(f'{record.seq}\n')
        output_file_list.append(output_file)


def subset_bed_file(input_file, bed_dfs, output_dir):
    """
    Subset the given BED files by clusters.
    :param input_file: str, the absolute path of the BED file
    :return: a list containing the clustered BED files
    """

    output_file_list = []

    for i, df in enumerate(bed_dfs):
        output_file = os.path.join(
            output_dir, f'{os.path.basename(input_file)}_g_{i + 1}.bed'
        )
        df.to_csv(output_file, sep='\t', header=False, index=False)
        output_file_list.append(output_file)

    return output_file_list


def clean_and_cluster_MSA(
    input_file,
    bed_file,
    output_dir,
    div_column_thr=0.8,
    clean_column_threshold=0.01,
    min_length_num=10,
    cluster_num=2,
    cluster_col_thr=500,
    muscle_ite_times=4,
    fast_mode=False,
    input_msa=None,
):
    """
    This function will cluster multiple sequence alignment files.
    :param input_file: str, The direct FASTA file derived from the BED file
    :param bed_file: The BED file used to generate pattern alignments
    :param output_dir: Output directory
    :param gap_threshold: num (0-1), default 0.8, columns with gap percentage higher than "gap_threshold" will be removed
    :param clean_column_threshold: num (0-1), default 0.08, nucleotide percentage (gaps not count)
    lower than threshold will be converted to "-"
    :param min_length_num: num default 10, the minimum line number for each cluster
    :param cluster_num: num default 2, the maximum cluster number for each MSA
    :return: A list of subset pattern alignments and BED files
    """

    # When the input file is after multiple sequence alignment, don't do multiple sequence alignment again here
    if input_msa is None:

        # Align_sequences will return the absolute file path of alignment file
        if fast_mode:
            muscle_ite_times = 2
        try:
            fasta_out_flank_mafft_file = muscle_align(
                input_file, output_dir, ite_times=muscle_ite_times
            )
        except Exception:
            fasta_out_flank_mafft_file = False

        # When muscle goes wrong, use mafft
        if not fasta_out_flank_mafft_file:
            fasta_out_flank_mafft_file = align_sequences(input_file, output_dir)


        #fasta_out_flank_mafft_file = famsa_align(input_file, output_dir)

        # Remove gaps. Return absolute path for gap removed alignment file
        fasta_out_flank_mafft_file_gap_filter = remove_gaps_with_similarity_check(
            fasta_out_flank_mafft_file,
            output_dir,
            gap_threshold=0.8,
            simi_check_gap_thre=0.4,
            similarity_threshold=0.85,
            min_nucleotide=5,
            return_map=False,
        )
    else:
        fasta_out_flank_mafft_file_gap_filter = input_msa

    # Extract columns that contain different alignment patter, use this for group separation.
    # This threshold will be used to replace nucleotides that are less than threshold with a gap character
    # for each column
    # CleanAndSelectColumn removes sequences that have too many gaps along the sequence
    pattern_alignment = CleanAndSelectColumn(
        fasta_out_flank_mafft_file_gap_filter, threshold=clean_column_threshold
    )

    # Clean_column() function will help MSA cluster and return the absolute path of column cleaned alignment file
    pattern_alignment.clean_column(output_dir)

    # Select_divergent_column() function will return a boolean value, true represents need cluster step
    if pattern_alignment.select_divergent_column(
        cluster_col_thr=cluster_col_thr, dis_col_threshold=div_column_thr
    ):
        # write_alignment_filtered() function return pattern_alignment absolute path
        pattern_alignment = pattern_alignment.write_alignment_filtered(output_dir)
        """
        "min_lines" will define the minimum line numbers for each cluster.
        "max_cluster" will define the maximum cluster numbers that will be returned.
        This function will return a list containing the absolute paths of all cluster files.
        """

        filtered_cluster_records, if_cluster = cluster_msa_iqtree_DBSCAN(
            pattern_alignment, min_cluster_size=min_length_num, max_cluster=cluster_num
        )
        # Turn on the next line of code when need the separated MSA
        # subset_alignment_dis(filtered_cluster_records, output_dir, pattern_alignment)
    # if false, meaning no enough divergent column for cluster
    # Do full length alignment cluster and only keep the biggest cluster. The reason for this is the MSA is consistent
    # In the other word, the divergent columns is less. It is not necessary to do clustering.
    # But there could be some sequences in the MSA that is distinct to the MSA, clustering step will eliminate
    # those outline sequences.
    else:
        # Remove sequences that contain many gaps. This can cause problem for iqtree clustering.
        filter_out_big_gap_seq(
            fasta_out_flank_mafft_file_gap_filter,
            f'{fasta_out_flank_mafft_file_gap_filter}_gr.fa',
            gap_threshold=0.9,
        )

        # Set max_cluster to 1. This will eliminate outline sequences from the original MSA
        filtered_cluster_records, if_cluster = cluster_msa_iqtree_DBSCAN(
            f'{fasta_out_flank_mafft_file_gap_filter}_gr.fa',
            min_cluster_size=min_length_num,
            max_cluster=1,
        )
        #if if_cluster:
            # keep the biggest cluster only
        #    filtered_cluster_records = [max(filtered_cluster_records, key=len)]
    if if_cluster:
        # Subset bed file and return a list contain all clustered bed absolute files
        bed_dfs = process_labels(bed_file, filtered_cluster_records)
        cluster_bed_files_list = subset_bed_file(bed_file, bed_dfs, output_dir)
        return cluster_bed_files_list, fasta_out_flank_mafft_file_gap_filter
    else:
        """
        if cluster = False, it means no cluster has line numbers greater than "min_lines" or only has a small '-1' cluster. In this case,
        it will be difficult to use the multiple sequence alignment method to define the consensus sequence.
        This does not mean the sequence is not a TE, but has few copies in the genome. Low-copy TEs will be checked later.
        """
        return False


#####################################################################################################
# Code block: Plot MSA
#####################################################################################################


def alignment_to_dataframe(alignment):
    """Convert the alignment to a DataFrame"""
    # "rec" represents an object for each element in alignment. It contains "seq", "id", and "additional annotations"
    # list(rec) is equal to list(rec.seq) it will convert sequence to a list
    # "np.array" will convert list to array.
    alignment_df = pd.DataFrame(
        # Convert sequence nucleotides into array (nested list, pandas regards each element in the nested list
        # as a row
        np.array([list(rec) for rec in alignment], dtype=str),
        # Define the column number, be aware the difference between the len() and range()
        columns=list(range(len(alignment[0]))),
        # Define the row names of the dataframe
        index=[rec.id for rec in alignment],
    )
    return alignment_df


def create_base_mapping(alignment_df):
    """Create a dictionary with bases mapped to unique integers"""

    """
    The ravel function is a NumPy method that converts a multi-dimensional NumPy array into a flattened 1D array.
    The 'K' argument is an optional order parameter that specifies the memory layout of the result.
    The "alignment_df.values" part of the code retrieves the underlying NumPy array from the dataframe.
    The "pd.unique" function is used to find the unique elements of an array or a series. Here, it is used to identify
    the unique nucleotide bases (or gaps) that occur in the alignment.
    """
    unique_bases = pd.unique(alignment_df.values.ravel('K'))
    """
    The expression "base == base" is a way of testing whether a base has the value "NaN". This works because, according
    to IEEE floating point standards (which Python follows), "NaN" is not equal to anything, including itself.
        "NaN" is used to represent missing or undefined values.
    """
    unique_bases = [base for base in unique_bases if base == base]
    """
    "zip(unique_bases, range(len(unique_bases)))" is creating a pairing of each base with a unique integer in the range
    of the total number of unique bases.
    "dict()" is converting these pairings into a dictionary.
    """
    # Example of base_mapping: {'A':0, 'G':1, 'C':2, 'T':3, '-':4}
    # This function will find all available unique characters and label them with number
    base_mapping = dict(zip(unique_bases, range(len(unique_bases))))
    return base_mapping


def replace_bases(alignment_df, base_mapping):
    # Replace bases with corresponding integers from base_mapping
    result_obj = alignment_df.replace(base_mapping)
    # Infer the best possible data types for object columns
    result_obj = result_obj.infer_objects()
    return result_obj


def highlight_columns(alignment_df, alignment_color_df, base_mapping):
    for col in alignment_df.columns:
        column_without_gaps = [
            base for base in alignment_df[col] if base != '-'
        ]  # exclude gaps
        count = Counter(column_without_gaps)
        total_bases = sum(count.values())
        if total_bases == 0:  # all gaps, nothing to do
            continue
        # count.most_common(1)[0] is used to get the most common element (nucleotide)
        # in a column of the alignment and its count
        most_common_base, freq = count.most_common(1)[0]
        # Check if the most common nucleotide proportion is greater than 60%
        if freq / total_bases >= 0.6:
            # By default, all nucleotide are colored
            # If the most common nucleotide occupies more than 60%, don't color the whole column
            alignment_color_df[col] = base_mapping['-']

            # If the most common nucleotide occupies more than 60%, check every nucleotide, when the proportion
            # of them is less than 20%, color that kind of nucleotide again
            for base, freq in count.items():
                if freq / total_bases < 0.2:
                    alignment_color_df.loc[alignment_df[col] == base, col] = (
                        base_mapping[base]
                    )


def plot_msa(
    alignment_df,
    alignment_color_df,
    unique_bases,
    start_point,
    end_point,
    output_file,
    sequence_len,
):
    """
    Plot the heatmap
    """
    color_map = {
        'A': '#00CC00',
        'a': '#00CC00',
        'G': '#949494',
        'g': '#949494',
        'C': '#6161ff',
        'c': '#6161ff',
        'T': '#FF6666',
        't': '#FF6666',
        '-': '#FFFFFF',
        np.nan: '#FFFFFF',
    }
    # The nucleotide order in unique_bases match with the nucleotide converted number
    palette = [color_map[base] for base in unique_bases]
    # cmap connect the number with color type
    cmap = ListedColormap(palette)

    # Compute figure size: a bit less than number of columns / 5 for width, and number of rows / 3 for height
    # In pandas, the DataFrame.shape attribute returns a tuple representing the dimensionality of the DataFrame.
    # The tuple (r, c), where r represents the number of rows and c the number of columns.
    # figsize = (max(1, self.alignment_df.shape[1] / 6), max(9, self.alignment_df.shape[0] / 3))
    # Set a fixed height per sequence (e.g., 0.4 units per sequence)
    height_per_sequence = 0.4
    total_height = alignment_color_df.shape[0] * height_per_sequence

    figsize = (max(1, alignment_color_df.shape[1] / 6), total_height)
    # This function allows you to adjust several parameters that determine the size of the margins
    # left: This adjusts the margin on the left side of the plot. A value of 0.1, for instance,
    # means that the left margin will take up 10% of the total figure width.
    # right: This adjusts the margin on the right side of the plot. A value of 0.95 means that the right
    # margin starts at 95% of the figure width from the left. Essentially, it leaves a 5% margin on the right side.
    plt.figure(figsize=figsize, facecolor='white')
    plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.8)
    plt.imshow(alignment_color_df, aspect='auto', cmap=cmap)
    plt.box(False)  # Removing frame
    plt.xticks([])  # Remove labels
    plt.yticks([])  # Remove labels

    """
    "plt.annotate()" is a function from the Matplotlib library that is used to add text annotations to a figure.
    "xy=(start_point, -0.5)" specifies the point (x, y) that you want your arrow annotation to point to.
    "xytext=(start_point, -3)" specifies the position (x, y) to place the text.
    "arrowprops=dict(facecolor='red', edgecolor='red', shrink=0.05)" is a dictionary containing properties of the
    arrow that is drawn between the text and the point.
    "ha='center'" specifies the horizontal alignment of the text annotation.
    "color='r'" sets the color of the text to red.
    """
    plt.annotate(
        'Start crop Point',
        xy=(start_point, -0.5),
        xytext=(start_point, -3),
        arrowprops={'facecolor': 'red', 'edgecolor': 'red', 'shrink': 0.05},
        ha='center',
        color='r',
    )
    plt.annotate(
        'End crop Point',
        xy=(end_point, -0.5),
        xytext=(end_point, -3),
        arrowprops={'facecolor': 'blue', 'edgecolor': 'blue', 'shrink': 0.05},
        ha='center',
        color='b',
    )
    # Use input file name as title
    # title = os.path.basename(self.input_file)
    # plt.title(title, y=1.05, fontsize=15)

    # Add sequence len information to the plot
    plt.text(
        start_point + 49,
        -1,
        f'MSA length = {str(sequence_len)}',
        ha='left',
        va='center',
        color='black',
        size=15,
        weight='bold',
    )
    # Add the base letters
    # (j, i) is a tuple representing the indices of each element in the array, and
    # label is the value of the element at that position.
    for (j, i), label in np.ndenumerate(alignment_df):
        label = str(label).upper()
        if palette[alignment_color_df.iloc[j, i]] == '#FFFFFF':
            if label == '-':
                text_color = 'black'
            else:
                text_color = color_map.get(
                    label, 'black'
                )  # Default to "black" if label is not found
        else:
            text_color = 'black'

        plt.text(i, j, label, ha='center', va='center', color=text_color, size=13)

    plt.savefig(output_file, format='pdf', dpi=200)
    plt.close()
    # Release memory for these dataframes
    del alignment_df
    del alignment_color_df


def process_msa(input_file, output_dir, start_point, end_point, sequence_len):
    alignment = read_msa(input_file)
    alignment_df = alignment_to_dataframe(alignment)
    base_mapping = create_base_mapping(alignment_df)
    alignment_color_df = replace_bases(alignment_df, base_mapping)
    # DataFrames are mutable objects, which means highlight_columns() function doesn't need to return the modified
    # dataframe and the modification will be kept
    highlight_columns(alignment_df, alignment_color_df, base_mapping)
    output_file = os.path.join(output_dir, f'{os.path.basename(input_file)}_plot.pdf')
    plot_msa(
        alignment_df,
        alignment_color_df,
        base_mapping.keys(),
        start_point,
        end_point,
        output_file,
        sequence_len,
    )

    return output_file
