import os.path
from Bio import AlignIO
from Bio import Phylo
import pandas as pd
import numpy as np
import subprocess
import click
import matplotlib.pyplot as plt
from collections import Counter
from matplotlib.colors import ListedColormap
import os
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
import re

from selectcolumns import CleanAndSelectColumn
from functions import prcyan, prgre, muscle_align, align_sequences, remove_gaps_with_similarity_check, \
    filter_out_big_gap_seq


def calculate_tree_dis(input_file):
    # Read the tree from a Newick file
    tree = Phylo.read(input_file, "newick")
    output_file = f"{input_file}.tdistance"
    # Get all terminals (leaves) in the tree
    terminals = tree.get_terminals()
    num_terminals = len(terminals)

    # Initialize an empty distance matrix with sequence names
    distance_matrix = [[""] * (num_terminals + 1) for _ in range(num_terminals + 1)]

    # Set the first row and column as sequence names
    for i in range(1, num_terminals + 1):
        distance_matrix[i][0] = terminals[i - 1].name
        distance_matrix[0][i] = terminals[i - 1].name

    # Calculate pairwise distances
    for i in range(1, num_terminals + 1):
        for j in range(1, num_terminals + 1):
            # If i and j are equal, set distance to 1
            if i == j:
                distance_matrix[i][j] = 1
            else:
                # Calculate the total branch length between terminal nodes i and j
                total_branch_length = tree.distance(terminals[i - 1], terminals[j - 1])

                # Assign the calculated distance to the matrix
                distance_matrix[i][j] = total_branch_length

    # Write the distance matrix to a file
    with open(output_file, "w") as f:
        for row in distance_matrix:
            f.write("\t".join(map(str, row)) + "\n")
    return output_file


def read_msa(input_file):
    """Read the MSA file"""
    return AlignIO.read(input_file, "fasta")


def cluster_msa_iqtree_DBSCAN(alignment, min_cluster_size=10, max_cluster=None):
    # this function uses iqtree to separate alignment into branches and then do dbscan_cluster on sequences
    # based on maximum likelihood tree distance
    # dependencies: calculate_tree_dis, dbscan_cluster
    # temp using K2P+I
    iqtree_command = ["iqtree",
                      "-m",
                      "K2P+I",
                      "--redo-tree",
                      "--seqtype",
                      "DNA",
                      "-s", 
                      alignment]
    try:
        subprocess.run(iqtree_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    except FileNotFoundError:
        prcyan("'iqtree' command not found. Please ensure 'iqtree' is correctly installed.")
        raise Exception

    except subprocess.CalledProcessError as e:
        prcyan(f"iqtree failed with error code {e.returncode}")
        prcyan(e.stderr)
        raise Exception("iqtree error")

    treefile = f"{alignment}.treefile"
    distance_file = calculate_tree_dis(treefile)
    cluster, sequence_names = dbscan_cluster(distance_file, pca=False)

    # map sequence_names name to cluster 
    sequence_cluster_mapping = dict(zip(sequence_names, cluster))

    # Use Counter to count occurrences
    counter = Counter(cluster)

    # Find cluster with size > min_cluster_size
    filter_cluster = [element for element, count in counter.items() if count > min_cluster_size]
    seq_cluster_list = []

    # if only -1 in filter_cluster and cluster size > 60% number of sequences in MSA:
    if filter_cluster == [-1]:
        if counter[-1] > 0.6 * len(sequence_names):
            #print (f"only -1 cluster > 10, size = {counter[-1]}, thr = {0.6 * len(sequence_names)}")
            seq_records = [sequence for sequence, cluster_label in sequence_cluster_mapping.items() if cluster_label == -1]
            seq_cluster_list.append(seq_records)
    else: 
        # has 0 cluster, or 1 cluster that's not -1, or multiple clusters, elimininate cluster that is -1
        filter_cluster = [c for c in filter_cluster if c != -1]
        top_cluster = Counter({element: count for element, count in counter.items() if element in filter_cluster}).most_common
        # Pick the top max_cluster
        if max_cluster is not None:
            top_cluster = top_cluster(max_cluster)
        else:
            top_cluster = top_cluster()
        # Find the seq name in selected cluster
        if len(top_cluster) > 0:
            for i in top_cluster:
                # Extract sequence id for each cluster
                seq_records = [sequence for sequence, cluster_label in sequence_cluster_mapping.items() if cluster_label == i[0]]
                seq_cluster_list.append(seq_records)
    
    if_cluster = (len(seq_cluster_list) > 0)
    # if all cluster size < 10, there is no meaningful cluster, if_cluster = False, skip this sequence
    # if only one cluster has size >10 and it is -1,
    #       if -1 cluster > 60% number of sequences in MSA, use this -1 cluster
    #       else, skip this sequence 
    # if at least one cluster has size >10 and it is not just -1, 
    #   only keep the max_cluster that are not -1
    return seq_cluster_list, if_cluster


def dbscan_cluster(input_file, pca=True):
    # Read the distance matrix from the file
    with open(input_file, "r") as f:
        lines = f.readlines()

    # Extract distance values
    distance_values = []
    sequence_names = []  # Add this line to store sequence names
    for line in lines[1:]:
        values = line.strip().split("\t")
        sequence_names.append(values[0])  # Assuming the sequence name is in the first column
        float_values = [float(value) for value in values[1:]]  # Adjust the index based on your data
        distance_values.append(float_values)
        
    # Convert distance values to a NumPy array
    distance_np = np.array(distance_values, dtype=np.float64)

    # Replace NaN values with a large number
    distance_np[np.isnan(distance_np)] = np.nanmax(distance_np) + 1

    # Apply DBSCAN clustering with the recommended eps value
    dbscan = DBSCAN(eps=0.1, min_samples=2, metric="precomputed")
    labels = dbscan.fit_predict(distance_np)
    if pca:
        title = os.path.basename(input_file)
        plot_pca(distance_np, labels, input_file, title=title)
        
    return labels, sequence_names


def plot_pca(distance_matrix, labels, output_file, title="PCA Plot"):
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
        plt.scatter(cluster_points['PC1'], cluster_points['PC2'], label=f'Cluster {label}')

        # Calculate the count of points in the cluster
        cluster_counts[label] = len(indices)

        # Annotate the cluster with its point count
        # The annotation is placed at the mean position of the cluster's points
        plt.annotate(f'Count: {len(indices)}',
                     (cluster_points['PC1'].mean(), cluster_points['PC2'].mean()),
                     textcoords="offset points", xytext=(0,10), ha='center')

    plt.title(f'{title} - PCA Plot')
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.legend()

    # Save the plot to a file
    output_file_n = f"{output_file}_PCA.pdf"
    plt.savefig(output_file_n)
    plt.close()


def process_labels(input_file, filtered_cluster_records):
    """
    during the iqtree process, the sequence name is modified to remove special characters
    this step extract bed info based on the index of sequence name in bedfile created by nameOnly
    :param input_file: str, the absolute path of bed file
    :return: a list contain the clustered bed files
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


def subset_bed_file(input_file, bed_dfs, output_dir):
    """
    Subset the give bed files by clusters
    :param input_file: str, the absolute path of bed file
    :return: a list contain the clustered bed files
    """

    output_file_list = []

    for i, df in enumerate(bed_dfs):
        output_file = os.path.join(output_dir, f"{os.path.basename(input_file)}_g_{i + 1}.bed")
        df.to_csv(output_file, sep='\t', header=False, index=False)
        output_file_list.append(output_file)

    return output_file_list


def clean_and_cluster_MSA(input_file, bed_file, output_dir, div_column_thr=0.8, clean_column_threshold=0.08,
                          min_length_num=10, cluster_num=2, cluster_col_thr=500, muscle_ite_times=4, fast_mode=False):
    """
    This function will cluster multiple sequence alignment file
    :param input_file: str, The direct fasta file derived from bed file
    :param bed_file: The bed file used to generate pattern alignment
    :param output_dir: Output directory
    :param gap_threshold: num (0-1) default 0.8, columns with gap percentage higher than "gap_threshold" will be removed
    :param clean_column_threshold: num (0-1) default 0.08, nucleotide percentage (gap not count)
    lower than threshold will be converted to "-"
    :param min_length_num: num default 10, the minimum line number for each cluster
    :param cluster_num: num default 2, the maximum cluster number for each MSA
    :return: A list of subset pattern alignment and bed files
    """

    # Align_sequences will return the absolute file path of alignment file
    if fast_mode:
        muscle_ite_times = 2
    try:
        fasta_out_flank_mafft_file = muscle_align(input_file, output_dir, ite_times=muscle_ite_times)
    except Exception as e:
        fasta_out_flank_mafft_file = False
        pass

    # When muscle goes wrong, use mafft
    if not fasta_out_flank_mafft_file:
        fasta_out_flank_mafft_file = align_sequences(input_file, output_dir)

    # Remove gaps. Return absolute path for gap removed alignment file
    fasta_out_flank_mafft_file_gap_filter = remove_gaps_with_similarity_check(
        fasta_out_flank_mafft_file, output_dir, gap_threshold=0.8, simi_check_gap_thre=0.4,
        similarity_threshold=0.85, min_nucleotide=5)

    # Extract columns that contain different alignment patter, use this for group separation.
    # This threshold will be used to replace nucleotides that are less than threshold with a gap character for each column
    pattern_alignment = CleanAndSelectColumn(fasta_out_flank_mafft_file_gap_filter, threshold=clean_column_threshold)

    # Clean_column() function will help MSA cluster and return the absolute path of column cleaned alignment file
    pattern_alignment.clean_column(output_dir)

    # Select_divergent_column() function will return a boolean value, true represents need cluster step
    if pattern_alignment.select_divergent_column(cluster_col_thr=cluster_col_thr, dis_col_threshold=div_column_thr):

        # write_alignment_filtered() function return pattern_alignment absolute path
        pattern_alignment = pattern_alignment.write_alignment_filtered(output_dir)
        """
        "min_lines" will define the minimum line numbers for each cluster
        "max_cluster" will define the maximum cluster numbers that will be returned
        this function will return a list contain all cluster file absolute path
        """

        filtered_cluster_records, if_cluster = cluster_msa_iqtree_DBSCAN(pattern_alignment,
                                                                         min_cluster_size=min_length_num,
                                                                         max_cluster=cluster_num)
    # if false, meaning no enough divergent column for cluster
    # Do full length alignment cluster and only keep the biggest cluster
    else:

        # Remove sequences that contain many gaps. This can cause problem for iqtree clustering.
        filter_out_big_gap_seq(fasta_out_flank_mafft_file_gap_filter, fasta_out_flank_mafft_file_gap_filter,
                               gap_threshold=0.9)
        filtered_cluster_records, if_cluster = cluster_msa_iqtree_DBSCAN(
            fasta_out_flank_mafft_file_gap_filter, min_cluster_size=min_length_num, max_cluster=cluster_num)
        if if_cluster:
            # keep the biggest cluster only
            filtered_cluster_records = [max(filtered_cluster_records, key=len)]
    if if_cluster:
        # Subset bed file and return a list contain all clustered bed absolute files
        bed_dfs = process_labels(bed_file, filtered_cluster_records)
        cluster_bed_files_list = subset_bed_file(bed_file, bed_dfs, output_dir)
        return cluster_bed_files_list
    else:
        """
        if cluster = False, it means no cluster has line numbers greater than "min_lines" or only has small -1 cluster. In this case, 
        it will be hard to still use multiple sequence alignment method to define consensus sequence.
        that isn't to say this won't be a TE, but with less copy numbers. Low copy TE will also be checked later
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
        np.array([list(rec) for rec in alignment], dtype=str),
        columns=[i for i in range(len(alignment[0]))],
        index=[rec.id for rec in alignment])
    return alignment_df


def create_base_mapping(alignment_df):
    """Create a dictionary with bases mapped to unique integers"""

    """
    The ravel function is a numpy method that converts a multi-dimensional numpy array into a flattened 1D array.
    The 'K' argument is an optional order parameter that specifies the memory layout of the result.
    The alignment_df.values part of the code gets the underlying numpy array of the dataframe.
    The pd.unique function is used to find the unique elements of an array or a series. Here, it's used to identify the
    unique nucleotide bases (or gaps) that occur in your alignment.
    """
    unique_bases = pd.unique(alignment_df.values.ravel('K'))
    """
    The expression base == base is a clever way of testing whether base is "NaN". This works because, according to IEEE
    floating point standards (which Python follows), "NaN" is not equal to anything, including itself.
        "NaN" is used to represent missing or undefined values.
    """
    unique_bases = [base for base in unique_bases if base == base]
    """
    zip(unique_bases, range(len(unique_bases))) is creating a pairing of each base with a unique integer in the range
    of the total number of unique bases.
    dict() is converting these pairings into a dictionary.
    """
    base_mapping = dict(zip(unique_bases, range(len(unique_bases))))
    return base_mapping


def replace_bases(alignment_df, base_mapping):
    """Replace bases with corresponding integers"""
    return alignment_df.replace(base_mapping)


def highlight_columns(alignment_df, alignment_color_df, base_mapping):
    for col in alignment_df.columns:
        column_without_gaps = [base for base in alignment_df[col] if base != '-']  # exclude gaps
        count = Counter(column_without_gaps)
        total_bases = sum(count.values())
        if total_bases == 0:  # all gaps, nothing to do
            continue
        # count.most_common(1)[0] is used to get the most common element (nucleotide)
        # in a column of the alignment and its count
        most_common_base, freq = count.most_common(1)[0]
        # Check if the most common nucleotide proportion is greater than 80%
        if freq / total_bases >= 0.6:
            # If the most common nucleotide occupies more than 80%, don't color the whole column
            alignment_color_df[col] = base_mapping['-']

            # If the most common nucleotide occupies more than 80%, check every nucleotides, when the proportion
            # of them is less than 20% color that kind of nucleotide again
            for base, freq in count.items():
                if freq / total_bases < 0.2:
                    alignment_color_df[col][alignment_df[col] == base] = base_mapping[base]


def plot_msa(alignment_df, alignment_color_df, unique_bases, start_point, end_point, output_file, sequence_len):
    """
    Plot the heatmap
    """
    color_map = {"A": "#00CC00", "a": "#00CC00", "G": "#949494", "g": "#949494", "C": "#6161ff", "c": "#6161ff",
                 "T": "#FF6666", "t": "#FF6666", "-": "#FFFFFF", np.nan: "#FFFFFF"}

    palette = [color_map[base] for base in unique_bases]
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
    plt.annotate() is a function from the matplotlib library that is used to add text annotations to a figure.
    xy=(start_point, -0.5) specifies the point (x, y) that you want your arrow annotation to point to. 
    xytext=(start_point, -3) specifies the position (x, y) to place the text.
    arrowprops=dict(facecolor='red', edgecolor='red', shrink=0.05) is a dictionary containing properties of the arrow 
    that is drawn from the text to the point. 
    ha='center' specifies the horizontal alignment of the text annotation. 
    color='r' sets the color of the text to red.
    """
    plt.annotate('Start crop Point', xy=(start_point, -0.5), xytext=(start_point, -3),
                 arrowprops=dict(facecolor='red', edgecolor='red', shrink=0.05), ha='center', color='r')
    plt.annotate('End crop Point', xy=(end_point, -0.5), xytext=(end_point, -3),
                 arrowprops=dict(facecolor='blue', edgecolor='blue', shrink=0.05), ha='center', color='b')
    # Use input file name as title
    # title = os.path.basename(self.input_file)
    # plt.title(title, y=1.05, fontsize=15)

    # Add sequence len information to the plot
    plt.text(start_point + 49, -1, f"MSA length = {str(sequence_len)}", ha='left', va='center', color='black', size=15,
             weight="bold")
    # Add the base letters
    # (j, i) is a tuple representing the indices of each element in the array, and
    # label is the value of the element at that position.
    for (j, i), label in np.ndenumerate(alignment_df):
        label = str(label).upper()
        if palette[alignment_color_df.iloc[j, i]] == "#FFFFFF":
            if label == "-":
                text_color = "black"
            else:
                text_color = color_map.get(label, "black")  # Default to "black" if label is not found
        else:
            text_color = "black"

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
    highlight_columns(alignment_df, alignment_color_df, base_mapping)
    output_file = os.path.join(output_dir, f"{os.path.basename(input_file)}_plot.pdf")
    plot_msa(alignment_df, alignment_color_df, base_mapping.keys(), start_point, end_point, output_file, sequence_len)

    return output_file
