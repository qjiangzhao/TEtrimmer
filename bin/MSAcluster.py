import os.path
from Bio import AlignIO
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from matplotlib.colors import ListedColormap
from Bio.Phylo.TreeConstruction import DistanceCalculator
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import os

from Class_select_ditinct_columns import CleanAndSelectColumn
import Function_blast_extension_mafft


def read_msa(input_file):
    """Read the MSA file"""
    return AlignIO.read(input_file, "fasta")


def clean_and_cluster_MSA(input_file, bed_file, output_dir, div_column_thr=0.8, clean_column_threshold=0.08,
                          min_length_num=10, cluster_num=2, cluster_col_thr=500, muscle_ite_times=4, fast_mode=False,
                          if_align=True):
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
    """
    if fast_mode:
        muscle_ite_times = 2
    try:
        fasta_out_flank_mafft_file = Function_blast_extension_mafft.muscle_align(input_file, output_dir,
                                                                                 ite_times=muscle_ite_times)
    except Exception as e:
        fasta_out_flank_mafft_file = False
        pass

    # When muscle goes wrong, use mafft
    if not fasta_out_flank_mafft_file:
        fasta_out_flank_mafft_file = Function_blast_extension_mafft.align_sequences(input_file, output_dir)
    """
    if if_align:
        fasta_out_flank_mafft_file = Function_blast_extension_mafft.align_sequences(input_file, output_dir)
    else:
        # In case the input file is aligned fasta file, it is not necessary to do it again
        fasta_out_flank_mafft_file = input_file
    # Remove gaps. Return absolute path for gap removed alignment file
    fasta_out_flank_mafft_file_gap_filter = Function_blast_extension_mafft.remove_gaps_with_similarity_check(
        fasta_out_flank_mafft_file, output_dir, gap_threshold=0.8, simi_check_gap_thre=0.4,
        similarity_threshold=0.85, min_nucleotide=5)

    # Extract columns that contain different alignment patter, use this for group separation.
    # This threshold will be used to replace nucleotides that are less than threshold with a gap character for each column
    pattern_alignment = CleanAndSelectColumn(fasta_out_flank_mafft_file_gap_filter, threshold=clean_column_threshold)

    # Clean_column() function will help MSA cluster and return the absolute path of column cleaned alignment file
    pattern_alignment.clean_column(output_dir)

    # write_alignment_filtered() function return pattern_alignment absolute path
    pattern_alignment = pattern_alignment.write_alignment_filtered(output_dir)

    # Select_divergent_column() function will return a boolean value, true represents need cluster step
    if pattern_alignment.select_divergent_column(cluster_col_thr=cluster_col_thr, dis_col_threshold=div_column_thr):

        """
        "min_lines" will define the minimum line numbers for each cluster
        "max_cluster" will define the maximum cluster numbers that will be returned
        this function will return a list contain all cluster file absolute path
        """

        # When the distinct occupied more than 35% of the MSA, lower silhouette score
        if len(pattern_alignment.dis_col_n) / pattern_alignment.alignment_seq_num >= 0.35:
            silhouette_score_thr = 0.5
        else:
            silhouette_score_thr = 0.6

        filtered_cluster_records, if_cluster = cluster_msa(pattern_alignment,
                                                           min_cluster_size=min_length_num,
                                                           max_cluster=cluster_num,
                                                           silhouette_score_thr=silhouette_score_thr)

        # Test if silhouette_scores is high enough to perform cluster
        if if_cluster:

            """
            Test if cluster number is 0, if so skip this sequence
            If cluster is 0, that means no cluster has line numbers greater than "min_lines". In this case, 
            it will be hard to still use multiple sequence alignment method to define consensus sequence.
            that isn't to say this won't be a TE, but with less copy numbers. Low copy TE will also be checked later
            """
            if len(filtered_cluster_records) == 0:
                return False
            else:
                # Subset bed file and return a list contain all clustered bed absolute files
                cluster_bed_files_list = subset_bed_file(bed_file, filtered_cluster_records, output_dir)
                return cluster_bed_files_list
        else:
            return True

    else:
        return True


def calculate_distance_matrix(alignment):
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)
    return np.array(dm)


def cluster_msa(alignment, min_cluster_size=10, max_cluster=False, silhouette_score_thr=0.6):
    alignment = read_msa(alignment)
    distance_matrix = calculate_distance_matrix(alignment)
    cluster_range = range(2, 6)
    # Initialize variables to track the best cluster
    best_silhouette_score = -1
    best_cluster_assignments = None
    best_cluster_num = None

    for n_clusters in cluster_range:
        kmeans = KMeans(n_clusters=n_clusters, random_state=0, n_init=10)
        cluster_assignments = kmeans.fit_predict(distance_matrix)
        silhouette_score_value = silhouette_score(distance_matrix, cluster_assignments)
        # Calculate the Silhouette Score for the current clustering
        # Update the best cluster if the current score is higher
        if silhouette_score_value > best_silhouette_score:
            best_silhouette_score = silhouette_score_value
            best_cluster_assignments = cluster_assignments
            best_cluster_num = n_clusters
    silhouette_score_thr = silhouette_score_thr

    if best_silhouette_score < silhouette_score_thr:
        return None, False

    cluster_recording_list = []
    for i in range(best_cluster_num):
        # Get line numbers for each clusters
        cluster_indices = [index for index, label in enumerate(best_cluster_assignments) if label == i]
        # Extract sequence id for each cluster
        cluster_records = [record.id for index, record in enumerate(alignment) if index in cluster_indices]
        cluster_recording_list.append(cluster_records)
    # Filter out cluster with less than min_cluster_size sequences
    filtered_cluster_records = [cluster for cluster in cluster_recording_list if len(cluster) >= min_cluster_size]
    # Sort the filtered_groups by the number of lines in each group, in descending order
    filtered_cluster_records.sort(key=len, reverse=True)

    if max_cluster:
        filtered_cluster_records = filtered_cluster_records[:max_cluster]

    return filtered_cluster_records, True


def subset_bed_file(input_file, filtered_cluster_records, output_dir):
    """
    Subset the give bed files by clusters
    :param input_file: str, the absolute path of bed file
    :return: a list contain the clustered bed files
    """
    bed_dfs = []
    bed_df = pd.read_csv(input_file, sep='\t', header=None)

    for cluster in filtered_cluster_records:
        ids = [seq_id.rstrip("+-()") for seq_id in cluster]  # Remove (-) or (+) from ids
        cluster_df = bed_df[bed_df.apply(lambda row: f"{row[0]}:{row[1]}-{row[2]}" in ids, axis=1)]
        bed_dfs.append(cluster_df)

    output_file_list = []

    for i, df in enumerate(bed_dfs):
        output_file = os.path.join(output_dir, f"{os.path.basename(input_file)}_g_{i + 1}.bed")
        df.to_csv(output_file, sep='\t', header=False, index=False)
        output_file_list.append(output_file)

    return output_file_list


def subset_alignment_intact(input_file, filtered_cluster_records, output_dir):
    """
    According to clusters, subset input_file to different clusters
    :param input_file: str, the absolute path of input MSA file, which is used to generate pattern MSA file.
    :return: a list contains absolute paths of cluster MSA files
    """
    output_file_list = []
    alignment_intact = AlignIO.read(input_file, "fasta")

    for i, cluster in enumerate(filtered_cluster_records):
        output_file = os.path.join(output_dir, f"{os.path.basename(input_file)}_g_{i + 1}.fa")
        with open(output_file, 'w') as file:
            for seq_id in cluster:
                record = [seq for seq in alignment_intact if seq.id == seq_id]
                if record:
                    record = record[0]
                    file.write(f'>{record.id}\n')
                    file.write(f'{record.seq}\n')
        output_file_list.append(output_file)

    return output_file_list

#####################################################################################################
# Code block: Plot MSA
#####################################################################################################


def alignment_to_dataframe(alignment):
    """Convert the alignment to a DataFrame"""
    # "rec" represents a object for each element in alignment. It contains "seq", "id", and "additional annotations"
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

def plot_msa(self, start_point, end_point):
    """
    Plot the heatmap
    """

    # Define your custom color map
    color_map = {"A": "#00CC00", "a": "#00CC00", "G": "#949494", "g": "#949494", "C": "#6161ff", "c": "#6161ff",
                 "T": "#FF6666", "t": "#FF6666", "-": "#FFFFFF", np.nan: "#FFFFFF"}

    # Create a custom color palette
    palette = [color_map[base] for base in self.unique_bases]
    cmap = ListedColormap(palette)

    # Compute figure size: a bit less than number of columns / 5 for width, and number of rows / 3 for height
    # In pandas, the DataFrame.shape attribute returns a tuple representing the dimensionality of the DataFrame.
    # The tuple (r, c), where r represents the number of rows and c the number of columns.
    #figsize = (max(1, self.alignment_df.shape[1] / 6), max(9, self.alignment_df.shape[0] / 3))

    # Set a fixed height per sequence (e.g., 0.4 units per sequence)
    height_per_sequence = 0.4
    total_height = self.alignment_df.shape[0] * height_per_sequence

    # Compute figure size
    figsize = (max(1, self.alignment_df.shape[1] / 6), total_height)

    # Plot the heatmap
    plt.figure(figsize=figsize, facecolor='white')

    # This function allows you to adjust several parameters that determine the size of the margins
    # left: This adjusts the margin on the left side of the plot. A value of 0.1, for instance,
    # means that the left margin will take up 10% of the total figure width.
    # right: This adjusts the margin on the right side of the plot. A value of 0.95 means that the right
    # margin starts at 95% of the figure width from the left. Essentially, it leaves a 5% margin on the right side.
    plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.8)
    plt.imshow(self.alignment_color_df, aspect='auto', cmap=cmap)
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
    plt.text(start_point + 49, -1, f"MSA length = {str(self.sequence_len)}",  ha='left',
                va='center', color='black', size=15, weight="bold")

    # Add the base letters
    # (j, i) is a tuple representing the indices of each element in the array, and
    # label is the value of the element at that position.
    for (j,i),label in np.ndenumerate(self.alignment_df):
        label = str(label).upper()  # Change the label to uppercase
        if palette[self.alignment_color_df.iloc[j, i]] == "#FFFFFF":
            if label == "-":
                text_color = "black"
            else:
                text_color = color_map.get(label, "black") # Default to "black" if label is not found
        else:
            text_color = "black"

        plt.text(i, j, label, ha='center', va='center', color=text_color, size=13)

    # Save the figure to a file
    self.output_file = os.path.join(self.output_dir, f"{os.path.basename(self.input_file)}_plot.pdf")
    plt.savefig(self.output_file, format='pdf', dpi=200)
    plt.close()

    # Release memory for these dataframes
    del self.alignment_df
    del self.alignment_color_df

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
    #figsize = (max(1, self.alignment_df.shape[1] / 6), max(9, self.alignment_df.shape[0] / 3))
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
    plt.text(start_point + 49, -1, f"MSA length = {str(sequence_len)}", ha='left', va='center', color='black', size=15, weight="bold")
    # Add the base letters
    # (j, i) is a tuple representing the indices of each element in the array, and
    # label is the value of the element at that position.
    for (j, i), label in np.ndenumerate(alignment_df):
        label = str(label).upper()
        if palette[alignment_color_df.iloc[j, i]] == "#FFFFFF":
            if label == "-":
                text_color = "black"
            else:
                text_color = color_map.get(label, "black") # Default to "black" if label is not found
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








