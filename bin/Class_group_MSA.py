from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import os
import pandas as pd


class MultipleSequenceAlignmentCluster:

    def __init__(self, input_file, output_dir, min_cluster_size=10, max_cluster=2):

        self.input_file = input_file
        self.output_dir = output_dir
        self.min_cluster_size = min_cluster_size
        self.max_cluster = max_cluster
        self.if_cluster = True
        self.filtered_cluster_records = []
        self.MSA_cluster()

    def MSA_cluster(self):

        alignment = AlignIO.read(self.input_file, "fasta")

        # Calculate pairwise distances
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(alignment)

        # Convert the distance matrix to a numpy array
        distance_matrix = np.array(dm)

        # Define a range of cluster numbers to test
        cluster_range = range(2, 6)

        # Initialize variables to track the best cluster
        best_silhouette_score = -1
        best_cluster_assignments = None
        best_cluster_num = None

        # Iterate through different cluster numbers
        for n_clusters in cluster_range:
            # K-Means clustering
            kmeans = KMeans(n_clusters=n_clusters, random_state=0, n_init=10)
            cluster_assignments = kmeans.fit_predict(distance_matrix)

            # Calculate the Silhouette Score for the current clustering
            silhouette_score_value = silhouette_score(distance_matrix, cluster_assignments)

            # Update the best cluster if the current score is higher
            if silhouette_score_value > best_silhouette_score:
                best_silhouette_score = silhouette_score_value
                best_cluster_assignments = cluster_assignments
                best_cluster_num = n_clusters

        if best_silhouette_score < 0.6:

            self.if_cluster = False

        else:

            cluster_recording_list = []

            for i in range(best_cluster_num):

                # Get line numbers for each clusters
                cluster_indices = [index for index, label in enumerate(best_cluster_assignments) if label == i]

                # Extract sequence id for each cluster
                cluster_records = [record.id for index, record in enumerate(alignment) if index in cluster_indices]

                # Notice: the sequence id contain (+) or (-), when you want to use this to subset bed file
                cluster_recording_list.append(cluster_records)

            # Filter out cluster with less than min_cluster_size sequences
            self.filtered_cluster_records = [cluster for cluster in cluster_recording_list if
                                             len(cluster) >= self.min_cluster_size]

            # Sort the filtered_groups by the number of lines in each group, in descending order
            self.filtered_cluster_records.sort(key=len, reverse=True)

    def subset_bed_file(self, input_file):
        """
        Subset the give bed files by clusters
        :param input_file: str, the absolute path of bed file
        :return: a list contain the clustered bed files
        """
        bed_dfs = []
        bed_df = pd.read_csv(input_file, sep='\t', header=None)

        for cluster in self.filtered_cluster_records:
            ids = [seq_id.rstrip("+-()") for seq_id in cluster]  # Remove (-) or (+) from ids
            cluster_df = bed_df[bed_df.apply(lambda row: f"{row[0]}:{row[1]}-{row[2]}" in ids, axis=1)]
            bed_dfs.append(cluster_df)

        # Return a list contain all bed file names
        output_file_list = []

        for i, df in enumerate(bed_dfs):
            output_file = os.path.join(self.output_dir, f"{os.path.basename(input_file)}_g_{i + 1}.bed")
            df.to_csv(output_file, sep='\t', header=False, index=False)
            output_file_list.append(output_file)

        return output_file_list


