from Bio import AlignIO
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os


class MultipleSequenceAlignmentCluster:
    """
    A class for performing clustering on multiple sequence alignments (MSA)
    """

    def __init__(self, input_align_dis, output_dir, min_lines=10, max_cluster=2):
        """
        :param input_align_dis: str, the absolute path of MSA file that only contains divergent columns.
        :param output_dir: str, the absolute path of output folder.
        :param min_lines: num default 5, minimum line numbers for a cluster.
        :param max_cluster: num (2-5) default 2 , maximum cluster numbers.
        """

        self.input_align_dis = input_align_dis
        self.alignment_dis = AlignIO.read(self.input_align_dis, 'fasta')  # MSA file only contain distinct columns
        self.output_dir = output_dir
        self.min_lines = min_lines  # the minimum number of lines for each cluster
        self.max_cluster = max_cluster  # there could have more clusters, use this to define the maximum number
        self.encoder = None
        self.df = None
        self.df_scaled = None
        self.pca = None
        self.principal_components = None
        self.kmeans = None
        self.labels = None  # this variable contains cluster number for each sequence
        self.best_k = None  # this is the final cluster number
        self.if_cluster = True
        self.bed_dfs = []  # define dataframe bed file
        self.cluster_records = []  # define list to store different cluster ids
        self.filtered_cluster_records = []  # define list to store filtered cluster ids
        self.encode_sequences()
        self.standardize_data()
        self.perform_pca()
        self.determine_optimal_clusters()


    def encode_sequences(self):
        """
        Encode the sequences in the alignment using LabelEncoder
        This function will convert alignment nucleotide to numbers. a g c t represent four different number. Gpa is 0
        """
        encoder = LabelEncoder()
        data = []
        for record in self.alignment_dis:
            encoded_seq = encoder.fit_transform(list(str(record.seq)))
            data.append(encoded_seq)
        self.encoder = encoder
        self.df = pd.DataFrame(data)

    def standardize_data(self):
        """
        Standardize the encoded sequence data using StandardScaler
        This function will do Z score transformation z = (x - μ) / σ
        z: standardized value of a feature
        x: the original value of the feature
        μ: the means of the feature values
        σ: the standard deviation of the feature values
        This formula ensures that the standardized values have a mean of 0 and a standard deviation of 1.
        """
        scaler = StandardScaler()
        self.df_scaled = scaler.fit_transform(self.df)

    def perform_pca(self, n_components=2):
        """
        Perform Principal Component Analysis (PCA) on the standardized data
        :param n_components: num default 2, Number of components for PCA
        """
        self.pca = PCA(n_components=n_components)
        self.principal_components = self.pca.fit_transform(self.df_scaled)

    def determine_optimal_clusters(self):
        """
        Determine the optimal number of clusters using silhouette score
        """

        n_sequences = self.principal_components.shape[0]
        min_clusters = 2
        max_clusters = min(5, n_sequences - 1)  # Ensure max_clusters is less than the number of samples

        # Higher silhouette scores indicate better-defined clusters.
        # define a list to store the silhouette scores for different cluster numbers
        silhouette_scores = []
        list_k = list(range(min_clusters, max_clusters + 1))

        for k in list_k:

            # K-means is an iterative algorithm that aims to partition a dataset into k clusters by minimizing
            # the within-cluster sum of squares. The n_clusters parameter specifies the desired number of clusters
            # to be created.
            km = KMeans(n_clusters=k, n_init=15)
            km.fit(self.principal_components)
            score = silhouette_score(self.principal_components, km.labels_)
            silhouette_scores.append(score)  # get score for each cluster number

        # find the index of the maximum value in silhouette_scores
        best_index = np.argmax(silhouette_scores)
        print(f"{self.input_align_dis}\n{silhouette_scores}")

        # Test if cluster is necessary.
        if silhouette_scores[best_index] < 0.6:
            self.if_cluster = False

        # use the index to retrieve the corresponding cluster number from 'list_k'
        self.best_k = list_k[best_index]

    def apply_kmeans(self, k_n_init=15):
        """
        Apply K-means clustering with the determined optimal number of clusters
        """
        # Setting n_init to 20 means that the K-means algorithm will be run 20 times with different initializations,
        # and the final result will be the one with the lowest within-cluster sum of squares.
        self.kmeans = KMeans(n_clusters=self.best_k, n_init=k_n_init)
        self.kmeans.fit(self.principal_components)
        self.labels = self.kmeans.labels_

    def plot_pca(self):
        """
        Plot the PCA results with clusters
        """
        plt.figure(figsize=(8, 6))
        plt.scatter(
            self.principal_components[:, 0],
            self.principal_components[:, 1],
            c=self.labels
        )
        plt.title('PCA plot with clusters')
        plt_file = os.path.join(self.output_dir, f"{os.path.basename(self.input_align_dis)}_PCA.png")
        plt.savefig(plt_file)
        plt.close()

    def generate_cluster_list(self):
        """
        Generate a list contain different clusters
        """
        for i in range(self.best_k):
            # get line numbers for each clusters
            cluster_indices = [index for index, label in enumerate(self.labels) if label == i]
            # extract sequence id for each cluster
            cluster_records = [record.id for index, record in enumerate(self.alignment_dis) if index in cluster_indices]
            # Notice: the sequence id contain (+) or (-), when you want to use this to subset bed file
            self.cluster_records.append(cluster_records)

    def apply_filter_cluster(self):
        """
        Filter cluster
        """
        self.filtered_cluster_records = [cluster for cluster in self.cluster_records if len(cluster) >= self.min_lines]

        # Sort the filtered_groups by the number of lines in each group, in descending order
        self.filtered_cluster_records.sort(key=len, reverse=True)

        if len(self.filtered_cluster_records) > 2:
            # If more than 3 groups meet the condition, take the 2 groups with the most lines
            self.filtered_cluster_records = self.filtered_cluster_records[:self.max_cluster]

    def subset_alignment_intact(self, input_file):
        """
        According to clusters, subset input_file to different clusters
        :param input_file: str, the absolute path of input MSA file, which is used to generate pattern MSA file.
        :return: a list contains absolute paths of cluster MSA files
        """
        output_file_list = []

        alignment_intact = AlignIO.read(input_file, "fasta")
        for i, cluster in enumerate(self.filtered_cluster_records):
            output_file = os.path.join(self.output_dir, f"{os.path.basename(input_file)}_g_{i+1}.fa")
            with open(output_file, 'w') as file:
                for seq_id in cluster:
                    record = [seq for seq in alignment_intact if seq.id == seq_id]
                    if record:
                        record = record[0]
                        file.write(f'>{record.id}\n')
                        file.write(f'{record.seq}\n')
            output_file_list.append(output_file)

        return output_file_list

    def subset_alignment_dis(self):
        """
        This function will subset the input file means the pattern MSA.
        """
        output_file_list = []
        for i, cluster in enumerate(self.filtered_cluster_records):
            output_file = os.path.join(self.output_dir, f"{os.path.basename(self.input_align_dis)}_pat_MSA_g_{i+1}.fa")
            with open(output_file, 'w') as file:
                for seq_id in cluster:
                    record = [seq for seq in self.alignment_dis if seq.id == seq_id]
                    if record:
                        record = record[0]
                        file.write(f'>{record.id}\n')
                        file.write(f'{record.seq}\n')
            output_file_list.append(output_file)

        return output_file_list

    def subset_bed_file(self, input_file):
        """
        Subset the give bed files by clusters
        :param input_file: str, the absolute path of bed file
        :return: a list contain the clustered bed files
        """
        bed_df = pd.read_csv(input_file, sep='\t', header=None)
        for cluster in self.filtered_cluster_records:
            ids = [seq_id.rstrip("+-()") for seq_id in cluster] # remove (-) or (+) from ids
            cluster_df = bed_df[bed_df.apply(lambda row: f"{row[0]}:{row[1]}-{row[2]}" in ids, axis=1)]
            self.bed_dfs.append(cluster_df)

        # this function will return a list contain all bed file names
        output_file_list = []
        for i, df in enumerate(self.bed_dfs):
            output_file = os.path.join(self.output_dir, f"{os.path.basename(input_file)}_g_{i + 1}.bed")
            df.to_csv(output_file, sep='\t', header=False, index=False)
            output_file_list.append(output_file)
        return output_file_list
