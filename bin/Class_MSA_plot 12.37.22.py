import os.path
from Bio import AlignIO
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from matplotlib.colors import ListedColormap


class MSAPainter:

    def __init__(self, input_file, output_dir, sequence_len):

        self.input_file = input_file
        self.output_dir = output_dir
        self.output_file = None
        self.sequence_len = sequence_len
        self.alignment = None
        self.alignment_df = None
        self.alignment_color_df = None
        self.unique_bases = None
        self.base_mapping = None

    def read_msa(self):
        """Read the MSA file"""
        self.alignment = AlignIO.read(self.input_file, "fasta")

    def alignment_to_dataframe(self):
        """Convert the alignment to a DataFrame"""
        # "rec" represents a object for each element in alignment. It contains "seq", "id", and "additional annotations"
        # list(rec) is equal to list(rec.seq) it will convert sequence to a list
        # "np.array" will convert list to array.

        self.alignment_df = pd.DataFrame(
            np.array([list(rec) for rec in self.alignment], dtype=str),
            columns=[i for i in range(len(self.alignment[0]))],
            index=[rec.id for rec in self.alignment])

    def create_base_mapping(self):
        """Create a dictionary with bases mapped to unique integers"""

        """
        The ravel function is a numpy method that converts a multi-dimensional numpy array into a flattened 1D array.
        The 'K' argument is an optional order parameter that specifies the memory layout of the result.
        The alignment_df.values part of the code gets the underlying numpy array of the dataframe.
        The pd.unique function is used to find the unique elements of an array or a series. Here, it's used to identify the
        unique nucleotide bases (or gaps) that occur in your alignment.
        """
        self.unique_bases = pd.unique(self.alignment_df.values.ravel('K'))
        """
        The expression base == base is a clever way of testing whether base is "NaN". This works because, according to IEEE
        floating point standards (which Python follows), "NaN" is not equal to anything, including itself.
         "NaN" is used to represent missing or undefined values.
        """
        self.unique_bases = [base for base in self.unique_bases if base == base] # exclude NaN
        """
        zip(unique_bases, range(len(unique_bases))) is creating a pairing of each base with a unique integer in the range
        of the total number of unique bases.
        dict() is converting these pairings into a dictionary.
        """
        self.base_mapping = dict(zip(self.unique_bases, range(len(self.unique_bases))))

    def replace_bases(self):
        """Replace bases with corresponding integers"""
        self.alignment_color_df = self.alignment_df.replace(self.base_mapping)

        #self.alignment_color_df.to_csv("G:\\TE_manual_curation\\Software_develop\\Test_for_new_gap_and_crop_by_gap\\test_color.csv")

    def highlight_columns(self):
        """
        In each column, if no nucleotide occupies more than 80%, color the whole column.
        If there is a nucleotide occupying more than 80%, color the nucleotides that occupy less than 20%.
        """
        for col in self.alignment_df.columns:
            column_without_gaps = [base for base in self.alignment_df[col] if base != '-']  # exclude gaps
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
                self.alignment_color_df[col] = self.base_mapping['-']

                # If the most common nucleotide occupies more than 80%, check every nucleotides, when the proportion
                # of them is less than 20% color that kind of nucleotide again
                for base, freq in count.items():
                    if freq / total_bases < 0.2:
                        self.alignment_color_df[col][self.alignment_df[col] == base] = self.base_mapping[base]

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
        figsize = (max(1, self.alignment_df.shape[1] / 6), max(9, self.alignment_df.shape[0] / 3))

        # Plot the heatmap
        plt.figure(figsize=figsize)

        # This function allows you to adjust several parameters that determine the size of the margins
        # left: This adjusts the margin on the left side of the plot. A value of 0.1, for instance,
        # means that the left margin will take up 10% of the total figure width.
        # right: This adjusts the margin on the right side of the plot. A value of 0.95 means that the right
        # margin starts at 95% of the figure width from the left. Essentially, it leaves a 5% margin on the right side.
        plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.85)
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
                     arrowprops=dict(facecolor='green', edgecolor='green', shrink=0.05), ha='center', color='g')

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


    def process(self, start_point, end_point):
        self.read_msa()
        self.alignment_to_dataframe()
        self.create_base_mapping()
        self.replace_bases()
        self.highlight_columns()
        self.plot_msa(start_point, end_point)

        return self.output_file


