import csv
import os.path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import SeqIO


class OrfDomainPlot:

    def __init__(self, input_file, output_dir, orf_filepath, domain_filepath):
        """
        :param input_file: str, the input fasta file that only contain one sequence.
        :param output_dir: str, the absolute path of output directory.
        :param orf_filepath: str, the path of orf table.
        :param domain_filepath: str, the path of domain table.
        """
        self.input_file = input_file
        self.output_dir = output_dir
        self.orf_filepath = orf_filepath
        self.domain_filepath = domain_filepath
        self.fasta_length = None
        self.calculate_fasta_length()
        self.load_orfs()
        self.load_domains()


    # Calculate input fasta file length
    def calculate_fasta_length(self):
        with open(self.input_file, "r") as f:
            sequence = SeqIO.read(f, "fasta")  #SeqIO.read() is suitable for fasta file only contain one sequence.
            self.fasta_length = len(sequence.seq)

    def load_orfs(self):
        """
        Loads the ORF data from the file at orf_filepath.
        """
        orfs = {}
        with open(self.orf_filepath, 'r') as file:
            reader = csv.DictReader(file, delimiter='\t')
            for row in reader:
                key = (int(row['orf_start']), int(row['orf_end']), row['direction'], row['orf_name'])
                if key not in orfs:
                    orfs[key] = {'start': key[0], 'end': key[1], 'direction': key[2], 'name': key[3]}
        return orfs

    def load_domains(self):
        """
        Loads the Domain data from the file at domain_filepath.
        """
        domains = {}
        with open(self.domain_filepath, 'r') as file:
            reader = csv.DictReader(file, delimiter='\t')
            for row in reader:
                key = (int(row['domain_start']), int(row['domain_end']), row['direction'], row['domain_name'])
                if key not in domains:
                    domains[key] = {'start': key[0], 'end': key[1], 'direction': key[2], 'name': key[3]}
        return domains

    def plot_features(self, features, base_level, color):
        """
        Plot orf and domain features along sequence.
        """
        levels = []

        # Sort features according to the start position.
        for feature in sorted(features, key=lambda x: x['start']):
            for i, level in enumerate(levels):

                # When two features are close than 200bp, plot them at different levels.
                if feature['start'] - level['end'] > 200:
                    levels[i] = feature
                    break
            else:
                levels.append(feature)
                i = len(levels) - 1

            # when i is positive, feature will be plotted above the sequence, otherwise it will be beneath the sequence.
            if base_level < 0:
                i = -i

            start = feature['start']
            end = feature['end']
            length = end - start
            mid = start + length / 2

            if feature['direction'] == '+':
                plt.arrow(start, base_level + i * 0.2, length, 0, color=color, head_width=0.08,
                          head_length=0.003 * self.fasta_length, linewidth=1.5, shape='full')
                plt.text(mid, base_level + i * 0.2 + 0.05, feature['name'], fontsize=11, ha='center')
            elif feature['direction'] == '-':
                plt.arrow(end, base_level + i * 0.2, -length, 0, color=color, head_width=0.08,
                          head_length=0.003 * self.fasta_length, linewidth=1.5, shape='full')
                plt.text(mid, base_level + i * 0.2 + 0.05, feature['name'], fontsize=11, ha='center')

    def plot(self):
        """
        Loads the data and plots the ORFs and Domains.
        """
        output_file = os.path.join(self.output_dir, f"{os.path.basename(self.input_file)}_orf_pfam.pdf")
        orfs = self.load_orfs()
        domains = self.load_domains()

        fig, ax = plt.subplots(figsize=(20, 8))

        # Plot input sequence
        plt.plot([0, self.fasta_length], [0, 0], color='black', linewidth=5)

        # Plot orf features above the sequence
        self.plot_features(list(orfs.values()), 0.2, 'blue')

        # Plot domain features beneath the sequence
        self.plot_features(list(domains.values()), -0.2, 'red')

        # Create legend for plot
        blue_patch = mpatches.Patch(color='blue', label='ORFs')
        red_patch = mpatches.Patch(color='red', label='Domains')
        plt.legend(handles=[blue_patch, red_patch])

        plt.xlim(0, self.fasta_length)
        plt.ylim(-1, 2)
        plt.gca().axes.get_yaxis().set_visible(False)  # hides the y-axis
        plt.title('Orf and domain plot')
        plt.savefig(output_file, format="pdf")
        # Explicitly close the figure to free up memory
        plt.close(fig)

        return output_file

