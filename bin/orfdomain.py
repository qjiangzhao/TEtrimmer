import os
import os.path

import click
import csv
import subprocess
import urllib.request
import gzip
import shutil
import requests
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import SeqIO
from functions import prcyan, prgre


def check_and_download(directory, filename, url):
    """
    Function to check if file exist, otherwise download and unzip it
    :param directory: str, directory your checked file should be
    :param filename: str, file name
    :param url: str, url address will be used to download file when file isn't found in the given folder
    :return: boolean, True represent file don't found but successfully downloaded. False file is found
    """

    file_path = os.path.join(directory, filename)

    # If Pfam database can't be found download it
    if not os.path.isfile(file_path):

        click.echo(f"\n{filename} not found. Downloading... This might take some time. Please be patient\n")

        # Check if the URL is valid
        try:
            response = requests.get(url, stream=True)
            response.raise_for_status()  # Raise an exception if the GET request was unsuccessful
        except (requests.exceptions.HTTPError, requests.exceptions.ConnectionError) as e:
            prcyan(f"\nFailed to reach the server at {url} for Pfam database downloading\n")
            return False

        try:
            gz_file_path = file_path + ".gz"  # Give a defined name for the file will be downloaded
            urllib.request.urlretrieve(url, gz_file_path)

            # Unzipping
            with gzip.open(gz_file_path, 'rb') as f_in:
                with open(file_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

            # Delete gz file after extraction
            os.remove(gz_file_path)
            click.echo(f"\n{filename} is downloaded and unzipped.\n")
        except Exception:
            prcyan("TE Trimmer can't properly unzip the downloaded Pfam file")
            return False

    # Check if download is successful
    if os.path.isfile(file_path):
        return True
    else:
        prcyan(f"{filename} found. Pfam can't be downloaded by TE Trimmer. Please check your internet connection "
               f"or download Pfam by yourself")
        return False


def check_pfam_index_files(directory):
    """
    Check if the Pfam index files exist in the provided directory.
    :param directory: str, The directory to check for the Pfam index files.
    :return: boolean, True if all Pfam index files exist, False otherwise.
    """
    # Define the list of Pfam index files
    pfam_files = ["Pfam-A.hmm.h3f", "Pfam-A.hmm.h3m", "Pfam-A.hmm.h3i", "Pfam-A.hmm.h3p"]

    # Check if all Pfam index files exist
    if all(os.path.isfile(os.path.join(directory, file)) for file in pfam_files):
        return True
    else:
        # If one or more files are missing, delete all files
        for file in pfam_files:
            try:
                os.remove(os.path.join(directory, file))
            except FileNotFoundError:
                pass
        return False


def prepare_pfam_database(pfam_database_dir):
    pfam_hmm_url = r"https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
    pfam_dat_url = r"https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz"

    try:
        # Check and download the Pfam database files
        if check_and_download(pfam_database_dir, r"Pfam-A.hmm", pfam_hmm_url):
            if not check_pfam_index_files(pfam_database_dir):

                # Creates binary files that allow for faster access by other HMMER
                hmmpress_pfam_command = ["hmmpress", os.path.join(pfam_database_dir, "Pfam-A.hmm")]
                try:
                    subprocess.run(hmmpress_pfam_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

                except FileNotFoundError:
                    prcyan("'hmmpress' command not found. Please ensure 'hmmpress' is correctly installed.")
                    return False

                except subprocess.CalledProcessError as e:
                    prcyan(f"\nhmmpress index files generation failed with error code {e.returncode}")
                    prgre("Please check if your 'hmmpress' has been correctly installed.\n"
                          "Try 'hmmpress <your_downloaded_pfam_file> Pfam-A.hmm")
                    return False
        else:
            return False

    except Exception:
        return False

    try:
        if check_and_download(pfam_database_dir, r"Pfam-A.hmm.dat", pfam_dat_url):
            return True
        else:
            return False

    except Exception:
        return False


def determine_sequence_direction(filename):
    # Initialize the sums
    sum_positive = 0
    sum_negative = 0

    with open(filename, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')

        for row in reader:
            # Extract necessary values
            domain_start = int(row['domain_start'])
            domain_end = int(row['domain_end'])
            direction = row['direction']

            # Calculate domain length
            domain_length = domain_end - domain_start + 1

            # Add to the corresponding sum
            if direction == '+':
                sum_positive += domain_length
            elif direction == '-':
                sum_negative += domain_length

    # Compute the difference
    difference = sum_positive - sum_negative

    # Return True if difference is greater or equal than 0, else return False
    return difference >= 0


class PlotPfam:
    def __init__(self, input_file, output_dir, pfam_database_dir, mini_orf=200):
        """
        :param input_file: str, input fasta file path.
        :param output_dir: str, output file directory.
        :param pfam_database_dir: str, path where store pfam database. Will automatically download when don't find.
        :param mini_orf: num default 200, minimum orf length.
        """
        self.input_file = input_file
        self.input_file_n = os.path.basename(self.input_file)
        self.output_dir = output_dir
        self.pfam_database_dir = pfam_database_dir
        self.mini_orf = mini_orf
        self.output_orf_file_name_modified = None
        self.output_orf_file_name_modified_table = None
        self.output_pfam_file_modified = None

    def run_getorf(self):
        """
        Run getorf to extract orfs from given sequence.
        """
        output_orf_file = os.path.join(self.output_dir, f"{self.input_file_n}_orf.txt")
        self.output_orf_file_name_modified = os.path.join(self.output_dir,
                                                          f"{self.input_file_n}_orfm.txt")
        self.output_orf_file_name_modified_table = os.path.join(self.output_dir, f"{self.input_file_n}_orf_modi_t.txt")

        get_orf_command = [
            "getorf",
            "-sequence", self.input_file,
            "-outseq", output_orf_file,
            "-minsize", str(self.mini_orf)
        ]

        try:
            subprocess.run(get_orf_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        except FileNotFoundError:
            prcyan("getorf command not found. Please ensure that getorf is installed and available in your PATH.")
            raise Exception

        except subprocess.CalledProcessError as e:
            prcyan(f"\ngetorf failed for {self.input_file_n} with error code {e.returncode}.")
            prcyan('\n' + e.stderr + '\n')
            raise Exception

        # Check if the output_orf_file is empty
        if os.path.getsize(output_orf_file) == 0:
            return False

        change_orf_name = f"cat {output_orf_file} | awk '{{if(/>/){{print $1$2$3$4}}else{{print}}}}' > {self.output_orf_file_name_modified}"

        try:
            subprocess.run(change_orf_name, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            prcyan(f"\nFilter ORF columns failed for {self.input_file_n} with error code {e.returncode}")
            prcyan('\n' + e.stderr + '\n')
            raise Exception

        # Generate orf table for orf plot
        output_orf_file_name_modified_table = (
            f"cat {self.output_orf_file_name_modified} | grep '>' | "
            f"sed -e 's/>//g' -e 's/\\[/\\t/g' -e 's/-/\\t/g' -e 's/\\]//g' | "
            f"awk 'BEGIN{{OFS=\"\\t\"; print \"TE_name\", \"orf_start\", \"orf_end\", \"direction\", \"orf_name\"}}{{n=split($1, arr, \"_\"); "
            f"if($3>$2){{print $1, $2, $3, \"+\", \"ORF\"arr[n]}}else{{print $1, $3, $2, \"-\", \"ORF\"arr[n]}}}}' "
            f"> {self.output_orf_file_name_modified_table}"
        )
        try:
            subprocess.run(output_orf_file_name_modified_table, shell=True, check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            return True
        except subprocess.CalledProcessError as e:
            prcyan(f"\nConvert ORF result to table failed for {self.input_file_n} with "
                   f"error code {e.returncode}")
            prcyan('\n' + e.stderr + '\n')
            raise Exception

    def run_pfam_scan(self):
        """
        Run pfam_scan.pl to check orfs against Pfam database
        """
        # Define output file name and the modified pfam output file name.
        output_pfam_file = os.path.join(self.output_dir, f"{self.input_file_n}_orfm_pf.txt")
        self.output_pfam_file_modified = os.path.join(self.output_dir, f"{self.input_file_n}_orfm_pfm.txt")

        # Delete the output file if it already exists
        if os.path.exists(output_pfam_file):
            os.remove(output_pfam_file)

        pfam_sacn_command = [
            "pfam_scan.pl",
            "-fasta", self.output_orf_file_name_modified,
            "-dir", self.pfam_database_dir,
            "-outfile", output_pfam_file,
            "-cpu", str("1")
        ]
        try:
            subprocess.run(pfam_sacn_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        except FileNotFoundError:
            prcyan("'pfam_scan.pl' command not found. Please ensure 'pfam_scan.pl' is correctly installed.")
            raise Exception

        except subprocess.CalledProcessError as e:
            prcyan(f"\npfam_scan.pl failed for {self.input_file_n} with error code {e.returncode}")
            prcyan('\n' + e.stderr + '\n')
            raise Exception

        # Modify pfam output file for plot

        modify_pfam_result = (
            f"cat {output_pfam_file} | grep -v '#' | grep -v '^$' | "
            f"sed -e 's/\\[/ /g' -e 's/-/ /g' -e 's/\\]//g' | "
            f"awk 'BEGIN{{OFS=\"\\t\";print \"TE_name\", \"orf_start\", \"orf_end\", \"domain_start\", \"domain_end\", \"direction\", \"domain_name\", \"domain_reference\"}} "
            f"{{if($3>$2){{print $1, $2, $3, $2+$4*3-3, $2+$5*3-1, \"+\", $9, $8}}else{{print $1, $3, $2, $2-$5*3+1, $2-$4*3+3, \"-\", $9, $8}}}}' "
            f"> {self.output_pfam_file_modified}"
        )

        try:
            subprocess.run(modify_pfam_result, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            prcyan(f"\nTransform Pfam result failed for {self.input_file_n} with error code {e.returncode}")
            prcyan('\n' + e.stderr + '\n')
            raise Exception

        # Check if the output file only have one line, one line means the Pfam prediction is None.
        with open(self.output_pfam_file_modified, 'r') as file:
            line_count = sum(1 for _ in file)
            # Return True when line number is greater than 1, otherwise False
            return line_count > 1, self.output_pfam_file_modified

    def calculate_fasta_length(self,input_file):
        with open(input_file, "r") as f:
            sequence = SeqIO.read(f, "fasta")
            return len(sequence.seq)

    def load_orfs(self, orf_filepath):
        """
        Loads the ORF data from the file at orf_filepath.
        """
        orfs = {}
        with open(orf_filepath, 'r') as file:
            reader = csv.DictReader(file, delimiter='\t')
            for row in reader:
                key = (int(row['orf_start']), int(row['orf_end']), row['direction'], row['orf_name'])
                if key not in orfs:
                    orfs[key] = {'start': key[0], 'end': key[1], 'direction': key[2], 'name': key[3]}
        return orfs

    def load_domains(self,domain_filepath):
        """
        Loads the Domain data from the file at domain_filepath.
        """
        domains = {}
        with open(domain_filepath, 'r') as file:
            reader = csv.DictReader(file, delimiter='\t')
            for row in reader:
                key = (int(row['domain_start']), int(row['domain_end']), row['direction'], row['domain_name'])
                if key not in domains:
                    domains[key] = {'start': key[0], 'end': key[1], 'direction': key[2], 'name': key[3]}
        return domains

    def plot_features(self, features, base_level, color, fasta_length):
        levels = []

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
                        head_length=0.003 * fasta_length, linewidth=1.5, shape='full')
                plt.text(mid, base_level + i * 0.2 + 0.05, feature['name'], fontsize=11, ha='center')
            elif feature['direction'] == '-':
                plt.arrow(end, base_level + i * 0.2, -length, 0, color=color, head_width=0.08,
                        head_length=0.003 * fasta_length, linewidth=1.5, shape='full')
                plt.text(mid, base_level + i * 0.2 + 0.05, feature['name'], fontsize=11, ha='center')

    def plot(self,input_file, output_dir, orf_filepath, domain_filepath):
        fasta_length = self.calculate_fasta_length(input_file)
        orfs = self.load_orfs(orf_filepath)
        domains = self.load_domains(domain_filepath)

        fig, ax = plt.subplots(figsize=(20, 8))
        # Plot input sequence
        plt.plot([0, fasta_length], [0, 0], color='black', linewidth=5)
        # Plot orf features above the sequence
        self.plot_features(list(orfs.values()), 0.2, 'blue', fasta_length)
        # Plot domain features beneath the sequence
        self.plot_features(list(domains.values()), -0.2, 'red', fasta_length)
        # Create legend for plot
        blue_patch = mpatches.Patch(color='blue', label='ORFs')
        red_patch = mpatches.Patch(color='red', label='Domains')
        plt.legend(handles=[blue_patch, red_patch])

        plt.xlim(0, fasta_length)
        plt.ylim(-1, 2)
        plt.gca().axes.get_yaxis().set_visible(False)# hides the y-axis
        plt.title('Orf and domain plot')

        output_file = os.path.join(output_dir, f"{os.path.basename(input_file)}_orf_pfam.pdf")
        plt.savefig(output_file, format="pdf")
        # Explicitly close the figure to free up memory
        plt.close(fig)
        return output_file

    def orf_domain_plot(self):
        orf_domain_plot = self.plot(self.input_file, self.output_dir,
                                               self.output_orf_file_name_modified_table, self.output_pfam_file_modified)
        return orf_domain_plot
