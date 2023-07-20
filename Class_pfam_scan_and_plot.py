import subprocess
import os
import urllib.request
import gzip
import shutil
import requests
from Class_orf_and_domain_plot import OrfDomainPlot

class PlotPfam:

    def __init__(self, input_file, output_dir, pfam_database_dir, mini_orf=200):
        """
        :param input_file: str, input fasta file path.
        :param output_dir: str, output file directory.
        :param pfam_database_dir: str, path where store pfam database. Will automatically download when don't find.
        :param mini_orf: num default 200, minimum orf length.
        """
        self.input_file = input_file
        self.output_dir = output_dir
        self.pfam_database_dir = pfam_database_dir
        self.mini_orf = mini_orf
        self.output_orf_file_name_modified = None
        self.output_orf_file_name_modified_table = None
        self.output_pfam_file_modified = None
        self.run_getorf()
        self.run_pfam_scan()


    def check_and_download(self, directory, filename, url):
        """
        Function to check if file exist, otherwise download and unzip it
        :param directory: str, directory your checked file should be
        :param filename: str, file name
        :param url: str, url address will be used to download file when file isn't found in the given folder
        :return: boolean, True represent file don't found but successfully downloaded. False file is found
        """
        # Check if the URL is valid
        try:
            response = requests.get(url, stream=True)
            response.raise_for_status()  # Raise an exception if the GET request was unsuccessful
        except (requests.exceptions.HTTPError, requests.exceptions.ConnectionError) as e:
            print(f"Failed to reach the server at {url}. Please check the URL.")
            print("Error:", str(e))
            return

        # Check if file exists
        file_path = os.path.join(directory, filename)
        if not os.path.isfile(file_path):
            print(f"{filename} not found. Downloading... This might take some time. Please be patient")
            gz_file_path = file_path + ".gz"  # Give a defined name for the file will be downloaded
            urllib.request.urlretrieve(url, gz_file_path)

            # Unzipping
            with gzip.open(gz_file_path, 'rb') as f_in:
                with open(file_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

            # Delete gz file after extraction
            os.remove(gz_file_path)
            print(f"{filename} downloaded and unzipped.")

            return True
        else:
            print(f"{filename} found.")
            return False

    #####################################################################################################
    # Code block: Get orf by "getorf" convert the output fasta sequence name and generate orf table
    #####################################################################################################

    def run_getorf(self):
        """
        Run getorf to extract orfs from given sequence.
        """
        output_orf_file = os.path.join(self.output_dir, f"{os.path.basename(self.input_file)}_orf.txt")
        self.output_orf_file_name_modified = os.path.join(self.output_dir, f"{os.path.basename(self.input_file)}_orf_modi.txt")
        self.output_orf_file_name_modified_table= os.path.join(self.output_dir, f"{os.path.basename(self.input_file)}_orf_modi_t.txt")

        get_orf_command = [
            "getorf",
            "-sequence", self.input_file,
            "-outseq", output_orf_file,
            "-minsize", str(self.mini_orf)
        ]

        try:
            subprocess.run(get_orf_command, check=True)
        except FileNotFoundError:
            print("getorf command not found. Please ensure that getorf is installed and available in your PATH.")
            return
        except subprocess.CalledProcessError:
            print("getorf command failed. Please check your input parameters.")
            return

        change_orf_name = f"cat {output_orf_file} | awk '{{if(/>/){{print $1$2$3$4}}else{{print}}}}' > {self.output_orf_file_name_modified}"

        subprocess.run(change_orf_name, shell=True, check=True)

        # Generate orf table for orf plot

        output_orf_file_name_modified_table = (
            f"cat {self.output_orf_file_name_modified} | grep '>' | "
            f"sed -e 's/>//g' -e 's/\\[/\\t/g' -e 's/-/\\t/g' -e 's/\\]//g' | "
            f"awk 'BEGIN{{OFS=\"\\t\"; print \"TE_name\", \"orf_start\", \"orf_end\", \"direction\", \"orf_name\"}}{{n=split($1, arr, \"_\"); "
            f"if($3>$2){{print $1, $2, $3, \"+\", \"ORF\"arr[n]}}else{{print $1, $3, $2, \"-\", \"ORF\"arr[n]}}}}' "
            f"> {self.output_orf_file_name_modified_table}"
        )

        subprocess.run(output_orf_file_name_modified_table, shell=True, check=True)

    #####################################################################################################
    # Code block: Check pfam database and run pfam_scan.pl
    #####################################################################################################

    def run_pfam_scan(self):
        """
        Run pfam_scan.pl to check orfs against Pfam database
        """
        # Check if pfam database is available from the give directory, otherwise download it.
        pfam_hmm_url = r"https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
        pfam_dat_url = r"https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz"

        if self.check_and_download(self.pfam_database_dir, r"Pfam-A.hmm", pfam_hmm_url):
            hmmpress_pfam_command = [
                "hmmpress" , os.path.join(self.pfam_database_dir, "Pfam-A.hmm")
            ]

            # Check if hmmpress tool is installed.
            try:
                subprocess.run(hmmpress_pfam_command, check=True)
            except FileNotFoundError:
                print(
                    "hmmpress command not found. Please ensure that hmmpress is installed and available in your PATH.")
                return

        # Check if Pfam-A.hmm.dat exists, which is required for tool pfam_scan.pl. Otherwise download it.
        self.check_and_download(self.pfam_database_dir, r"Pfam-A.hmm.dat", pfam_dat_url)

        # Define output file name and the modified pfam output file name.
        output_pfam_file = os.path.join(self.output_dir, f"{os.path.basename(self.input_file)}_orf_modi_pfam.txt")
        self.output_pfam_file_modified = os.path.join(self.output_dir, f"{os.path.basename(self.input_file)}_orf_modi_pfam_modi.txt")
        pfam_sacn_command = [
            "pfam_scan.pl",
            "-fasta", self.output_orf_file_name_modified,
            "-dir", self.pfam_database_dir,
            "-outfile", output_pfam_file
        ]
        try:
            subprocess.run(pfam_sacn_command, check=True)
        except FileNotFoundError:
            print("pfam_scal.pl command not found. Please ensure that tool is installed and available in your PATH.")

        # Modify pfam output file for plot

        modify_pfam_result = (
            f"cat {output_pfam_file} | grep -v '#' | grep -v '^$' | "
            f"sed -e 's/\\[/ /g' -e 's/-/ /g' -e 's/\\]//g' | "
            f"awk 'BEGIN{{OFS=\"\\t\";print \"TE_name\", \"orf_start\", \"orf_end\", \"domain_start\", \"domain_end\", \"direction\", \"domain_name\", \"domain_reference\"}} "
            f"{{if($3>$2){{print $1, $2, $3, $2+$4*3-3, $2+$5*3-1, \"+\", $9, $8}}else{{print $1, $3, $2, $2-$5*3+1, $2-$4*3+3, \"-\", $9, $8}}}}' "
            f"> {self.output_pfam_file_modified}"
        )

        subprocess.run(modify_pfam_result, shell=True, check=True)

    #####################################################################################################
    # Code block: Plot orf and domain features
    #####################################################################################################

    def orf_domain_plot(self):
        orf_domain_plot_object = OrfDomainPlot(self.input_file, self.output_dir,
                                               self.output_orf_file_name_modified_table,self.output_pfam_file_modified )
        orf_domain_plot = orf_domain_plot_object.plot()

        return orf_domain_plot







