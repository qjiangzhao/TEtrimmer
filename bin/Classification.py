import subprocess
import os
from Bio import SeqIO
import click


# Classify single fasta
def classify_single(consensus_fasta):
    """
    Run RepeatClassifier with the provided parameters.
    """

    # Store the current working directory
    original_dir = os.getcwd()

    # Change the working directory to the directory of the consensus_fasta
    os.chdir(os.path.dirname(consensus_fasta))

    # Define RepeatClassifier command, the output file will store at the same directory of the consensus_fasta
    command = ["RepeatClassifier", "-consensi", consensus_fasta]

    # Run RepeatClassifier using subprocess
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Change the working directory back to the original directory
    os.chdir(original_dir)

    if result.returncode != 0:
        click.echo(f"RepeatClassifier error for {os.path.basename(consensus_fasta)}\n{result.stderr.decode('utf-8')}")
        return False

    classified_file = f'{consensus_fasta}.classified'

    # Get the first record of classified file
    record = next(SeqIO.parse(classified_file, "fasta"))
    seq_name = record.id.split("#")[0]
    seq_TE_type = record.id.split("#")[-1]

    return seq_TE_type
