import subprocess
import os
from Bio import SeqIO
import click


# Classify single fasta
def classify_single(consensus_fasta):
    """
    Run RepeatClassifier with the provided parameters.
    """

    # Define RepeatClassifier command, the output file will store at the same directory of the consensus_fasta
    command = ["RepeatClassifier", "-consensi", consensus_fasta]

    # Run RepeatClassifier using subprocess
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if result.returncode != 0:
        click.echo(f"RepeatClassifier error for {os.path.basename(consensus_fasta)}")
        return False

    classified_file = f'{consensus_fasta}.classified'

    # Get the first record of classified file
    record = next(SeqIO.parse(classified_file, "fasta"))
    # seq_name = record.id.split("#")[0]
    seq_TE_type = record.id.split("#")[-1]

    return seq_TE_type
