import subprocess
import os
from Bio import SeqIO
import matplotlib.pyplot as plt

# classify uses the final consensus.fasta
def classify(consensus_file):
    """
    Run RepeatClassifier with the provided parameters.
    """
    # Construct the RepeatClassifier command
    command = ["RepeatClassifier",
               "-consensi", consensus_file
               ]
    # Run RepeatClassifier using subprocess
    subprocess.run(command, check=True)
    classified_file = f'{consensus_file}.classified'
    output_file, unknown_count, classified_count = summary(classified_file)
    plot_classification(output_file, unknown_count, classified_count)

def summary(classified_file):
    classification = {}
    unknown_count = 0
    with open(classified_file, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq_name = record.id.split("#")[0]
            seq_type = record.id.split("#")[-1]
            # is this necessary? It depends on if the input file checks redundancy
            if seq_name not in classification:
                classification[seq_name] = seq_type
                if "Unknown" in seq_type:
                    unknown_count += 1  
    total_count = len(classification)
    classified_count = total_count - unknown_count
    classified_pct = "{:.2%}".format(classified_count/total_count)
    # write to a summary file
    # TODO: if already has a summary file, need to modify
    output_file =f'{classified_file}.summary'
    with open(output_file, "w") as summary_output:
        summary_output.write(f"Total Count: {total_count}\n")
        summary_output.write(f"Classified Count: {classified_count}\n")
        summary_output.write(f"Unknown Count: {unknown_count}\n")
        summary_output.write(f"Classified Percentage: {classified_pct}\n")
        summary_output.write("\nClassification:\n")
        for seq_name, seq_type in classification.items():
            summary_output.write(f"{seq_name}: {seq_type}\n")
    
    return output_file, unknown_count, classified_count

# Visulization
def plot_classification(output_file, unknown_count, classified_count):
    # Data to plot
    output_file = f'{output_file}.png'
    labels = 'Unknown', 'Classified'
    sizes = [unknown_count, classified_count]
    colors = ['lightcoral', 'lightskyblue']
    explode = (0.1, 0)  # explode unknown slice

    # Create a pie chart
    # TODO: if doing comparision with the original classification, need to modify
    plt.pie(sizes, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', shadow=True, startangle=140)
    plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

    # Save the pie chart as a PNG file
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()

# Classify single fasta
def classify_single(dir, consensus_fasta):
    """
    Run RepeatClassifier with the provided parameters.
    """
    # Construct the RepeatClassifier command
    os.chdir(dir)
    stdout = open("stdout.txt", "w")
    stderr = open("stderr.txt", "w")
    command = ["RepeatClassifier",
            #    "-debug",
               "-consensi", consensus_fasta
               ]
    # Run RepeatClassifier using subprocess
    subprocess.run(command, check=True, stdout= stdout, stderr= stderr)
    classified_file = f'{consensus_fasta}.classified'
    for record in SeqIO.parse(classified_file, "fasta"):
        seq_name = record.id.split("#")[0]
        seq_TE_type = record.id.split("#")[-1]
    return seq_name, seq_TE_type