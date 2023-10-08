from Bio import SeqIO
import os
from Class_seq_object import SeqObject


def separate_sequences(input_file, output_dir, continue_analysis=False):
    """
    separates input file into single fasta file and creates object for each input sequence
    """
    os.makedirs(output_dir, exist_ok=True)
    seq_list = []

    if not continue_analysis:

        print(
            "TE Trimmer is modifying sequence names, '/', '-', ':' and empty space before '#' will be converted to '_'")
        detected_pound = False
        with open(input_file, 'r') as fasta_file:
            for record in SeqIO.parse(fasta_file, 'fasta'):
                if len(record.id.split("#")) > 1:
                    detected_pound = True
                    sanitized_id = record.id.split("#")[0].replace('/', '_').replace(' ', '_').\
                        replace('-', '_').replace(':', '_')
                    te_type = record.id.split("#")[-1]

                # Check if # is in the seq.id. If # is present, the string before # is the seq_name, and the string after # is the seq_TE_type
                # TODO need to check if assigning seq_name in this way will create duplicates
                else:
                    sanitized_id = record.id.replace('/', '_').replace(' ', '_').replace('-', '_').replace(':', '_')
                    te_type = "Unknown"
                    # modify header to add #Unknown 
                    record.id = f"{record.id}#{te_type}"
                    record.description = record.id

                # Define output file name
                output_filename = os.path.join(output_dir, f"{sanitized_id}.fasta")
                seq_obj = SeqObject(str(sanitized_id), str(output_filename), len(record.seq), te_type)

                # Store all input file information (object manner) to seq_list
                seq_list.append(seq_obj)

                # Write single fasta file with sanitized name
                with open(output_filename, 'w') as output_file:
                    SeqIO.write(record, output_file, 'fasta')

            if detected_pound:
                # TODO, user can disable this annotation
                print("TE Trimmer detects # in your input fasta sequence. The string before # is denoted as the "
                        "seq_name, and the string after # is denoted as the TE type\n")

        print("\nFinish to generate single sequence files.\n")

    elif continue_analysis:

        # When continue_analysis is true, generate seq_list based on single fasta files
        for filename in os.listdir(output_dir):

            file = os.path.join(output_dir, filename)
            with open(file, 'r') as fasta_file:
                for record in SeqIO.parse(fasta_file, 'fasta'):

                    # Get sanitized_id from single fasta file name
                    sanitized_id = os.path.splitext(filename)[0]
                    # Note: single fasta file name is different with record.id
                    te_type = record.id.split("#")[-1]
                    seq_obj = SeqObject(str(sanitized_id), str(file), len(record.seq), te_type)
                    seq_list.append(seq_obj)
        print("\nFinish to read in single sequence files generated from previous analysis.\n")
    return seq_list
