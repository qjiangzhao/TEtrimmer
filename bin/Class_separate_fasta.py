from Bio import SeqIO
import os
from Class_seq_object import SeqObject


class FastaSequenceSeparator:
    """
    Separate fasta file to multiple file, which only contain one sequence.
    """

    def __init__(self, input_fasta, output_dir):
        self.input_fasta = input_fasta
        self.output_dir = output_dir

    def separate_sequences(self, continue_analysis=False):
        """
        separates input file into single fasta file and creates object for each input sequence
        """
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        seq_list = []

        if not continue_analysis:

            print(
                "'/', '-', ':' and ' ' before '#' of input sequence names will be converted to '_'\n")
            detected_pound = False
            with open(self.input_fasta, 'r') as fasta_file:
                for record in SeqIO.parse(fasta_file, 'fasta'):
                    if len(record.id.split("#")) > 1:
                        detected_pound = True
                        sanitized_id = record.id.split("#")[0].replace('/', '_').replace(' ', '_').\
                            replace('-', '_').replace(':', '_')
                        te_type = record.id.split("#")[1]

                    # Check if # is in the seq.id. If # is present, the string before # is the seq_name,
                    # and the string after # is the seq_TE_type
                    # If # isn't found, use Unknown as TE type
                    # TODO need to check if assigning seq_name in this way will create duplicates
                    else:
                        sanitized_id = record.id.replace('/', '_').replace(' ', '_').replace('-', '_').replace(':', '_')
                        te_type = "Unknown"

                    # Define output file name
                    output_filename = os.path.join(self.output_dir, f"{sanitized_id}.fasta")
                    seq_obj = SeqObject(str(sanitized_id), str(output_filename), len(record.seq), te_type)

                    # Store all input file information (object manner) to seq_list
                    seq_list.append(seq_obj)

                    # Write single fasta file with sanitized name
                    with open(output_filename, 'w') as output_file:
                        SeqIO.write(record, output_file, 'fasta')

                if detected_pound:
                    # TODO, user can disable this annotation
                    print("The input sequence name before '#' is denoted as the sequence name, "
                          "and the string after '#' is denoted as the TE type\n")

            print("\nFinish to generate single sequence files.\n")

        elif continue_analysis:

            # When continue_analysis is true, generate seq_list based on single fasta files
            for filename in os.listdir(self.output_dir):

                file = os.path.join(self.output_dir, filename)

                record = SeqIO.read(file, "fasta")

                if len(record.id.split("#")) > 1:
                    sanitized_id = record.id.split("#")[0].replace('/', '_').replace(' ', '_'). \
                        replace('-', '_').replace(':', '_')
                    te_type = record.id.split("#")[-1]

                else:
                    sanitized_id = record.id.replace('/', '_').replace(' ', '_').replace('-', '_').replace(':', '_')
                    te_type = "Unknown"

                seq_obj = SeqObject(str(sanitized_id), str(file), len(record.seq), te_type)
                seq_list.append(seq_obj)

            print("\nFinish to read in single sequence files generated from previous analysis.\n")

        return seq_list
