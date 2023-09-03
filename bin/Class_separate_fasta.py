from Bio import SeqIO
import os
from Class_seq_object import Seq_object


class FastaSequenceSeparator:
    """
    Separate fasta file to multiple file, which only contain one sequence.
    """

    def __init__(self, input_fasta, output_dir):
        self.input_fasta = input_fasta
        self.output_dir = output_dir

    def separate_sequences(self):
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        print("TE Trimmer is modifying sequence names, '/', '#', '-', ':' and empty space will be converted to '_'")

        with open(self.input_fasta, 'r') as fasta_file:
            for record in SeqIO.parse(fasta_file, 'fasta'):
                sanitized_id = record.id.replace('/', '_').replace(' ', '_').replace('#', '_').\
                    replace('-', '_').replace(':', '_')
                output_filename = os.path.join(self.output_dir, f"{sanitized_id}.fasta")
                with open(output_filename, 'w') as output_file:
                    SeqIO.write(record, output_file, 'fasta')

        print("Finish to generate single sequence files.")
    
    def create_seq_obj(self):
        seq_dict = {}

        with open(self.input_fasta, 'r') as fasta_file:
            print("TE Trimmer is modifying sequence names, '/', '-', ':' and empty space will be converted to '_'")
            detected_pound = False
            for record in SeqIO.parse(fasta_file, 'fasta'):
                # check if # is in the seq.id. If # is present, the string before # is the seq_name, and the string after # is the seq_TE_type
                # TODO need to check if assigning seq_name in this way will create duplicates
                if len(record.id.split("#")) > 0:
                    detected_pound = True
                    sanitized_id = record.id.split("#")[0].replace('/', '_').replace(' ', '_').replace('#', '_').\
                    replace('-', '_').replace(':', '_')
                    te_type = record.id.split("#")[-1]
                else:
                    sanitized_id = record.id.replace('/', '_').replace(' ', '_').replace('#', '_').\
                    replace('-', '_').replace(':', '_')
                    te_type = "Unknown"
                seq_obj = Seq_object(str(sanitized_id), str(record.seq),len(record.seq), te_type)
                seq_dict[sanitized_id] = seq_obj
            if detected_pound:
                print("TE Trimmer detects # in your input fasta sequence. The string before # is denoted as the seq_name, and the string after # is denoted as the TE type")
        return seq_dict