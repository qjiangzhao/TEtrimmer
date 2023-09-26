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

    def separate_sequences(self, continue_analysis = False):
        """
        separates input file into single fasta file and creates object for each input sequence
        """
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        seq_list = []
        if not continue_analysis:
            print("TE Trimmer is modifying sequence names, '/', '-', ':' and empty space will be converted to '_'")
            detected_pound = False
            with open(self.input_fasta, 'r') as fasta_file:
                for record in SeqIO.parse(fasta_file, 'fasta'):
                    if len(record.id.split("#")) > 1:
                        detected_pound = True
                        sanitized_id = record.id.split("#")[0].replace('/', '_').replace(' ', '_').replace('#', '_').\
                        replace('-', '_').replace(':', '_')
                        te_type = record.id.split("#")[-1]
                    # check if # is in the seq.id. If # is present, the string before # is the seq_name, and the string after # is the seq_TE_type
                    # TODO need to check if assigning seq_name in this way will create duplicates
                    else:
                        sanitized_id = record.id.replace('/', '_').replace(' ', '_').replace('#', '_').\
                        replace('-', '_').replace(':', '_')
                        te_type = "Unknown"
                    output_filename = os.path.join(self.output_dir, f"{sanitized_id}.fasta")
                    seq_obj = Seq_object(str(sanitized_id), str(output_filename),len(record.seq), te_type)
                    seq_list.append(seq_obj)
                    with open(output_filename, 'w') as output_file:
                        SeqIO.write(record, output_file, 'fasta')
                if detected_pound:
                    # TODO, user can disable this annotation
                    print("TE Trimmer detects # in your input fasta sequence. The string before # is denoted as the seq_name, and the string after # is denoted as the TE type")
            print("Finish to generate single sequence files.")
        elif continue_analysis:
            for filename in os.listdir(self.output_dir):
                file = os.path.join(self.output_dir, filename)
                with open(file, 'r') as fasta_file:
                    for record in SeqIO.parse(fasta_file, 'fasta'):
                            sanitized_id = filename.split(".")[0]
                            te_type = record.id.split("#")[-1]
                            seq_obj = Seq_object(str(sanitized_id), str(file),len(record.seq), te_type)
                            seq_list.append(seq_obj)
            print("Finish to read in single sequence files generated from previous analysis.")
        return seq_list