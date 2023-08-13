from Bio import SeqIO
import os


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

        with open(self.input_fasta, 'r') as fasta_file:
            for record in SeqIO.parse(fasta_file, 'fasta'):
                sanitized_id = record.id.replace('/', '_').replace(' ', '_').replace('#', '_')
                output_filename = os.path.join(self.output_dir, f"{sanitized_id}.fasta")
                with open(output_filename, 'w') as output_file:
                    SeqIO.write(record, output_file, 'fasta')