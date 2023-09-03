import os

class Seq_object:
    """
    create object corresponding to one sequence.
    """

    def __init__(self, name, input_fasta, length, TE_type):
        self.name = name
        self.input_fasta = input_fasta
        self.old_length = length
        self.new_length = length
        self.TE_type = TE_type
        self.single_copy = False
        self.status = "not_processed" # "not_processed","processed", "skipped"

    def get_type(self):
        return self.TE_type
    
    def change_in_length(self):
        return self.new_length - self.old_length
    
    def check_unknown(self):
        if "unknown" in self.TE_type.lower():
            return True
        else:
            return False