import os


class Seq_object:

    """
    create object for each input sequence.
    """

    def __init__(self, name, path_input_fasta_file, length, TE_type):
        self.name = str(name)
        self.input_fasta = str(path_input_fasta_file)
        self.old_length = int(length)
        self.old_TE_type = str(TE_type)
        self.low_copy = False
        self.consi_obj_list = []  # an input seq can end up with multiple consensus seq
        self.blast_hit_n = 0
        self.status = "unprocessed"  # "unprocessed","processed", "skipped"

    def get_seq_name(self):
        return self.name
    
    def get_old_TE_type(self):
        return self.old_TE_type
    
    def get_length(self):
        return self.old_length
    
    def get_input_fasta(self):
        return self.input_fasta
    
    def check_unknown(self):
        if "unknown" in self.old_TE_type.lower():
            return True
        else:
            return False
    
    def create_consi_obj(self, consi_name):
        consi_obj = ConsensusObject(self, consi_name)  # Pass the parent seq_object as an argument
        self.consi_obj_list.append(consi_obj)
        return consi_obj

    # update_status function will write object information to progress file, use it when the analysis is finished
    def update_status(self, new_status, progress_file):

        self.status = new_status
        with open(progress_file, "a") as f:

            # "input_name,consensus_name, blast_hit_n, cons_MSA_seq_n, input_length, cons_length, "
            # "input_TE_type, reclassified_type, cons_flank_repeat, low_copy, status\n"
            if len(self.consi_obj_list) > 0:

                for consi_obj in self.consi_obj_list:
                    TE_type = consi_obj.get_TE_type_for_file()
                    f.write(f"{str(self.name)},{str(consi_obj.consensus_name)}#{str(TE_type)},"  # name
                            f"{str(self.blast_hit_n)},{str(consi_obj.new_TE_MSA_seq_n)},"  # sequence number
                            f"{str(self.old_length)},{str(consi_obj.new_length)},"  # sequence length
                            f"{str(self.old_TE_type)},{str(consi_obj.get_TE_type_for_file())},"  # TE type
                            f"{str(consi_obj.new_TE_flank_repeat)},{str(self.low_copy)},{str(self.status)}\n")
            else:

                f.write(f"{str(self.name)},NaN,"  # name
                        f"{str(self.blast_hit_n)},NaN,"  # sequence number
                        f"{str(self.old_length)},NaN,"  # sequence length
                        f"{str(self.old_TE_type)},NaN,"  # TE type
                        f"NaN,{str(self.low_copy)},{str(self.status)}\n")

    def update_low_copy(self, low_copy):
        self.low_copy = low_copy
        return self.low_copy

    def update_blast_hit_n(self, blast_hit_n):
        self.blast_hit_n = blast_hit_n


class ConsensusObject:
    def __init__(self, parent_seq_object, consensus_name):
        self.parent_seq_object = parent_seq_object
        self.consensus_name = str(consensus_name)
        self.new_length = "NaN"
        self.new_TE_type = "Unknown"
        self.new_TE_MSA_seq_n = "NaN"
        self.new_TE_flank_repeat = "None"
        self.cons_seq = "NaN"
    
    def set_new_lenth(self, new_length):
        self.new_length = new_length
    
    def set_new_TE_type(self, new_TE_type):
        self.new_TE_type = new_TE_type

    # Store multiple sequence alignment sequence number to object
    def set_cons_MSA_n(self, new_TE_MSA_seq_n):
        self.new_TE_MSA_seq_n = new_TE_MSA_seq_n

    # Store flank type ("LTR" or "TIR") to object
    def set_new_flank_repeat(self, flank_repeat):
        self.new_TE_flank_repeat = flank_repeat

    # Store consensus sequence to object
    def set_cons_seq(self, cons_seq):
        self.cons_seq = cons_seq
    
    def get_consi_name(self):
        return self.consensus_name
    
    def get_new_TE_type(self):
        return self.new_TE_type
    
    def get_TE_type_for_file(self):
        # For writing the final consensus file, the TE_type can be either the input TE type or the new
        # classified one if it is newly classified
        if "Unknown" in self.new_TE_type:
            return self.parent_seq_object.get_old_TE_type()
        else:
            return self.new_TE_type