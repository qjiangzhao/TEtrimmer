import os


class Seq_object:
    """
    create object for each input sequence.
    """

    def __init__(self, name, path_input_fasta_file, length, TE_type):
        self.name = name
        self.input_fasta = path_input_fasta_file
        self.old_length = length
        self.old_TE_type = TE_type
        self.low_copy = False
        self.consi_obj_list = []  # an input seq can end up with multiple consensus seq
        self.blast_hit_n = 0
        self.status = "unprocessed"  # "unprocessed","processed", "skipped"

    def get_type(self):
        return self.old_TE_type
    
    def get_length(self):
        return self.old_length
    
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
            # input_name, consensus_name, old_length, blast_hit_n, low_copy, old_TE_type, status, new_length, new_TE_type
            if len(self.consi_obj_list) > 0:
                for consi_obj in self.consi_obj_list:
                    TE_type = consi_obj.get_TE_type_for_file()
                    f.write(f"{str(self.name)},{str(consi_obj.consensus_name)}#{str(TE_type)},{str(self.old_length)},{str(self.blast_hit_n)},{str(self.low_copy)},{str(self.old_TE_type)},{str(self.status)},{str(consi_obj.new_length)},{str(consi_obj.new_TE_type)}\n")
            else:
                f.write(f"{str(self.name)},N/A,{str(self.old_length)},{str(self.blast_hit_n)},{str(self.low_copy)},{str(self.old_TE_type)},{str(self.status)},N/A,N/A\n")

    def update_low_copy(self, check_80, found_match):
        if check_80 and found_match:
            self.low_copy = True
        return self.low_copy

    def update_blast_hit_n(self, blast_hit_n):
        self.blast_hit_n = blast_hit_n


class ConsensusObject:
    def __init__(self, parent_seq_object, consensus_name):
        self.parent_seq_object = parent_seq_object
        self.consensus_name = consensus_name
        self.new_length = "N/A"
        self.new_TE_type = "N/A"
    
    def set_new_lenth(self, new_length):
        self.new_length = new_length
    
    def set_new_TE_type(self, new_TE_type):
        self.new_TE_type = new_TE_type
    
    def get_TE_type_for_file(self):

        # For writing the final consensus file, the TE_type can be either the input TE type or the new
        # classified one if it is newly classified
        if "N/A" in self.new_TE_type or "Unknow" in self.new_TE_type:
            return self.parent_seq_object.get_type()
        else:
            return self.new_TE_type