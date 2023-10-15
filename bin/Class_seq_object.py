class SeqObject:

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
        self.old_terminal_repeat = "None"
        self.old_blast_full_n = "None"

    def get_seq_name(self):
        return self.name
    
    def get_old_TE_type(self):
        return self.old_TE_type
    
    def get_length(self):
        return self.old_length
    
    def get_input_fasta(self):
        return self.input_fasta

    def set_old_terminal_repeat(self, terminal_repeat):
        self.old_terminal_repeat = terminal_repeat

    def set_old_blast_full_n(self, blast_full_length_n):
        self.old_blast_full_n = blast_full_length_n

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

            # "input_name,consensus_name, blast_hit_n, cons_MSA_seq_n, cons_full_blast_n, input_length, cons_length, "
            # "input_TE_type, reclassified_type, cons_flank_repeat, evaluation, low_copy, status\n"
            if len(self.consi_obj_list) > 0:

                for consi_obj in self.consi_obj_list:

                    f.write(f"{str(self.name)},{str(consi_obj.consensus_name)},"  # name
                            f"{str(self.blast_hit_n)},{str(consi_obj.new_TE_MSA_seq_n)},"  # sequence number
                            f"{str(consi_obj.new_TE_blast_full_length_n)},"  # blast full length number
                            f"{str(self.old_length)},{str(consi_obj.new_length)},"  # sequence length
                            f"{str(self.old_TE_type)},{str(consi_obj.get_TE_type_for_file())},"  # TE type
                            f"{str(consi_obj.new_TE_terminal_repeat)},{str(self.low_copy)},"
                            f"{str(consi_obj.get_evaluation())},{str(self.status)}\n")  # Evaluation

            elif self.low_copy:

                f.write(f"{str(self.name)},{str(self.name)},"  # name
                        f"{str(self.blast_hit_n)},NaN,"  # sequence number
                        f"{str(self.old_blast_full_n)},"  # blast full length number for low copy elements
                        f"{str(self.old_length)},{str(self.old_length)},"  # sequence length
                        f"{str(self.old_TE_type)},{str(self.old_TE_type)},"  # TE type
                        f"{str(self.old_terminal_repeat)},{str(self.low_copy)},"
                        f"NaN,{str(self.status)}\n")
            else:

                f.write(f"{str(self.name)},NaN,"  # name
                        f"{str(self.blast_hit_n)},NaN,"
                        f"NaN,"  # sequence number
                        f"{str(self.old_length)},NaN,"  # sequence length
                        f"{str(self.old_TE_type)},NaN,"  # TE type
                        f"NaN,{str(self.low_copy)},"
                        f"NaN,{str(self.status)}\n")

    def update_low_copy(self, check_blast, found_match):
        if check_blast and found_match:
            self.low_copy = True
        return self.low_copy

    def update_blast_hit_n(self, blast_hit_n):
        self.blast_hit_n = blast_hit_n


class ConsensusObject:
    def __init__(self, parent_seq_object, consensus_name):
        self.parent_seq_object = parent_seq_object
        self.consensus_name = str(consensus_name)
        self.proof_annotation_file = "None"
        self.hmm_file = "None"
        self.proof_pdf = "None"
        self.proof_fasta = "None"
        self.proof_anno = "None"
        self.new_length = "NaN"
        self.new_TE_type = "NaN"
        self.new_TE_MSA_seq_n = "NaN"
        self.new_TE_terminal_repeat = "None"
        self.new_TE_blast_full_length_n = "NaN"
        self.cons_seq = "NaN"
        self.cons_pfam = False
        self.cons_evaluation = "Need_check"

    def set_new_length(self, new_length):
        self.new_length = int(new_length)
    
    def set_new_TE_type(self, new_TE_type):
        self.new_TE_type = new_TE_type

    # Store multiple sequence alignment sequence number to object
    def set_cons_MSA_n(self, new_TE_MSA_seq_n):
        self.new_TE_MSA_seq_n = int(new_TE_MSA_seq_n)

    # Store flank type ("LTR" or "TIR") to object
    def set_new_terminal_repeat(self, terminal_repeat):
        self.new_TE_terminal_repeat = terminal_repeat

    # Store consensus sequence to object
    def set_cons_seq(self, cons_seq):
        self.cons_seq = cons_seq

    def set_cons_pfam(self, if_pfam):
        self.cons_pfam = if_pfam

    def set_cons_evaluation(self, level):
        self.cons_evaluation = str(level)

    def get_evaluation(self):
        return self.cons_evaluation

    def get_consi_name(self):
        return self.consensus_name
    
    def get_new_TE_type(self):
        return self.new_TE_type

    # Store full length blast number
    def set_blast_full_n(self, blast_full_length_n):
        self.new_TE_blast_full_length_n = int(blast_full_length_n)

    def set_proof_annotation_file(self):

        proof_TE_type = self.new_TE_type.replace("/", "__")
        self.proof_pdf = f"{self.consensus_name}#{proof_TE_type}.pdf"
        self.proof_fasta = f"{self.consensus_name}#{proof_TE_type}.fa"
        self.proof_anno = f"{self.consensus_name}#{proof_TE_type}.anno.fa"

    def set_hmm_file(self):

        proof_TE_type = self.new_TE_type.replace("/", "__")
        self.hmm_file = f"{self.consensus_name}#{proof_TE_type}.hmm"

    def get_TE_type_for_file(self):

        # For writing the final consensus file, the TE_type can be either the input TE type or the new
        # classified one if it is newly classified
        if "NaN" in self.new_TE_type or "unknown" == self.new_TE_type.lower():
            return self.parent_seq_object.old_TE_type

        else:
            return self.new_TE_type
