import traceback
import logging

class SeqObject:
    """
    Create object for each input sequence.
    """

    def __init__(self, name, path_input_fasta_file, length, TE_type):
        self.name = str(name)
        self.input_fasta = str(path_input_fasta_file)
        self.old_length = int(length)
        self.old_TE_type = str(TE_type)
        self.low_copy = False
        self.consi_obj_list = []  # an input sequence can result in multiple consensus sequences
        self.blast_hit_n = 0
        self.status = 'unprocessed'  # "unprocessed", "processed", "skipped"
        self.old_terminal_repeat = 'NaN'
        self.old_blast_full_n = 'NaN'
        self.genome_length = 'NaN'
        self.input_genome_cov_len = 'NaN'

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

    def set_input_genome_cov_len(self, input_genome_cov_len):
        self.input_genome_cov_len = input_genome_cov_len


    def check_unknown(self):
        if 'unknown' in self.old_TE_type.lower():
            return True
        else:
            return False

    def create_consi_obj(self, consi_name):
        consi_obj = ConsensusObject(
            self, consi_name
        )  # Pass the parent seq_object as an argument
        self.consi_obj_list.append(consi_obj)
        return consi_obj

    def update_low_copy(self, check_blast, found_match):
        if check_blast and found_match:
            self.low_copy = True
        return self.low_copy

    def update_blast_hit_n(self, blast_hit_n):
        self.blast_hit_n = blast_hit_n


    # update_status function will write object information to progress file to be used when the analysis is complete
    def update_status(self, new_status, progress_file):
        try:
            self.status = new_status
            with open(progress_file, 'a') as f:
                if len(self.consi_obj_list) > 0:  # This doesn't include skipped and low_copy instances
                    for consi_obj in self.consi_obj_list:
                        f.write(
                            f'{str(self.name)},'  # input_name
                            f'{str(consi_obj.consensus_name)},'  # output_name
                            f'{str(self.blast_hit_n)},'  # input_blast_n
                            f'{str(self.old_blast_full_n)},'  # input_full_blast_n
                            f'{str(consi_obj.cons_blast_n)},'  # output_blast_n
                            f'{str(consi_obj.new_TE_blast_full_length_n)},'  # output_full_blast_n
                            f'{str(self.input_genome_cov_len)},'  # input_genome_cov_len
                            f'{str(consi_obj.output_genome_cov_len)},'  # output_genome_cov_len
                            f'{str(int(consi_obj.new_TE_MSA_seq_n))},'  # output_MSA_seq_n
                            f'{str(self.old_length)},'  # input_length
                            f'{str(consi_obj.new_length)},'  # output_length
                            f'{str(consi_obj.in_out_identity)},'  # identity
                            f'{str(consi_obj.input_coverage)},'  #  input_coverage
                            f'{str(consi_obj.output_coverage)},'  # output_coverage
                            f'{str(self.old_TE_type)},'  # input_TE_type
                            f'{str(consi_obj.get_TE_type_for_file())},'  # output_TE_type
                            f'{str(self.old_terminal_repeat)},'  # input_terminal_repeat
                            f'{str(consi_obj.new_TE_terminal_repeat)},'  # output_terminal_repeat
                            f'{str(self.low_copy)},'  # low_copy
                            f'{str(consi_obj.cons_tsd)},'  # TSD
                            f'{str(consi_obj.get_evaluation())},'  # evaluation
                            f'{str(self.status)}\n'  # status
                        )

                else:
                    f.write(
                        f'{str(self.name)},'  # input_name
                        f'{str(self.name)},'  # output_name
                        f'{str(self.blast_hit_n)},'  # input_blast_n
                        f'{str(self.old_blast_full_n)},'  # input_full_blast_n
                        f'{str(self.blast_hit_n)},'  # output_blast_n
                        f'{str(self.old_blast_full_n)},'  # output_full_blast_n
                        f'{str(self.input_genome_cov_len)},'  # input_genome_cov_len
                        f'NaN,'  # output_genome_cov_len
                        f'NaN,'  # output_MSA_seq_n
                        f'{str(self.old_length)},'  # input_length
                        f'{str(self.old_length)},'  # output_length
                        f'NaN,'  # identity
                        f'NaN,'  #  input_coverage
                        f'NaN,'  # output_coverage
                        f'{str(self.old_TE_type)},'  # input_TE_type
                        f'{str(self.old_TE_type)},'  # output_TE_type
                        f'{str(self.old_terminal_repeat)},'  # input_terminal_repeat
                        f'NaN,'  # output_terminal_repeat
                        f'{str(self.low_copy)},'  # low_copy
                        f'NaN,'  # TSD
                        f'NaN,'  # evaluation
                        f'{str(self.status)}\n'  # status
                    )
        except Exception as e:
            tb_content = traceback.format_exc()
            logging.error(f'seqclass error:\n{e}\n{tb_content}\n')



class ConsensusObject:
    def __init__(self, parent_seq_object, consensus_name):
        self.parent_seq_object = parent_seq_object
        self.consensus_name = str(consensus_name)
        self.proof_curation_file = 'None'
        self.hmm_file = 'None'
        self.proof_pdf = 'None'
        self.proof_fasta = 'None'
        self.proof_raw = 'None'
        self.proof_cluster = 'None'
        self.new_length = 'NaN'
        self.new_TE_type = 'NaN'
        self.new_TE_MSA_seq_n = 'NaN'
        self.new_TE_terminal_repeat = 'False'
        self.new_TE_blast_full_length_n = 'NaN'
        self.cons_seq = 'NaN'
        self.cons_pfam = False
        self.cons_evaluation = 'Need_check'
        self.cons_tsd = 'False'
        self.output_genome_cov_len = 'NaN'
        self.cons_blast_n = 'NaN'
        self.in_out_identity = 'NaN'
        self.input_coverage = 'NaN'
        self.output_coverage = 'NaN'

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

    # Store number of full-length BLAST hits
    def set_blast_full_n(self, blast_full_length_n):
        self.new_TE_blast_full_length_n = int(blast_full_length_n)

    # Store consensus sequence to object
    def set_cons_seq(self, cons_seq):
        self.cons_seq = cons_seq

    def set_cons_pfam(self, if_pfam):
        self.cons_pfam = if_pfam

    def set_cons_evaluation(self, level):
        self.cons_evaluation = str(level)

    def set_tsd(self, if_tsd):
        self.cons_tsd = if_tsd

    def set_output_genome_cov_len(self, output_genome_cov_len):
        self.output_genome_cov_len = output_genome_cov_len

    def set_cons_blast_n(self, cons_blast_n):
        self.cons_blast_n = cons_blast_n

    def set_in_out_identity(self, in_out_identity):
        self.in_out_identity = in_out_identity

    def set_input_coverage(self, input_coverage):
        self.input_coverage = input_coverage

    def set_output_coverage(self, output_coverage):
        self.output_coverage = output_coverage


    # Defining get functions
    def get_tsd(self):
        return self.cons_tsd

    def get_cons_blast_n(self):
        return self.cons_blast_n

    def get_evaluation(self):
        return self.cons_evaluation

    def get_consi_name(self):
        return self.consensus_name

    def get_new_TE_type(self):
        return self.new_TE_type

    def set_proof_curation_file(self):
        proof_TE_type = self.new_TE_type.replace('/', '__')
        self.proof_pdf = f'{self.consensus_name}#{proof_TE_type}.pdf'
        self.proof_fasta = f'{self.consensus_name}#{proof_TE_type}.fa'
        self.proof_raw = f'{self.consensus_name}#{proof_TE_type}.raw.fa'
        self.proof_cluster = f'{self.consensus_name}#{proof_TE_type}.cluster.fa'

    def set_hmm_file(self):
        proof_TE_type = self.new_TE_type.replace('/', '__')
        self.hmm_file = f'{self.consensus_name}#{proof_TE_type}.hmm'

    def get_TE_type_for_file(self):
        # For writing the final consensus file, the TE_type can be either the input TE type or the new
        # type if it has been newly classified
        if 'NaN' in self.new_TE_type or 'unknown' == self.new_TE_type.lower():
            return self.parent_seq_object.old_TE_type

        else:
            return self.new_TE_type
