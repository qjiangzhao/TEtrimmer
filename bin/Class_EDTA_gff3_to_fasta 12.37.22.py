from Bio import SeqIO


class EDTAGffSequenceExtractor:
    """
    A class that extract EDTA modified GFF sequences
    """
    def __init__(self, gff_file, genome_file):
        """
        :param gff_file: str, modified EDTA intact gff3 file
        :param genome_file: str, genome file
        """
        self.gff_file = gff_file
        self.genome_file = genome_file
        self.sequences = {}
        self.extract_sequences()

    def extract_sequences(self):
        """
        A method to extract sequences
        """
        with open(self.gff_file, "r") as gff_file:
            for line in gff_file:
                cols = line.strip().split("\t")
                # get rid of lines contain "repeat_region", "target_site_duplication", "long_terminal_repeat"
                if cols[2] in ["repeat_region", "target_site_duplication", "long_terminal_repeat"]:
                    continue
                # store position and strand information to variables
                start, end, strand = int(cols[3]), int(cols[4]), cols[6]

                attrs = cols[8].split(";")
                name_attr = [x for x in attrs if x.startswith("ID=")][0]
                name = name_attr.split("=")[1].split("__")[0] # __ separator can refer to class modify_EDTA_intact
                # read genome fasta file
                genome_seq = SeqIO.parse(self.genome_file, "fasta")
                # extract sequences
                for record in genome_seq:
                    if record.id == cols[0]:
                        if strand == "+":
                            sequence = record.seq[start - 1:end]
                        else:
                            sequence = record.seq[start - 1:end].reverse_complement()
                # store sequences into self.sequences = {} dictionary
                if name in self.sequences:
                    self.sequences[name] += sequence
                else:
                    self.sequences[name] = sequence

    def write_sequences(self, output_file):
        with open(output_file, "w") as output_file:
            for name, sequence in self.sequences.items():
                output_file.write(f">{name}\n{sequence}\n")


"""
test = EDTAGffSequenceExtractor("/work/ur376715/TE_analysis/Second_time_EDTA_result_based_on_bgh_gene_model_v4.2/bgh_dh14_v4.fa.mod.EDTA.intact_modified.gff3",
                                "/home/ur376715/Bgh/Bgh_seq/Genome/bgh_dh14_v4.fa")

test.write_sequences("/work/ur376715/TE_analysis/Second_time_EDTA_result_based_on_bgh_gene_model_v4.2/bgh_dh14_v4.fa.mod.EDTA.intact_modified.fa")
"""