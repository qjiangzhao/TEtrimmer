from Bio import SeqIO

def get_genome_length(fasta_file):
    total_length = sum(len(record.seq) for record in SeqIO.parse(fasta_file, "fasta"))
    return total_length

test = get_genome_length("/Users/panstrugamacbook/Documents/TE_Trimmer/bgh_genome/bgh_dh14_v4.fa")
print(test)