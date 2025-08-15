import sys
import os
import pyhmmer
import io

import re

import datetime

import collections


class pyhmmer_manager:
    def __init__(self, hmm_models = None, hmm_dat = None, cpus = 1):
        self.model_file = hmm_models
        self.dat_file = hmm_dat

        self.hmm_models = None

        self.threads = cpus

        self.nt_easel = pyhmmer.easel.Alphabet.dna()
        self.amino_easel = pyhmmer.easel.Alphabet.amino()

        self.pfam_operator = None

    def prep_pfam_op(self):
        self.pfam_operator = pfam_filterer(self.dat_file)
        self.pfam_operator.load_pfam_data()

    def load_hmm_database(self):
        print("")
        print("Loading DB")
        ar = open(self.model_file, mode = 'rb')
        hmm_text = ar.read()
        #hmm_io = io.BytesIO(hmm_text.encode(encoding = "ascii"))
        hmm_io = io.BytesIO(hmm_text)

        print("Interpreting models")
        with pyhmmer.plan7.HMMFile(hmm_io) as fh:
            hmms = list(fh)

        #I really don't think we need to do this, but this is WAY faster if we do.
        #print("Creating press db")
        #pyhmmer.hmmer.hmmpress(hmms, "test_pyhmmer_press")

        print("Optimizing HMMs for search")
        hmms = [h.to_profile().to_optimized() for h in hmms]

        self.hmm_models = hmms

        print("HMM search prepared!")
        print("")

    #Function for prepping a sequence for use with pyHMMER. Reads the sequence into memory, then converts it for a search
    def prepare_nucleotide_seq_for_search(self, nucleotide_fasta):
        ar = open(nucleotide_fasta)
        fasta = ar.read()
        ar.close()

        total_seqlen = 0

        seqs = []
        fasta = fasta.split(">")
        fasta = fasta[1:] #The first split is on the first character and returns an empty string at pos 0
        for record in fasta:
            record = record.split("\n")
            seqid = record[0]
            seqid = seqid.encode(encoding = "ascii")
            seq = ''.join(record[1:])

            total_seqlen += len(seq)

            #Pyhmmer digitization of sequences for searching.
            easel_seq = pyhmmer.easel.TextSequence(name = seqid, sequence = seq)
            easel_seq = easel_seq.digitize(self.amino_easel)
            seq = None

            seqs.append(easel_seq)

        return seqs, total_seqlen

    #Actually execute the HMM search in memory
    def execute_search(self, sequences, seq_e_value_cutoff = 1e-2, domain_e_value_cutoff = 1e-2, output_file = None, cutoff_method = "evalue"):

        #Use of Z to set the number of sequences for E-value calculations so it matches the E values produced by an HMMscan.
        #The difference in E values is caused by treating the HMM models as queries in the hmmscan and as targets in the hmmsearch (or maybe backwards from that)
        #Alignments produced are identical, only E values differ.
        if cutoff_method == "evalue":
            res = pyhmmer.hmmsearch(self.hmm_models, sequences, cpus = 1, E=seq_e_value_cutoff, domE=domain_e_value_cutoff, Z=len(self.hmm_models))
        #Gathering thresholds are manually set bitscore cutoffs assigned by the curator of each pfam model; they are generally conservative
        if cutoff_method == "gathering":
            res = pyhmmer.hmmsearch(self.hmm_models, sequences, bit_cutoffs = "gathering", cpus = 1, Z=len(self.hmm_models))
        #Trusted cutoffs are a calculated value representing the lowest bitscore that filters out all false positive hits for an HMM model within pfams data
        if cutoff_method == "trusted":
            res = pyhmmer.hmmsearch(self.hmm_models, sequences, bit_cutoffs = "trusted", cpus = 1, Z=len(self.hmm_models))

        #The above options should be increasingly conservative.

        #Cleanup.
        res = [r for r in res if len(r) > 0]


        formatted = self.pfam_operator.parse_hmmsearch_output(res)
        cleaned = self.pfam_operator.resolve_overlapping_domains(formatted)

        #Check if there's any results
        has_any_domains = False
        for c in cleaned:
            if len(cleaned[c]) > 0:
                has_any_domains = True
                break

        self.pfam_operator.pretty_print(cleaned, output_file)

        return has_any_domains

    def prepare(self):
        self.prep_pfam_op()
        self.load_hmm_database()

    def run(self, args):
        infile, outfile = args[0], args[1]
        seqs, seqlen = self.prepare_nucleotide_seq_for_search(infile)
        #print("HMMsearching", infile, "start", seqlen, "bp")
        domains_detected = self.execute_search(sequences = seqs, output_file = outfile)
        #print("HMMsearching", infile, "end")

        return (outfile, domains_detected,)



class pfam_filterer:
    def __init__(self, pfam_data_file):
        self.pfam_file = pfam_data_file
        self.pfam_data = None
        self.load_pfam_data()
        self.te_name_regex = re.compile(r'(.+)\[(\d+)-(\d+)]$')


    #The .dat file here is a separate download. We only use it for the clan data
    def load_pfam_data(self):
        #Function for collecting annotation data from a pfam database
        """Reads the Pfam data file to dictionary.

        Args:
            filename: Name/Path of the Pfam data file (Pfam-A.hmm.dat).

        Returns:
            A dict mapping HMM profile name to the corresponding information.
            For example:

            {'1-cysPrx_C': Data(type='Domain', clan=None, ga_seq=21.1, ga_dom=21.1),
             'RRM': Data(type='Domain', clan=None, ga_seq=21.0, ga_dom=21.0),
             'SOXp': Data(type='Family', clan=None, ga_seq=22.1, ga_dom=22.1)}
        """
        self.pfam_data = {}
        Data = collections.namedtuple('Data', ['type', 'clan', 'ga_seq', 'ga_dom'])
        with open(self.pfam_file) as fh:
            clan = None   # Not all domains have clan assigned.
            for line in fh:
                if line.startswith('#=GF ID'):
                    hmm_name = line[10:-1]
                elif line.startswith('#=GF TP'):
                    typ = line[10:-1]
                elif line.startswith('#=GF CL'):
                    clan = line[10:-1]
                elif line.startswith('#=GF GA'):
                    scores = line[10:-1].strip().rstrip(';').split(';')
                    ga_seq = float(scores[0])
                    ga_dom = float(scores[1])
                elif line.startswith('//'):
                    self.pfam_data[hmm_name] = Data(typ, clan, ga_seq, ga_dom)
                    clan = None


    #TE_name	orf_start	orf_end	domain_start	domain_end	direction	domain_name	domain_reference
    #rnd_1_family_470_1	261	1628	861	1082	+	Retrotrans_gag	PF03732.24
    def parse_hmmsearch_output(self, hmmsearch_records):
        DOMAIN_INFO = [
            'te_name',
            'te_start',
            'te_end',
            'aln_start',
            'aln_end',
            'domain_start',
            'domain_end',
            'hmm_acc',
            'hmm_name',
            'direction',
            'evalue',
            'clan'
        ]
        Domain = collections.namedtuple('Domain', DOMAIN_INFO)

        results = {}
        for tophits in hmmsearch_records:
            hmm_accession = tophits.query.accession.decode(encoding="ascii")
            hmm_name = tophits.query.name.decode(encoding="ascii")
            clan = self.pfam_data[hmm_name].clan

            for hit in tophits:
                name_long = hit.name.decode(encoding="ascii")
                regex_match = self.te_name_regex.search(name_long).groups()
                seq_name = regex_match[0]
                te_start = int(regex_match[1])
                te_end = int(regex_match[2])
                align_from = hit.best_domain.alignment.target_from
                align_to = hit.best_domain.alignment.target_to

                evalue = hit.evalue

                if te_end > te_start:
                    direction = "+"
                    domain_start = te_start + (3*align_from)-3
                    domain_end = te_start + (3*align_to) - 1
                else:
                    direction = "-"
                    domain_start = te_start - (3*align_to)+1
                    domain_end = te_start - (3*align_from)+3

                #print(seq_name, te_start, te_end, domain_start, domain_end, direction, hmm_name, clan, hit.evalue, align_from, align_to)

                dom = Domain(seq_name,
                            te_start,
                            te_end,
                            align_from,
                            align_to,
                            domain_start,
                            domain_end,
                            hmm_accession,
                            hmm_name,
                            direction,
                            evalue,
                            clan
                        )

                if seq_name not in results:
                    results[seq_name] = []
                results[seq_name].append(dom)

        # For each protein, sort domains by their start position.
        for seq_name in results:
            results[seq_name].sort(key=lambda x: x.aln_start)

        return results

    def resolve_overlapping_domains(self, results):
        """Resolves overlapping domains belonging to the same clan.

        When a protein sequence region has overlapping matches to more than one
        domains within the same clan, only the best scoring domain match is shown.

        Args:
            results: a dict with the domain results
        """

        # TODO: Consider resolving overlapping domains belonging to different clans.
        for seq_id, domains in results.items():
            is_overlap = True
            while is_overlap:
                is_overlap = False
                indexes = set()
                for i in range(len(domains)-1):
                    j = i + 1
                    di = domains[i]  # Domain i (previous domain)
                    dj = domains[j]  # Domain j (next domain)
                    if di.clan:
                        # Overlapping domains within the same clan
                        if di.clan == dj.clan and dj.aln_start <= di.aln_end:
                            idx = j if di.evalue < dj.evalue else i
                            indexes.add(idx)
                            is_overlap = True
                # Filter out overlapping domains that have weak scores.
                domains = [d for i, d in enumerate(domains) if i not in indexes]
            results[seq_id] = domains

        return results

    def pretty_print(self, results, outfile = None):
        with open(outfile, "w") as out:
            header = ['TE_name', 'orf_start', 'orf_end', 'domain_start', 'domain_end', 'direction', 'domain_name', 'domain_reference']
            print(*header, sep = "\t", file = out)
            for seq in sorted(results):
                for domain in results[seq]:
                    if domain.direction == "+":
                        print(seq,
                            domain.te_start,
                            domain.te_end,
                            domain.domain_start,
                            domain.domain_end,
                            domain.direction,
                            domain.hmm_name,
                            domain.hmm_acc,
                            sep= "\t",
                            file = out)
                    else:
                        print(seq,
                            domain.te_end,
                            domain.te_start,
                            domain.domain_start,
                            domain.domain_end,
                            domain.direction,
                            domain.hmm_name,
                            domain.hmm_acc,
                            sep= "\t",
                            file = out)


'''
db = "TEtrimmer-main/TEtrimmer-main/pfam_database/Pfam-A.hmm"
dat = "TEtrimmer-main/TEtrimmer-main/pfam_database/Pfam-A.hmm.dat"

dir_base = "TE_OG_out/Multiple_sequence_alignment/"
orfs = os.listdir(dir_base)
orfs = [dir_base + o for o in orfs if o.endswith(".fasta_orfm.txt")]


thds = int(sys.argv[1])
mn = pyhmmer_manager(threads = thds)
mn.prep_pfam_op(dat)

mn.load_hmm_database(db)

for o in orfs:
    print("Searching", o)
    #start = datetime.datetime.now()

    seq = mn.prepare_nucleotide_seq_for_search(o)
    printout = mn.execute_search(seq)
    end = datetime.datetime.now()
    #print("Search took", (end-start).seconds, "seconds")
    #print("Search took", (end-start).microseconds, "microseconds")
    print("")


#check_outputs actually used by TEtrimmer
#rnd_1_family_188.fasta_orfm_pfm.txt
outs = os.listdir(dir_base)
outs = [dir_base + o for o in outs if o.endswith(".fasta_orfm_pfm.txt")]
'''
