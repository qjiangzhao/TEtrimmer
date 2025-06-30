import sqlite3
import sys
import os

import pyfastx

import polars

from functions import blast_to_database
from seqclass import SeqObject

#Just let me access values by name
class sequence_params:
	def __init__(self, seq, genome_file, MSA_dir, min_blast_len, min_seq_num, max_msa_lines,
		 top_msa_lines, max_cluster_num, cons_thr, ext_thr, ext_step, classification_dir,
		 max_ext, gap_thr, gap_nul_thr, crop_end_div_thr, crop_end_div_win, crop_end_gap_thr, crop_end_gap_win,
		 start_patterns, end_patterns, output_dir, pfam_dir, mini_orf, single_fasta_n, hmm, hmm_dir,
		 ext_check_win, debug, progress_file, classify_unknown, classify_all,
		 final_con_file, final_con_file_no_low_copy, final_unknown_con_file, final_classified_con_file, low_copy_dir,
		 fast_mode, error_files, plot_skip, skipped_dir, plot_query, engine, proof_curation_dir, poly_patterns, poly_len,
		 perfect_seq_num):
		 
		self.seq = seq
		self.genome_file = genome_file
		self.MSA_dir = MSA_dir
		self.min_blast_len = min_blast_len
		self.min_seq_num = min_seq_num
		self.max_msa_lines = max_msa_lines
		self.top_msa_lines = top_msa_lines
		self.max_cluster_num = max_cluster_num
		self.cons_thr = cons_thr
		self.ext_thr = ext_thr
		self.ext_step = ext_step
		self.classification_dir = classification_dir
		self.max_ext = max_ext
		self.gap_thr = gap_thr
		self.gap_nul_thr = gap_nul_thr
		self.crop_end_div_thr = crop_end_div_thr
		self.crop_end_div_win = crop_end_div_win
		self.crop_end_gap_thr = crop_end_gap_thr
		self.crop_end_gap_win = crop_end_gap_win
		self.start_patterns = start_patterns
		self.end_patterns = end_patterns
		self.output_dir = output_dir
		self.pfam_dir = pfam_dir
		self.mini_orf = mini_orf
		self.single_fasta_n = single_fasta_n
		self.hmm = hmm
		self.hmm_dir = hmm_dir
		self.ext_check_win = ext_check_win
		self.debug = debug
		self.progress_file = progress_file
		self.classify_unknown = classify_unknown
		self.classify_all = classify_all
		self.final_con_file = final_con_file
		self.final_con_file_no_low_copy = final_con_file_no_low_copy
		self.final_unknown_con_file = final_unknown_con_file
		self.final_classified_con_file = final_classified_con_file
		self.low_copy_dir = low_copy_dir
		self.fast_mode = fast_mode
		self.error_files = error_files
		self.plot_skip = plot_skip
		self.skipped_dir = skipped_dir
		self.plot_query = plot_query
		self.engine = engine
		self.proof_curation_dir = proof_curation_dir 
		self.poly_patterns = poly_patterns
		self.poly_len = poly_len
		self.perfect_seq_num = perfect_seq_num

class database:
	def __init__(self, dbpath):
		self.path = dbpath
		self.exists = os.path.exists(dbpath)
		self.is_ok = False
		
		self.conn = None
		self.curs = None
		
	
	def open(self):
		if self.conn is None:
			self.conn = sqlite3.connect(self.path)
			self.curs = self.conn.cursor()
		else:
			print("Database connection is already open!")
		
	def close(self):
		if self.conn is not None:
			self.curs.close()
			self.conn.close()
			self.conn = None
			self.curs = None
		else:
			print("There is no open database to close.")
		
	def initialize(self):
		#The way this is going to be used is one sequence per fasta but thats also ugly
		
		grname = "genome_reference"
		
		genome_records = ', '.join([
						"subject_sequence_name TEXT PRIMARY KEY",
						"subject_sequence_id INTEGER",
						"subject_sequence_length INTEGER"
						])
						
		genome_index = ', '.join(['subject_sequence_id'])
		
		srname = "query_seq_reference"
		
		sequence_records_and_metadata = ', '.join([
										"sequence_name TEXT PRIMARY KEY",
										"query_sequence_id INTEGER",
										"query_sequence_length INTEGER",
										"blast_hit_count INTEGER",
										])
										
		sequence_index = ', '.join([
								"query_sequence_id",
								"blast_hit_count"
								])
										
		
		blname = "blast_results"
										
		blast_results = ', '.join([
						"query_sequence_id INTEGER",
						"subject_sequence_id INTEGER",
						"pident REAL",
						"length INTEGER",
						"mismatch INTEGER",
						"qstart INTEGER",
						"qend INTEGER",
						"sstart INTEGER",
						"send INTEGER",
						"sstrand BOOLEAN",
						"evalue REAL",
						"qcovhsp REAL"
						])
						
		blast_idx = ', '.join([
						"query_sequence_id",
						'subject_sequence_id'
						])	
		
		self.curs.execute("CREATE TABLE IF NOT EXISTS {tblname} ({schema})".format(tblname = grname, schema = genome_records))
		self.curs.execute("CREATE TABLE IF NOT EXISTS {tblname} ({schema})".format(tblname = srname, schema = sequence_records_and_metadata))
		self.curs.execute("CREATE TABLE IF NOT EXISTS {tblname} ({schema})".format(tblname = blname, schema = blast_results))
		
		self.curs.execute("CREATE INDEX IF NOT EXISTS {tblname}_index ON {tblname} ({fields})".format(tblname = grname, fields = genome_index))
		self.curs.execute("CREATE INDEX IF NOT EXISTS {tblname}_index ON {tblname} ({fields})".format(tblname = srname, fields = sequence_index))
		self.curs.execute("CREATE INDEX IF NOT EXISTS {tblname}_index ON {tblname} ({fields})".format(tblname = blname, fields = blast_idx))
		
		self.conn.commit()
		
	#Load as polars df for super fast joins
	def load_genome_index(self):
		pass
		
	def load_sequence_ref_index(self):
		pass
		
		
	def add_genome_records(self, records):
		if self.curs is not None:
			values = ', '.join(['?']*3)
			self.curs.executemany("INSERT OR REPLACE INTO genome_reference VALUES ({values})".format(values=values), records)
			self.conn.commit()
		
	def add_query_records(self, records):
		if self.curs is not None:
			values = ', '.join(['?']*4)
			self.curs.executemany("INSERT OR REPLACE INTO query_seq_reference VALUES ({values})".format(values=values), records)
			self.conn.commit()
		
	def add_blast_records(self, records):
		if self.curs is not None:
			values = ', '.join(['?']*12)
			self.curs.executemany("INSERT OR REPLACE INTO blast_results VALUES ({values})".format(values=values), records)
			self.conn.commit()
		
class blaster:
	def __init__(self,
				seq_obj, 
				min_blast_len, 
				genome_file, 
				MSA_dir,  
				search_type,
				error_files,
				task = 'blastn',
				query_ids = {},
				subject_ids = {}):
				
		self.seq_obj = seq_obj
		self.min_blast_len = min_blast_len
		self.genome_file = genome_file
		self.MSA_dir = MSA_dir
		self.task = task
		self.engine = search_type
		
		self.error_files = error_files
		
		self.seq_name = self.seq_obj.get_seq_name()
		self.seq_type = self.seq_obj.get_old_TE_type()
		self.seq_file = self.seq_obj.get_input_fasta()
		self.seq_length = self.seq_obj.get_length()
		
		self.bed_out_file_dup = None
		self.blast_hits_count = None
		self.blast_out_file = None
		
		#self.file_ids = file_ids
		self.query_ids =  query_ids
		self.subject_ids = subject_ids
	
	def set_params(self):
		#####################################################################################################
		# Code block: Set different elongation number for different elements and do BLAST search
		#####################################################################################################

		try:
			# Since DNA element are significantly shorter than LTR and LINE elements, adjust default parameters
			if "DNA" in self.seq_type:
				self.ex_step = 500
				self.max_extension = 7000
				self.crop_end_gap_win = 100
				self.check_extension_win = 50

			# The average length of SINE elements is around 500 bp, adjust default parameters
			if "SINE" in self.seq_type:
				self.ex_step = 200
				self.max_extension = 1400
				self.min_blast_len = 80
				self.min_blast_len = min(self.min_blast_len, 50)
				self.check_extension_win = 50

			if "helitron" in self.seq_type.lower():
				self.ex_step = 500
				self.max_extension = 7000
				self.crop_end_gap_win = 100
				self.check_extension_win = 50

			if "MITE" in self.seq_type:
				self.ex_step = 100
				self.max_extension = 500
				self.min_blast_len = min(self.min_blast_len, 50)
				self.crop_end_gap_win = 40
				self.check_extension_win = 50

		except Exception as e:
			with open(self.error_files, "a") as f:
				# Get the traceback content as a string
				tb_content = traceback.format_exc()
				f.write(f"Error while running blast for sequence: {self.seq_name}\n")
				f.write(tb_content + '\n\n')
			prcyan(f"Error while running blast for sequence: {self.seq_name}. Main Error: {str(e)}. \n"
				   f"Trace back content: {tb_content}\n")
			
		return

	def blastup(self):
		try:
			# run BLAST search for each FASTA file and return a BED file absolute path
			self.blast_out_file = blast_to_database(self.seq_file, 
													self.genome_file, 
													self.MSA_dir,
													search_type = self.engine,
													task = self.task,
													seq_obj = self.seq_obj)
		
		except Exception as e:
			with open(self.error_files, "a") as f:
				# Get the traceback content as a string
				tb_content = traceback.format_exc()
				f.write(f"Error while running blast for sequence: {self.seq_name}\n")
				f.write(tb_content + '\n\n')
			prcyan(f"Error while running blast for sequence: {self.seq_name}. Main Error: {str(e)}. \n"
				   f"Trace back content: {tb_content}\n")
			
		return
			
	def blast_to_database_format(self):
		blast_records = []
		with open(self.blast_out_file) as fh:
			dat = fh.read()
		
		#outfmt 6 qseqid sseqid pident length mismatch qstart qend sstart send sstrand evalue qcovhsp
		dat = dat.splitlines()
		for line in dat:
			segs = line.split("\t")
			query_name = segs[0]
			qseqid   = self.query_ids[query_name]
			sseqid   = self.subject_ids[segs[1]]
			pident   = float(segs[2])
			length   = int(segs[3])
			mismatch = int(segs[4])
			qstart   = int(segs[5])
			qend     = int(segs[6])
			sstart   = int(segs[7])
			send     = int(segs[8])
			sstrand  = segs[9]
			evalue   = float(segs[10])
			qcovhsp  = float(segs[11])
			
			#Sanity, please
			if sstrand == 'plus':
				sstrand = True
			else:
				sstrand = False
			
			next_blast_record = (qseqid, sseqid, pident, length, mismatch, qstart, qend, sstart, send, sstrand, evalue, qcovhsp,)
			blast_records.append(next_blast_record)
		
		self.blast_hits_count = len(blast_records)
		
		query_record = [(query_name, self.query_ids[query_name], self.seq_length, self.blast_hits_count,)]
		
		#bed_format = self.blast_to_bed(blast_records)
			
		return blast_records, query_record
			
	def blast_to_bed(self, formatted_blast_data):
		'''
		if search_type == "blast":
			bed_cmd = (f"awk 'BEGIN{{OFS=\"\\t\"; counter=0}} !/^#/ {{counter+=1; "
					   f"if ($10~/plus/){{print $2, $8, $9, counter, $3, \"+\", $4, $1}} "
					   f"else {{print $2, $9, $8, counter, $3, \"-\", $4, $1}}}}' < {blast_out_file} > {bed_out_file}")
		'''
		bed_records = []
		counter = 1
		for record in formatted_blast_data:
			if record[9]: #this is strand
				next_bed_record = (record[1],
									record[8],
									record[7],
									counter,
									record[2],
									"+",
									record[3],
									record[0],)
			else:
				next_bed_record = (record[1],
									record[8],
									record[7],
									counter,
									record[2],
									"-",
									record[3],
									record[0],)
				
				
			bed_records.append(next_bed_record)
			
		return bed_records
					
	def run_blast(self):
		self.set_params()
		self.blastup()
		blast_records, sequence_record_update = self.blast_to_database_format()
			
		return blast_records, sequence_record_update
		
	
import multiprocessing

def run_one_blast(args):
	seq_obj, subject_ids, query_ids = args[0], args[1], args[2]
	mn = blaster(seq_obj.seq,
				min_blast_len = seq_obj.min_blast_len, 
				genome_file = seq_obj.genome_file, 
				MSA_dir = seq_obj.MSA_dir,  
				search_type = seq_obj.engine,
				error_files = seq_obj.error_files,
				query_ids = query_ids,
				subject_ids = subject_ids)

	blast_records, seq_records = mn.run_blast()
	
	return blast_records, seq_records

def run_blast_search(sequence_objects, genome_file, database_path, threads = 8):
	ref_genome_index = {}
	te_sequence_index = {}
	
	genome_insert = []
	query_insert = []
	
	gidx = 0
	fa = pyfastx.Fasta(genome_file, build_index = True)
	for seq in fa:
		full = seq.description
		full = full.split()[0]
		l = len(seq)
		ref_genome_index[full] = gidx
		next_record = (full, gidx, l,)
		genome_insert.append(next_record)
		gidx += 1
		
	seqidx = 0
	for item in sequence_objects:
		full = item.seq.name
		path = item.seq.get_input_fasta()
		te_sequence_index[full] = seqidx
		next_query_insert = (full, seqidx, item.seq.get_length(), -1,) #Use -1 to indicate a non-updated value
		query_insert.append(next_query_insert)
		seqidx += 1
		
	#print(te_sequence_index)
	
	#return
	
	args = [(item, ref_genome_index, te_sequence_index,) for item in sequence_objects]
		
	db = database(database_path)
	db.open()
	db.initialize()
	
	db.add_genome_records(genome_insert)
	db.add_query_records(query_insert)
	
	pool = multiprocessing.Pool(threads)
	for result in pool.imap_unordered(run_one_blast, args):
		db.add_blast_records(result[0])
		db.add_query_records(result[1])
	
	db.close()
