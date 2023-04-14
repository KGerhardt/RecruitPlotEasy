from .supports.bam_parser import bam_parser
from .supports.sam_parser import sam_parser
from .supports.blast_parser import blast_parser

from .supports.file_checker import file_checker
from .supports.fasta_reader import read_fasta
from .supports.protein_predictor import pyrodigal_manager

'''
from supports.bam_parser import bam_parser
from supports.sam_parser import sam_parser
from supports.blast_parser import blast_parser

from supports.file_checker import file_checker
from supports.fasta_reader import read_fasta
from supports.protein_predictor import pyrodigal_manager
'''

import sqlite3 as sq
import numpy as np

import argparse
import sys
import os

def convert_array(bytestring):
	return np.frombuffer(bytestring, dtype = np.int32)


class rpe_database_builder:
	def __init__(self, database_path):
		self.db = database_path
		self.name = None
		self.sqlname = None
		self.conn = None
		self.curs = None
		
		self.cur_mag_ids = None
		
	def open(self):
		sq.register_converter("array", convert_array)
		self.conn = sq.connect(self.db)
		self.curs = self.conn.cursor()
		
	def close(self):
		self.curs.close()
		self.curs = None
		self.conn.close()
		self.conn = None
		
	def add_genomes_to_db_from_file(self, file, is_mag = False, predict = False):
		check = file_checker(file)
		check.run_checks()
		
		if check.format != "fasta":
			print("Genomes format didnt't parse.")
			return None
		
		contents, deflines, self.cur_mag_ids = read_fasta(file, is_mag)
				
		self.curs.execute("CREATE TABLE IF NOT EXISTS genome_reference (genome_name TEXT PRIMARY KEY, mag_group TEXT, genome_size INTEGER)")
		to_add = []
		for genome in contents:
			to_add.append(('"'+genome+'"', '"'+self.cur_mag_ids[genome]+'"', len(contents[genome])))

		self.curs.executemany("INSERT OR REPLACE INTO genome_reference VALUES (?, ?, ?)", to_add)
		to_add = None
		
		self.curs.execute("CREATE INDEX IF NOT EXISTS genome_reference_index ON genome_reference (genome_name, mag_group)")

		self.conn.commit()
		
		if predict:
			print("Predicting proteins.")
			self.predict_proteins(contents)
			
		self.cur_mag_ids = None
		
	def predict_proteins(self, contents):
		mn = pyrodigal_manager()
		mn.convert_sequences(contents)
		mn.train_manager()
		mn.predict_genes()
		mn.translate()
		
		to_add = []
		for g in mn.predicted_genes:
			source = mn.source_gen[g]
			mag_id = self.cur_mag_ids[source]
			annot = mn.annot[g]
			start, end = annot[1], annot[2]
			strand = annot[0]
			annot = annot[3]
			
			to_add.append(('"'+g+'"', '"'+mag_id+'"', '"'+source+'"', strand, start, end, annot,))
		
		self.curs.execute("CREATE TABLE IF NOT EXISTS protein_reference (gene_name TEXT PRIMARY KEY, mag_group TEXT, genome_name TEXT, strand INTEGER, protein_start INTEGER, protein_end INTEGER, annotation TEXT)")
		self.curs.executemany("INSERT OR REPLACE INTO protein_reference VALUES (?, ?, ?, ?, ?, ?, ?)", to_add)
		self.curs.execute("CREATE INDEX IF NOT EXISTS protein_reference_index ON protein_reference (gene_name, mag_group)")
		self.conn.commit()
			
	def add_genomes_to_db_from_header(self, header):
		self.cur_mag_ids = None
		to_add = []
		for line in header:
			if line.startswith("@SQ"):
				segs = line.strip().split("\t")
				genome = segs[1][3:]
				seqlen = int(segs[2][3:])
				this_row = ('"'+genome+'"', '"'+genome+'"', seqlen,)
				to_add.append(this_row)
	
		self.curs.execute("CREATE TABLE IF NOT EXISTS genome_reference (genome_name TEXT PRIMARY KEY, mag_group TEXT, genome_size INTEGER)")
		#Do NOT overwrite an existing record.
		self.curs.executemany("INSERT OR IGNORE INTO genome_reference VALUES (?, ?, ?)", to_add)
		
		to_add = None
		
		self.curs.execute("CREATE INDEX IF NOT EXISTS genome_reference_index ON genome_reference (genome_name, mag_group)")

		self.conn.commit()

	def add_genomes_from_reads(self, gens):
		self.curs.execute("CREATE TABLE IF NOT EXISTS genome_reference (genome_name TEXT PRIMARY KEY, mag_group TEXT, genome_size INTEGER)")
		
		to_add = []
		for g in gens:
			row = (g, g, int(gens[g]),)
			to_add.append(row)
			
		self.curs.executemany("INSERT OR IGNORE INTO genome_reference VALUES (?, ?, ?)", to_add)

		self.curs.execute("CREATE INDEX IF NOT EXISTS genome_reference_index ON genome_reference (genome_name, mag_group)")

		self.conn.commit()
		
	def add_reads_to_db(self, file):
		check = file_checker(file)
		check.run_checks()
		
		table_name = check.sqlname
		
		if check.format == "sam":
			parser = sam_parser(file)
			self.add_genomes_to_db_from_header(parser.header)
		if check.format == "bam":
			parser = bam_parser(file)
			self.add_genomes_to_db_from_header(parser.header)
		if check.format == "blast":
			parser = blast_parser(file)
						
		self.curs.execute("DROP TABLE IF EXISTS {tab}".format(tab=table_name))
		self.curs.execute("DROP TABLE IF EXISTS {tab}_query_reference".format(tab=table_name))
		self.curs.execute("DROP TABLE IF EXISTS {tab}_target_reference".format(tab=table_name))
		self.curs.execute("CREATE TABLE {tab}_query_reference (query TEXT, query_id INTEGER)".format(tab=table_name))
		self.curs.execute("CREATE TABLE {tab}_target_reference (target TEXT, target_id INTEGER)".format(tab=table_name))
		self.curs.execute("CREATE TABLE {tab} (qid INTEGER, tid INTEGER, local REAL, global REAL, pct_aln REAL, start INTEGER, stop INTEGER)".format(tab=table_name))
		
		self.conn.commit()
		
		to_add = []
		seen_genomes = {}
		
		query_ref = {}
		qr_num = 0
		target_ref = {}
		tr_num = 0
		
		print("Adding reads.")
		reads_added = 0
				
		for grouping in parser:
			if grouping is None:
				continue
				
			next_set = []
			#One result looks like so
			#[query, target, pct_id_local, pct_id_global, pct_alignment, covered_ranges]
			
			query_gen = '"'+grouping[0]+'"'
			tgt_gen = '"'+grouping[1]+'"'
			
			#grouping[5] is a list of tuples with starts, stops
			max_aln = grouping[5][len(grouping[5])-1][1]
			
			if tgt_gen not in target_ref:
				target_ref[tgt_gen] = tr_num
				tr_num += 1
			if query_gen not in query_ref:
				query_ref[query_gen] = qr_num
				qr_num += 1
				
			curq, curt = query_ref[query_gen], target_ref[tgt_gen]

			if tgt_gen not in seen_genomes:
				seen_genomes[tgt_gen] = max_aln
			else:
				if max_aln > seen_genomes[tgt_gen]:
					seen_genomes[tgt_gen] = max_aln
			
			local, glob, pctaln = grouping[2], grouping[3], grouping[4]
			
			for tuple in grouping[5]:
				start, stop = tuple[0], tuple[1]
				next_set = (curq, curt, local, glob, pctaln, start, stop)
				to_add.append(next_set)
				
			
			if len(to_add) > 10**6:
				reads_added += len(to_add)
				self.curs.executemany("INSERT INTO {tab} VALUES (?, ?, ?, ?, ?, ?, ?)".format(tab=table_name), to_add)
				self.conn.commit()
				to_add = []
				print(reads_added, "reads added so far.")
		
		if len(to_add) > 0:
			reads_added += len(to_add)
			self.curs.executemany("INSERT INTO {tab} VALUES (?, ?, ?, ?, ?, ?, ?)".format(tab=table_name), to_add)
			self.conn.commit()
			to_add = []
			
		print(reads_added, "reads added in total.")

		to_add = []
		for tgt in target_ref:
			to_add.append((tgt, target_ref[tgt],))
		self.curs.executemany("INSERT OR REPLACE INTO {tab}_target_reference VALUES (?, ?)".format(tab = table_name), to_add)

		to_add = []
		for q in query_ref:
			to_add.append((q, query_ref[q],))
		self.curs.executemany("INSERT OR REPLACE INTO {tab}_query_reference VALUES (?, ?)".format(tab = table_name), to_add)
		
		self.conn.commit()
		to_add = None	
		
		self.curs.execute("CREATE INDEX IF NOT EXISTS {tab}_index ON {tab} (qid, tid)".format(tab=table_name))
		self.curs.execute("CREATE INDEX IF NOT EXISTS {tab}_query_reference_index ON {tab}_query_reference (query)".format(tab=table_name))
		self.curs.execute("CREATE INDEX IF NOT EXISTS {tab}_target_reference_index ON {tab}_target_reference (target)".format(tab=table_name))
		self.conn.commit()
		
		if check.format == "blast":
			#We've already added sequences as headers from sam/bam
			self.add_genomes_from_reads(seen_genomes)
		
	
def build_opts():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
			description='''
	''')
	parser.add_argument('-g', '--genome',  dest = 'genome', default = None, 
	help =  'A path to a genome in nt FASTA format.')
	
	parser.add_argument('-gf', '--genome_file',  dest = 'gf', default = None, 
	help =  'A file containing a list of FASTA-format genomes, one per line. Do not use with --genome at the same time.')
	
	parser.add_argument('--mag',  dest = 'is_mag', action = 'store_true', 
	help =  'Any genomes to add are MAGs: all sequences in each genome file will be treated as the contigs of one genome.')
	
	parser.add_argument('--predict',  dest = 'predict', action = 'store_true', 
	help =  'Predict proteins for each genome file.')
	
	parser.add_argument('-r', '--reads',  dest = 'reads', default = None, 
	help =  'A path to a set of aligned nt reads in SAM, BAM, or tabular BLAST format.')
	
	parser.add_argument('-rf', '--reads_file',  dest = 'rf', default = None, 
	help =  'A file containing a list of aligned nt reads in SAM, BAM, or tabular BLAST format, one set of reads per line. Do not use with --reads at the same time.')
	
	parser.add_argument('-d', '--database',  dest = 'db', default = None, 
	help =  'Path to the database to create. Required.')
	
	args, unknown = parser.parse_known_args()
	
	return parser, args

def load_file_paths(list_file):
	file_paths = []
	with open(list_file) as fh:
		for line in fh:
			cleaned = line.strip()
			if os.path.exists(cleaned):
				file_paths.append(cleaned)
			else:
				print("Can't find file:", cleaned, "skipping this file.")
	return file_paths
	
def run_build():
	parser, opts = build_opts()
	if len(sys.argv) < 3:
		parser.print_help()
		sys.exit()
	
	genome = opts.genome
	is_mag = opts.is_mag
	pred = opts.predict
	reads = opts.reads
	db = opts.db
	
	gf = opts.gf
	rf = opts.rf
	
	if gf is not None:
		gf_paths = load_file_paths(gf)
		if len(gf_paths) == 0:
			gf_paths = None
	else:
		gf_paths = None
	if rf is not None:
		rf_paths = load_file_paths(rf)
		if len(rf_paths) == 0:
			rf_paths = None
	else:
		rf_paths = None
	
	if reads is not None and rf_paths is not None:
		sys.exit("Use either --reads or --reads_file, not both.")
		
	if genome is not None and gf_paths is not None:
		sys.exit("Use either --genomes or --genome_file, not both.")
	
	if db is None:
		sys.exit("You need to specify a database.")
		
	if reads is None and genome is None and gf_paths is None and rf_paths is None:
		sys.exit("You need to supply at least one of --reads, --genome, --reads_file, or --genome_file")
	
	
	mn = rpe_database_builder(db)
	mn.open()
	
	if reads is not None:
		mn.add_reads_to_db(reads)
	
	if rf_paths is not None:
		for read_file in rf_paths:
			mn.add_reads_to_db(read_file)

	if genome is not None:
		mn.add_genomes_to_db_from_file(file = genome, 
										is_mag = is_mag, 
										predict = pred)	
										
	if gf_paths is not None:
		for genome_file in gf_paths:
			mn.add_genomes_to_db_from_file(file = genome_file, 
										is_mag = is_mag, 
										predict = pred)
			
	mn.close()
	
