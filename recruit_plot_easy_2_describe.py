import plotly

import plotly.graph_objects as go

import sys
import os
import numpy as np

import sqlite3 as sq

import argparse

#Modified from the plot codebase
class rpdb_descriptor:
	def __init__(self, db = None, mf = None, sf = None, hm = None):
		self.db = db
		self.conn = None
		self.curs = None
		
		self.db_schema = None
		self.mag_list = None
		self.mag_dict = None
		self.genome_sizes = None
		self.samples = None
		
		self.samples_info = None
		
		self.mags_and_genomes_file = mf
		self.per_samples_file = sf
		self.presence_plot_file = hm
		
	def open(self):
		self.conn = sq.connect(self.db)
		self.curs = self.conn.cursor()
		
	def close(self):
		self.curs.close()
		self.conn.close()
		self.conn = None
		self.curs = None
	
	def check_db(self):
		sql = "SELECT * FROM sqlite_master"
		self.db_schema = self.curs.execute(sql).fetchall()
	
	def collect_mags(self):
		self.mag_list = []
		sql = "SELECT * FROM genome_reference"
		self.mag_dict = {}
		self.genome_sizes = {}
		for result in self.curs.execute(sql).fetchall():
			genome, mag, size = result[0], result[1], result[2]
			mag = mag[1:-1] #remove enclosing quotes
			genome = genome[1:-1] #remove enclosing quotes
			self.mag_dict[genome] = mag
			self.genome_sizes[genome] = size
			self.mag_list.append(mag)
		
		self.mag_list = set(self.mag_list)
		
	def collect_samples(self):
		self.samples = []
		for t in self.db_schema:
			type = t[0]
			if type == "table":
				name = t[1]
				if name.endswith("_query_reference"):
					sample_name = t[1].split("_query_reference")[0]
					self.samples.append(sample_name)			
		
	def collect_mags_per_sample(self, sample):
		#collect target idss
		sql = "SELECT * FROM '{sample}_target_reference'".format(sample=sample)
		sample_dict = {}
		for result in self.curs.execute(sql).fetchall():
			genome, id = result[0][1:-1], result[1]
			sample_dict[id] = genome
		
		#this_sample = []
		sql = 'SELECT tid, COUNT(qid), SUM(stop-start) FROM "{sample}" GROUP BY tid'.format(sample=sample)
		for result in self.curs.execute(sql).fetchall():
			genome_id, read_count, bp_count = result[0], result[1], result[2]
			genome_name = sample_dict[genome_id]
			mag = self.mag_dict[genome_name]
			next_row = (sample, mag, genome_name, read_count, bp_count,)
			
			self.samples_info.append(next_row)
		
	def write_mags_file(self):
		if self.mags_and_genomes_file is not None:
			fh = open(self.mags_and_genomes_file, "w")
			print("MAG_group", sep = "\t", file = fh)
			for mag in self.mag_list:
				print(mag, sep = "\t", file = fh)
			fh.close()
		
	def write_samples_file(self):
		if self.per_samples_file is not None:
			fh = open(self.per_samples_file, "w")
			print("sample", "MAG_group", "contig_name", "num_reads", "num_bp", file = fh)
			for row in self.samples_info:
				print(*row, sep = "\t", file = fh)
			fh.close()
			
	def prepare_plot(self):
		if self.presence_plot_file is not None:
	
			mag_enum = {} #get mag to rows
			mag_list = [] #ensure name ordering
			idx = 0
			for mag in sorted(self.mag_dict):
				mag_list.append(mag)
				mag_enum[mag] = idx
				idx += 1
			
			sample_enum = {} #get to columns
			samp_list = [] #ensure name ordering
			idx = 0
			for sample in sorted(self.samples):
				samp_list.append(sample)
				sample_enum[sample] = idx
				idx += 1
				
			reads_per_samp_by_mag = np.empty(shape = (len(mag_enum), len(sample_enum)), dtype = np.float_)
			reads_per_samp_by_mag[:] = np.nan
					
			for row in self.samples_info:
				sample = row[0]
				mag = row[1]
				num_reads = row[3]
				
				mag_idx = mag_enum[mag]
				samp_idx = sample_enum[sample]
				
				reads_per_samp_by_mag[mag_idx, samp_idx] = num_reads
			
			hover_format = 	"Sample Name: %{x:c}<br>" +\
							"MAG Name: %{y:c}<br>" +\
							"Number of Reads: %{z:d}<br>" +\
							"<extra></extra>"
			
			fig = go.Figure(
					go.Heatmap(x = samp_list,
								y = mag_list,
								z = reads_per_samp_by_mag,
								hovertemplate = hover_format)
							)
			
			html_config = {
				'scrollZoom': True,
					'toImageButtonOptions': {
					'format': "svg", # one of png, svg, jpeg, webp
					'height': 1080,
					'width': 1920,
					'scale': 1 # Multiply title/legend/axis/canvas sizes by this factor
			  }
			}
				
			fig.write_html(self.presence_plot_file, config = html_config)
		
	def run(self):
		self.open()
		self.check_db()
		self.collect_mags()
		self.collect_samples()
		
		#Adds the info to samples_info
		self.samples_info = []
		for sample in self.samples:
			self.collect_mags_per_sample(sample)
			
		self.close()
			
		self.write_mags_file()
		self.write_samples_file()
		self.prepare_plot()
		
		
#Should just be 3 file inputs + database
def options():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
			description='''Module for assessing a RecruitPlotEasy database and selecting MAGs to plot.
			
			* To select only some MAGs to plot for the RecruitPlotEasy plotting module, use the --mags_file argument.
			
			* To get a text file showing what MAGs are in which samples and how many reads they have, use --samples_file
			
			* To get a heatmap of MAG presence/absence with read counts in each sample, use --heatmap
			
			This module produces only the outputs you supply. If you omit an output argument, nothing will be produced for that argument.''')
			
	parser.add_argument('-d', '--database',  dest = 'db', default = None, 
	help =  'Path to the RecruitPlotEasy database you want to make plots from. Required.')
	
	parser.add_argument('-m', '--mags_file',  dest = 'mags_file', default = None,
	help =  "An output file containing a list the names of MAGs in the database.\nRemove rows from this file before giving it to RecruitPlotEasy's plot module with the --mag_file argument to plot only the remaining MAGs")
	
	parser.add_argument('-s', '--samples_file',  dest = 'samps_file', default = None,
	help =  "An output file containing the MAGs and contigs in each sample with a list of read and base pair counts.")

	parser.add_argument('--heatmap',  dest = 'heatmap', default = None,
	help =  "An output HTML interactive heatmap showing the presence of each MAG per sample and their read counts. You should probably end this file path with '.html'")
	
	args, unknown = parser.parse_known_args()
	
	return parser, args

def run_descriptor():
	parser, opts = options()
	db = opts.db
	if db is None:
		print("You must supply a database. Quitting.")
		parser.print_help()
		sys.exit()
	else:
		if not os.path.exists(db):
			print("Database", db, "could not be found.")
			print("Have you created one with RecruitPlotEasy build yet? Quitting.")
			sys.exit()
	
	mf = opts.mags_file
	sf = opts.samps_file
	hm = opts.heatmap
	
	if mf is None and sf is None and hm is None:
		print("You must supply at least one output file. Quitting.")
		parser.print_help()
		sys.exit()
		
	mn = rpdb_descriptor(db = db, mf = mf, sf = sf, hm = hm)
	mn.run()
		
#run_descriptor()