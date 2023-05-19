import plotly

import plotly.graph_objects as go
from plotly.subplots import make_subplots

import sys
import os
import numpy as np

import sqlite3 as sq

import argparse

#import multiprocessing

class rpdb:
	def __init__(self, db, id_cut = 95, id_step = 0.5, gen_step = 1000,
				criteria = "local", id_measure = "local", do_prot = False,
				output_base = "recruitment_plots"):
		
		self.db = db
		self.conn = None
		self.curs = None
		
		self.samples = None
		self.genomes = None
		self.proteins = None
		self.do_prot = do_prot
		
		self.current_sample = None
		
		self.mags_in_sample = None
		
		self.current_mag = None
		self.contigs_in_mag = None
		self.mag_contig_ids = None
		self.contig_lens_by_id = None
		
		self.best_hit_criteria = criteria #or "global or "aln"
		self.most_recent_query = None
		
		self.pct_id_metric = id_measure # local or global
		
		self.do_genes = False
		self.raw_data = None
		self.cur_gen_size = None
		self.x_binstarts = None
		self.x_binends = None
		self.y_bins = None
		
		#Hard coded on purpose.
		self.y_min = 70
		self.y_max = 100
		
		self.id_cutoff = id_cut
		
		self.y_step = id_step
		self.x_step = gen_step
		
		self.outbase = output_base
		self.recplot = None

	def open(self):
		self.conn = sq.connect(self.db)
		self.curs = self.conn.cursor()
		
	def close(self):
		self.curs.close()
		self.conn.close()
		self.conn = None
		self.curs = None
		
	def parse_db(self):
		self.samples = []
		self.genomes = None
		self.proteins = None
		#Tables
		tables = self.curs.execute("SELECT name FROM sqlite_master").fetchall()
		for t in tables:
			t = t[0]
			if "index" not in t:
				if t == "genome_reference":
					self.genomes = {}
					for result in self.curs.execute("SELECT * FROM genome_reference").fetchall():
						genome = result[0][1:-1]
						mag = result[1][1:-1]
						size = int(result[2])
						if mag not in self.genomes:
							self.genomes[mag] = []
						self.genomes[mag].append((genome, size))
						
				elif t == "protein_reference":
					self.proteins = {}
					#protein_reference (gene_name, mag_group, genome_name, strand INTEGER, protein_start INTEGER, protein_end INTEGER, annotation TEXT)")

					for result in self.curs.execute("SELECT * FROM protein_reference").fetchall():
						gene = result[0][1:-1]
						
						mag = result[1][1:-1]
						
						genome = result[2][1:-1]
						strand = int(result[3])
						start = int(result[4])
						end = int(result[5])
						annot = result[6]
						if mag not in self.proteins:
							self.proteins[mag] = {}
						if genome not in self.proteins[mag]:
							self.proteins[mag][genome] = []
						self.proteins[mag][genome].append((gene, strand, start, end, annot))
				else:
					if t.endswith("_query_reference"):
						continue
					elif t.endswith("_target_reference"):
						continue
					else:
						self.samples.append(t)
	
	def set_sample(self, sample):
		self.current_sample = sample
		self.mags_in_sample = []
		#load mags in this sample
		for r in self.curs.execute("SELECT mag_group FROM genome_reference WHERE genome_name IN (SELECT target FROM {sample}_target_reference)".format(sample=sample)).fetchall():
			self.mags_in_sample.append(r[0])
			
		self.mags_in_sample = list(set(self.mags_in_sample))
			
	def set_mag(self, mag):
		self.current_mag = mag[1:-1]
	
	def craft_query(self):
		mag_to_contigs = "SELECT genome_name FROM genome_reference WHERE mag_group='\"{mag}\"'"
		mag_to_contigs = mag_to_contigs.format(mag = self.current_mag)
		contigs_to_ids = "SELECT target_id FROM {sample}_target_reference WHERE {sample}_target_reference.target IN ({contigs})"
		contigs_to_ids = contigs_to_ids.format(sample = self.current_sample, contigs = mag_to_contigs)
		
		self.mag_contig_ids = {}
		#print(contigs_to_ids.replace("SELECT target_id FROM", "SELECT target, target_id FROM"))
		for tid in self.curs.execute(contigs_to_ids.replace("SELECT target_id FROM", "SELECT target, target_id FROM")).fetchall():
			target, id = tid[0], tid[1]
			self.mag_contig_ids[target[1:-1]] = id
			
		ids_to_selection = "SELECT * FROM {sample} WHERE {sample}.{lg} >= {lgc} AND {sample}.tid IN ({ids})"
		ids_to_selection = ids_to_selection.format(sample = self.current_sample, ids = contigs_to_ids, lg = self.pct_id_metric, lgc = self.y_min)
		
		self.most_recent_query = ids_to_selection
		
	def convert_to_bins(self, starts, ends, target):
		#Starts, ends are lists of ints from the read that correspond to the read's [start, end) pos in the genome.
		#bin_starts, bin_ends are genome coordinates indicating cutoffs where bases are to be counted. Rules:
		#Every start has a matching end. 
		#Overlaps are resolved by eliminating any overlap that is a subsection and by choosing the midpoint of the overlap starts/stops otw.
		#This function expects that to have already been done. i.e., each base can and must fall into EXACTLY one bin.
		#Gotta figure out hbins, counts now.
		#starts, ends = [], []
		
		returns_by_bin = {}
		
		bin_starts, bin_ends = self.x_binstarts[target], self.x_binends[target]
			
		for s, e in zip(starts, ends):
			start_bin = np.searchsorted(bin_ends, s, side = 'right')
			end_bin = np.searchsorted(bin_ends, e, side = 'right')
			#print("----------------")
			#print(s, e, start_bin, end_bin, bin_starts[start_bin], bin_ends[end_bin])
			
			if start_bin == end_bin:
				#start and end bin are the same; the whole read falls into this bin.
				if start_bin in returns_by_bin:
					returns_by_bin[start_bin] += (e-s)+1
				else:
					returns_by_bin[start_bin] = (e-s)+1
			else:
				#Figure out how many bases need adding in total as a running tracker of how many are left.
				total_bases = e-s + 1
				#How far can we go before we run out of bin
				current_end = bin_ends[start_bin]
				#How many do we add to this bin, then?
				to_add = current_end - s
				#remove those from the tracker
				total_bases -= to_add
				#Final position is not counted in the same bin
				total_bases += 1
				#Add as needed.
				if start_bin in returns_by_bin:
					returns_by_bin[start_bin] += to_add -1 
				else:
					returns_by_bin[start_bin] = to_add - 1
				#Update start to start of the next bin	
				s = current_end
					
				#Repeat until done.
				while total_bases > 0:
					#Move to the next bin
					start_bin += 1
					#Next end
					if start_bin == len(bin_ends) or end_bin == len(bin_ends):
						break
					
					current_end = bin_ends[start_bin]
					
					if current_end > e:
						#This is the last bin for this s-e window, so we use the remaining bases as the count instead of the bases to the end of the bin
						if start_bin in returns_by_bin:
							returns_by_bin[start_bin] += total_bases
						else:
							returns_by_bin[start_bin] = total_bases
						#Doesn't matter what this is, as long as it's less than 0
						total_bases = -1
					else:
						#This is not the last bin to add to, so we use distance to the bin as the marker.
						to_add = current_end - s
						total_bases -= to_add
						#We don't count the final position in the bin, as it's non-inclusive
						total_bases += 1
						if start_bin in returns_by_bin:
							returns_by_bin[start_bin] += to_add - 1
						else:
							returns_by_bin[start_bin] = to_add - 1
						
						s = current_end
						
		return returns_by_bin
	
	def load_sample(self):
		if not os.path.exists(self.outbase):
			os.mkdir(self.outbase)
		this_samp = os.path.normpath(self.outbase+"/"+self.current_sample)
		if not os.path.exists(this_samp):
			os.mkdir(this_samp)
	
		num_ybins = int((self.y_max-self.y_min) / self.y_step)+1
		
		self.y_bins = np.linspace(self.y_min, self.y_max, num = num_ybins, dtype = np.float_)
		
		self.contig_lens_by_id = {}
		
		self.raw_data = {}
		
		#actually, we should just move bin finding to the recplot class.
		for tuple in self.genomes[self.current_mag]:
			contig = tuple[0]
			contig_len = tuple[1]
			
			if contig not in self.mag_contig_ids:
				continue
			
			contig_id = self.mag_contig_ids[contig]

			self.contig_lens_by_id[contig_id] = contig_len
				
			self.raw_data[contig_id] = {}
		
		print("Loading", self.current_mag)
		
		#execute read pull
		self.curs.execute(self.most_recent_query)
		
		#1 million rows at a time max
		chunk_size = 10**6
		
		next_group = self.curs.fetchmany(chunk_size)
		keep_going = len(next_group) > 0
		while keep_going:
			for read in next_group:
				#target = read[1][1:-1]
				target = read[1]
				local = read[2]
				glob = read[3]
				aln = read[4]
				
				start = read[5]
				end = read[6]
				
				#Maybe this should be pushed to the query step...
				#Skip insuff. pct. ID reads
				if self.pct_id_metric == "local":
					current_y_bin = (local - self.y_min) // self.y_step
				else:
					current_y_bin = (glob - self.y_min) // self.y_step
					
				current_y_bin = int(current_y_bin)
				
				#Only add the ybins we need to add.
				if current_y_bin not in self.raw_data[target]:
					this_len = self.contig_lens_by_id[target]
					self.raw_data[target][current_y_bin] = np.zeros(this_len, dtype = np.int32)

				self.raw_data[target][current_y_bin][start:end] += 1
				

			next_group = self.curs.fetchmany(chunk_size)
			keep_going = len(next_group) > 0
		
		if self.proteins is not None:
			if self.current_mag in self.proteins:
				protein_subset = self.proteins[self.current_mag]
			else:
				protein_subset = None
		if not self.do_prot:
			protein_subset = None
		
		self.recplot = recplot(data = self.raw_data, 
							#x1 = self.x_binstarts, 
							#x2 = self.x_binends, 
							y = self.y_bins,
							#cut = self.id_cutoff,
							id_step = self.y_step,
							genome_step = self.x_step,
							contig_sizes = self.contig_lens_by_id,
							protein_info = protein_subset,
							contig_name_dict = self.mag_contig_ids,
							mag_name = self.current_mag,
							sample = self.current_sample,
							outdir = self.outbase) 
		
		self.recplot.build()
	
class recplot:
	def __init__(self, data, y, 
				#cut, 
				id_step, genome_step, contig_sizes, protein_info, 
				contig_name_dict, tad = 80, mag_name = None, sample = None,
				outdir = "recruitment_plots", selected_mags = None,
				font_size = 18, in_group_col = '#ff7f0e', 
				out_group_col = '#1f77b4', overlay_alpha = 0.25,
				scroll_zoom = False, save_static_img_as = "svg"):
		
		self.raw_data = data
		self.contig_sizes = contig_sizes
		self.gen_step = genome_step
		#Flip k-v
		self.ct_dict = dict([value, key] for key, value in contig_name_dict.items())
		self.ct_names = None
		
		self.tad_level = tad
		self.tad_values = None
		self.breadths = None
		self.contig_ends = None
		self.tad_ticks = None
		self.tad_names = None
		self.tad_contigs = None
		self.tad_max = 0
		
		self.prot = protein_info
		self.do_prot = (protein_info is not None)
		self.protein_labels = None

		self.bin_left = None
		self.bin_right = None
		
		self.bin_mids = None
	
		self.data = None

		#sizes are implicitly given by x1 and x2
		#self.contig_sizes = sizes
		
		self.base_pd = None
		
		self.max_count = 0
		
		self.y = y
		self.id_step = id_step
		self.cutoff = None
		self.selected = None
		
		self.upper_left_in = None
		self.upper_left_out = None
		self.tad_in = None
		self.cov_in = None
		
		self.depth_hist_breaks = None
		self.upper_right_in = None
		self.upper_right_out = None
		#In-grp histogram local maxima here
		self.peaks = None
		
		self.lower_right_data = None
		
		self.output_base = outdir

		self.sample = sample
		self.mag = mag_name
		
		if self.do_prot:
			self.plot_name = os.path.normpath(self.output_base + "/" +self.sample + "/" + mag_name + "_proteins_recruitment_plot.html")
		else:
			self.plot_name = os.path.normpath(self.output_base + "/" +self.sample + "/" + mag_name + "_recruitment_plot.html")
		

		self.save_img_fmt = save_static_img_as
		self.scrollable = scroll_zoom
		
		self.html_config = {
			'scrollZoom': self.scrollable,
				'toImageButtonOptions': {
				'format': self.save_img_fmt, # one of png, svg, jpeg, webp
				'filename': self.plot_name,
				'height': 1080,
				'width': 1920,
				'scale': 1 # Multiply title/legend/axis/canvas sizes by this factor
		  }
		}
		
		self.bot_right_line_col = 'rgb(66,146,198)'
		self.main_plot_fill_colorscale = ['rgb(247,251,255)', 'rgb(222,235,247)', 'rgb(198,219,239)', 'rgb(158,202,225)', 'rgb(107,174,214)', 'rgb(66,146,198)', 'rgb(33,113,181)', 'rgb(8,81,156)', 'rgb(8,48,107)']
		
		self.main_plot_highlight_col = in_group_col
		
		self.main_plot_highlight_alpha = overlay_alpha
		
		#'#1f77b4' is a medium blue
		#'#ff7f0e' is a burnt orange
		self.in_group_depth_col = in_group_col
		self.out_group_depth_col = out_group_col
		
		self.axis_font_size = font_size
		
	#Data comes in as a per-base count 
	def bin_raw(self):
		print("Processing", self.mag)
		self.data = {}
		self.bin_left = {}
		self.bin_right = {}
		
		#self.tad_values = {}
		self.breadths = {}
		self.contig_ends = []
		previous_end = 0
		next_end = 0
		
		#tad_steps = 20

		total_tads = 11
		#These are the cut sizes.
		tads = np.arange(0.05, 0.50, 0.05)
		
		tad_labels = ["TAD-"+str(int(v)) for v in np.linspace(90, 10, num = 9)]
		raw = ["Average Depth"]
		raw.extend(tad_labels)
		raw.append("Median Depth")
		tad_labels = raw
		self.tad_names = raw
		self.tad_contigs = []
		
		#print(self.tad_names)
		
		tad_data = {}
		num_contigs = len(self.raw_data)
		ct_count = 0
		for ybin in range(len(self.y)-1, -1, -1):
			tad_data[ybin] = []
			
		if self.do_prot:
			self.protein_labels = []
		
		for contig in self.raw_data:
			#print("Contig", contig)
			ct_name = self.ct_dict[contig]
			self.tad_contigs.append(ct_name)
			
			contig_len = self.contig_sizes[contig]
			
			next_end = previous_end+contig_len
			prev_record = previous_end
			
			next_indices = np.linspace(previous_end, next_end, num = total_tads, dtype = np.int32).tolist()
			
			self.contig_ends.extend(next_indices)
			
			previous_end = next_end+1
			
			self.contig_ends.append(previous_end)
			
			#Okay, so we need a dict of per-y-bin tads from high to low arranged in a tad row, contig col arrangement 
			#Then we heatmap that below.
			
			if self.do_prot:
				if ct_name not in self.prot:
					print("No proteins detected for", ct_name)
					print("This will be omitted!")
					continue
					
				this_genome = self.prot[ct_name]
				
				#Bin left, bin right
				bin_breaks = [0]
				first_start = this_genome[0][2]
				
				loc = 0
				for protein_tuple in this_genome:
					last_end = bin_breaks[loc]
					start, end = protein_tuple[2], protein_tuple[3]
					
					self.protein_labels.append("Contig: " + ct_name + "<br>" +
											   "Range: " + str(last_end+1)+ "-" +str(start)+"<br>"+
											   "Intergenic")
					
					
					gene = protein_tuple[0]
					strand = str(protein_tuple[1])
					annotation = "<br>".join(protein_tuple[4].split(";"))
					
					self.protein_labels.append("Contig: " + ct_name + "<br>" +
											   "Gene: " + gene + " <br>"+
											   "Range: " + str(start) + "-" + str(end)+"<br>"+
											   "Strand: " + strand +"<br>"+
											   "Annot: " + annotation)
					
					#This is a panic case for when genes overlap.
					#Split the overlap down the middle.
					if start <= last_end:
						move = int((last_end + start)/2)
						bin_breaks[loc] = move
						start = move + 1
						
					bin_breaks.append(start)
					bin_breaks.append(end)
					
					loc += 2
				
				#Cap.
				last_end = bin_breaks[loc]
				self.protein_labels.append("Contig: " + ct_name + "<br>" +
										   "Range: " + str(last_end+1)+ "-" +str(contig_len)+"<br>"+
										   "Intergenic")
				
				bin_breaks.append(contig_len)
				bin_breaks = np.array(bin_breaks, dtype = np.int32)
				bin_count = len(bin_breaks)
				
				self.data[contig] = np.zeros(shape = (self.y.shape[0], bin_count-1), dtype = np.int32)
				
			else:
				#math-less ceil function
				bin_count = -(-contig_len//self.gen_step) + 1
				
				self.data[contig] = np.zeros(shape = (self.y.shape[0], bin_count-1), dtype = np.int32)
				
				#This works within a single contig, but resets X axis to zero when adding bins
				bin_breaks = np.linspace(0, contig_len, dtype = np.int32, num = bin_count)
				
			
			self.bin_left[contig] = bin_breaks[:-1]
			self.bin_right[contig] = bin_breaks[1:]
				
			breadth_and_depth = np.zeros(contig_len, dtype = np.int32)
			
			self.breadths[contig] = np.zeros(self.y.shape[0], dtype = np.float_)
			#self.tad_values[contig] = np.zeros(shape = (self.y.shape[0], total_tads), dtype = np.float_)
			
			last_added_breadth = 0
			last_added_depth = np.zeros(total_tads, dtype = np.float_)
		
			starts = np.multiply(tads, contig_len)
			ends = np.subtract(contig_len, starts)
				
			starts = starts.astype(int).tolist()
			ends = ends.astype(int).tolist()
				
			is_even = (contig_len % 2 == 0)
			if is_even:
				median_idx = contig_len/2 - 1
			else:
				median_idx = ((contig_len - 1)/2) -1
			
			median_idx = int(median_idx)
			
			#Descending, considers all possible y values for main plot, not just observed.
			for ybin in range(len(self.y)-1, -1, -1):
				pct = self.y[ybin]
				#We have data for this y
				if ybin in self.raw_data[contig]:
					breadth_and_depth += self.raw_data[contig][ybin]
				
					#calc here
					last_added_breadth = np.count_nonzero(breadth_and_depth) / contig_len
					
					depths = np.sort(breadth_and_depth)
					
					next_tads = []
					
					#Non-truncated average, or TAD-100
					next_tads.append(np.mean(depths))
					#TADs for  90 to 10, by 10 pct at a time.
					for s, e in zip(starts, ends):
						next_tads.append(np.mean(depths[s:e]))
					
					#Median depth, or TAD-50
					if is_even:
						next_tads.append((depths[(median_idx-1)] + depths[(median_idx+1)])/2)
					else:
						next_tads.append(depths[median_idx])
					
					last_added_depth = next_tads
											
					depths = None
					
					self.breadths[contig][ybin] = last_added_breadth
					#self.tad_values[contig][ybin] = last_added_depth
					
					tad_data[ybin].append(last_added_depth)
					
					#Histogram handles binning no matter where the edges are.
					next_row = np.histogram(np.arange(0, contig_len, dtype = np.int32), 
											bins = bin_breaks,
											weights = self.raw_data[contig][ybin],)
					#self.raw_data[contig][ybin].shape					
					
					self.data[contig][ybin, :] += next_row[0]
				else:
					#We don't have data for this y
					self.breadths[contig][ybin] = last_added_breadth
					#self.tad_values[contig][ybin] = last_added_depth
					tad_data[ybin].append(last_added_depth)
		
			self.bin_left[contig] = np.add(self.bin_left[contig], prev_record)
			self.bin_right[contig] = np.add(self.bin_right[contig], prev_record)
		#The final contig end is illusory and needs trimmed.
		self.contig_ends = self.contig_ends[:len(self.contig_ends)-1]
		
		self.tad_values = tad_data
		for ybin in self.tad_values:
			combined = np.vstack(self.tad_values[ybin])
			self.tad_values[ybin] = np.transpose(combined)
			self.tad_max = max([self.tad_max, combined.max()])
			#print(self.tad_values[ybin].shape)
		
		#print(self.tad_values.keys())
		
	#Join up the windows from the incoming data into one big df.
	#each window corresp. to a single contig
	def concatenate(self):
		self.ct_names = []
		reshape = []
		rx1 = []
		rx2 = []
		bd = []
		#td = []
		#pad = 0
		for contig in self.data:
			reshape.append(self.data[contig])
			
			rx1.append(self.bin_left[contig])
			rx2.append(self.bin_right[contig])
			
			bd.append(self.breadths[contig])
			#td.append(self.tad_values[contig])

			contig_name = self.ct_dict[contig]
			#print(contig_name, np.sum(self.data[contig]))
			
			next_names = [contig_name]*len(self.bin_left[contig])
			
			self.ct_names.extend(next_names)
			
		reshape = np.hstack(reshape)
		
		rx1 = np.hstack(rx1)
		rx2 = np.hstack(rx2)
		bd = np.hstack(bd)
		
		#self.tad_names = self.tad_names * len(td)
		
		#td = np.hstack(td)
		
		self.data = reshape
		
		self.bin_left = rx1
		self.bin_right = rx2
		self.breadths = bd
		#self.tad_values = td
		
		reshape = None
		rx1 = None
		rx2 = None
		bd = None
		td = None
		
		
		#We want to do this before we divide the main data by binwdiths
		self.lower_right_data = np.sum(self.data, axis = 1)
		
		self.binwidths = np.subtract(self.bin_right, self.bin_left)
		
		self.bin_mids = ((self.bin_right + self.bin_left) / 2).astype(int)
		
		#print(self.data.shape)
		#print(len(self.protein_labels))
		
		#print(self.bin_mids)
		#print(np.sum(self.data, axis = 0))
		#print(np.nonzero(np.sum(self.data, axis = 0)))
		
		self.data = np.divide(self.data, self.binwidths[None, :])
		self.data[np.isnan(self.data)] = 0
		
	#Top two chart data
	def top_half(self):
		#Cut the main data into above/below pct ID threshold; summarize
		self.selected = self.y >= self.cutoff
		
		#In group is the area at or above the (default) 95% ID cutoff
		self.upper_left_in = self.data[self.selected, :]
		self.upper_left_in = np.sum(self.upper_left_in, axis = 0)
		#self.upper_left_in = np.divide(self.upper_left_in, self.binwidths)
		#out group is the area below the (default) 95% ID cutoff
		self.upper_left_out = self.data[np.logical_not(self.selected), :]
		self.upper_left_out = np.sum(self.upper_left_out, axis = 0)
		#self.upper_left_out = np.divide(self.upper_left_out, self.binwidths)
		
		#minimum index of true corresponds to maximum y.
		lowest_avail = np.argmax(self.selected)
		#self.tad_in = self.tad_values[lowest_avail]

		self.cov_in = self.breadths[lowest_avail]
		self.tad_in = self.tad_values[lowest_avail]

		#top right chart
		max_obs = max([max(self.upper_left_in), max(self.upper_left_out)])
		
		self.depth_hist_breaks = np.linspace(0, max_obs, num = 199)
		#reuse depth chart info
		
		self.upper_right_in =  np.histogram(self.upper_left_in, bins = self.depth_hist_breaks)[0]
		self.upper_right_out = np.histogram(self.upper_left_out, bins = self.depth_hist_breaks)[0]
		
		#self.upper_right_in = np.histogram(np.log10(self.upper_left_in), bins = self.depth_hist_breaks)[0]
		#self.upper_right_out = np.histogram(self.upper_left_out, bins = self.depth_hist_breaks)[0]
		#self.upper_right_in[self.upper_right_in == 0] = 0.01
		#self.upper_right_out[self.upper_right_out==0] = 0.01
		#self.upper_right_in = np.log10(self.upper_right_in)
		#self.upper_right_out = np.log10(self.upper_right_out)
		
	def make_plots(self):
		print("Plotting", self.mag)
		max_bin_count = self.data.max()
		if max_bin_count < 1: #i.e. max bin == 0, guarantee this with floats by using 1 as the inequality
			#min_bin_count = 1
			print("No best-hit reads were found for", self.mag)
			print("There is no information to plot for this genome under your current read filtering settings. This genome will be skipped.")
			
		else:
			min_bin_count = self.data[np.nonzero(self.data)].min()

			#order of magnitude
			#smallest oom is this
			current_oom = -len(str(int(1/min_bin_count))) + 1
			
			#current_oom = 0
			tick_positions = []
			
			while min_bin_count < max_bin_count:
				tick_positions.append(current_oom)
				min_bin_count *= 10
				current_oom += 1

			tick_labels = []
			for oom in tick_positions:
				tick_labels.append("10e"+str(oom))
			
			overall_plot =  make_subplots(rows=3, cols=2, 
										column_widths=[0.66, 0.34], 
										row_heights = [0.30, 0.10, 0.60],

										shared_xaxes=True,
										shared_yaxes=True,
										horizontal_spacing = 0.025,
										vertical_spacing = 0.035)
			
			#This one is strange
			#The protein hovertext needs multiple lines, but the customdata
			#arg doesn't appear to work for customdata of more than 1 row
			#Therefore we just <br>.join() a list of text as the labels.
			annot_hov = "{text}"
				
			#Protein independent hover templates
			bot_left_hov = 	"Position in Genome: %{x:0d}<br>" +\
							"Percent Identity: %{y:}%<br>" +\
							"Log 10 Base Count: %{z:.4f}<br>" +\
							"<extra></extra>"
									
			#hover data templates
			in_dep_hov = "Position in Genome: %{x:0d}<br>" +\
						"Avg. Within-Pop. Depth: %{y:.2f}" +\
						"<extra></extra>"

			out_dep_hov = "Position in Genome: %{x:0d}<br>" +\
						"Avg. Outside-Pop. Depth: %{y:.2f}" +\
						"<extra></extra>"
										
			#Protein independent hover templates
			bot_right_hov = "Total Bases: %{x:0d}<br>" +\
								"Percent Identity: %{y:}%<br>" + \
								"<extra></extra>"
								
								
			#TAD *could* have protein info, but I don't think it bears a fourth repetition
			tad_hov = "Contig: %{x:c}<br>" +\
						"TAD Level: %{y:c}<br>" +\
						"Depth: %{z:.2f}" +\
						"<extra></extra>"
						
			
			top_right_hov_in = "Count of Obs.: %{x:0d}<br>" +\
			"Avg. Within-Pop. Depth: %{y:.4f}<br>" +\
			"<extra></extra>"
			
			top_right_hov_out = "Count of Obs.: %{x:0d}<br>" +\
			"Avg. Outside-Pop. Depth: %{y:.4f}<br>" +\
			"<extra></extra>"
			
			#bot left main plot
			#Needs to have labels added.
			overall_plot.add_trace(go.Heatmap(z=np.log10(self.data), 
									x = self.bin_mids,
									y = self.y,
									colorbar = dict(
										x = -0.1,
										tickvals = tick_positions,
										ticktext = tick_labels,
										tickfont = {"size": self.axis_font_size}
										),
									colorscale=self.main_plot_fill_colorscale,
									hovertemplate = bot_left_hov
									),
									row = 3, col = 1)
			
			#bot right - only one of these
			overall_plot.add_trace(go.Bar(x = self.lower_right_data, 
									y = self.y,
									orientation='h',
									marker = dict(color = self.bot_right_line_col),
									hovertemplate = bot_right_hov
									),
									row = 3, col = 2)
			
			#text_rec = []
			
			#We just use the protein labels as a repo - normally it would be filled out
			if not self.do_prot:
				self.protein_labels = []
				for c, s, e in zip(self.ct_names, self.bin_left, self.bin_right):
					next_label = "<br>".join(["Contig: " + c,
											"Genome Region: " + str(s) + "-" + str(e)
											])
					self.protein_labels.append(next_label)
			
			overall_plot.add_trace(go.Scatter(x = self.bin_mids,
									y = [1] * len(self.bin_mids),
									text = self.protein_labels,
									marker = dict(color = self.in_group_depth_col),
									hoverinfo = "text"
									),
									row = 2, col = 1)
									
			which_viz = None
			which_id = None
			has_default = False
			group = 0
			
			starting_step = len(overall_plot.data)
			
			step_groups = {}
			id_grp = {}
			
			tad_ticks = np.linspace(0, self.tad_max, num = 5, dtype = int).tolist()
			#tad_labels = 
			tad_colorbar = dict(tickvals = tad_ticks, 
								tickfont = {"size" : self.axis_font_size})
			
			#Here, we iterate over the in-groups to add traces for each pct id in-group.
			for pct_id_cutoff in self.y:
				self.cutoff = pct_id_cutoff
				#recalculate data
				self.top_half()
				
				idname = str(pct_id_cutoff)
				
				#print(pct_id_cutoff, self.tad_in, self.cov_in)
				
				step_groups[pct_id_cutoff] = [starting_step,
											starting_step + 1,
											starting_step + 2,
											starting_step + 3,
											starting_step + 4]
											#starting_step + 5,]
											
				starting_step += 5
												
				#Set the only default line *near* 95 pct id.
				if self.cutoff >= 95.0 and not has_default:
					which_viz = group
					which_id = pct_id_cutoff
					has_default = True
				else:
					group += 1
				
				line_height = self.cutoff-(self.id_step/2)
				

				#top left out group							
				overall_plot.add_trace(go.Scatter(x = self.bin_mids,
												y = self.upper_left_out, 
												marker = dict(color = self.out_group_depth_col),
												visible = False,
												hovertemplate = out_dep_hov	
												), 
												row = 1, col = 1)
												
				#top left in group - this goes second so that it's plotted on top in the layers
				overall_plot.add_trace(go.Scatter(x = self.bin_mids, 
												y = self.upper_left_in,
												marker = dict(color = self.in_group_depth_col),
												visible = False,
												hovertemplate = in_dep_hov
												),
												row = 1, col = 1)

				#mid right tad
				#Update fixed colorbar, add the coverage values and this is done.
				overall_plot.add_trace(go.Heatmap(z = self.tad_in,
												x = self.tad_contigs,
												y = self.tad_names,
												visible = False,
												hovertemplate = tad_hov,
												colorbar = tad_colorbar,
												zmin = 0,
												zmax = self.tad_max*1.05
												),
												row = 2, col = 2)
				
				#Prevent 0-depth outgroups from running away with the scale
				#max_in_grp = self.upper_right_in.max() + 1
				#self.upper_right_out[self.upper_right_out > max_in_grp] = max_in_grp
				

				#top right out group default
				overall_plot.add_trace(go.Scatter(x = self.upper_right_out, 
				y = self.depth_hist_breaks, 
				marker = dict(color = self.out_group_depth_col),
				visible = False,
				hovertemplate = top_right_hov_out
				), 
				row = 1, col = 2)
				
				#top right in group default - same ordering issue as above.
				overall_plot.add_trace(go.Scatter(x = self.upper_right_in, 
				y = self.depth_hist_breaks, 
				marker = dict(color = self.in_group_depth_col),
				visible = False,
				hovertemplate = top_right_hov_in
				), 
				row = 1, col = 2)
				
			#Set default lines as nearest to 90
			for i in step_groups[which_id]:
				overall_plot.data[i].visible = True
				
			#in-group highlight
			overall_plot.add_hrect(y0 = which_id-(self.id_step/2), 
									y1 = 100+(self.id_step/2),
									row = 3, col = 1,
									line_width=0,
									#line_width=5,
									fillcolor=self.main_plot_highlight_col,
									#line_color=self.main_plot_highlight_col,
									#line_dash="dash")
									opacity=self.main_plot_highlight_alpha)
									#,layer="below")
									
			overall_plot.add_hrect(y0 = which_id-(self.id_step/2), 
									y1 = 100+(self.id_step/2),
									row = 3, col = 2,
									line_width=0,
									#line_width=5,
									fillcolor=self.main_plot_highlight_col,
									#line_color=self.main_plot_highlight_col,
									#line_dash="dash")
									opacity=self.main_plot_highlight_alpha)
									#,layer="below")
									
			#print(overall_plot.layout['shapes'])
				
			
			left_box = {'type': 'rect', 
						'x0': 0, 
						'x1': 1, 
						'xref': 
						'x3 domain', 
						'y0': 0, 
						'y1': 1, 
						'yref': 'y3 domain'}
			
			right_box = {'type': 'rect', 
						'x0': 0, 
						'x1': 1, 
						'xref': 
						'x4 domain', 
						'y0': 0, 
						'y1': 1, 
						'yref': 'y4 domain'}
						
						
						
			initial_shapes = list(overall_plot.layout["shapes"])
			initial_shapes.append(left_box)
			initial_shapes.append(right_box)
			overall_plot.layout["shapes"] = tuple(initial_shapes)
			
			steps = []
			for group in step_groups:
				step = dict(
						method="update",
						args=[{"visible": [False] * len(overall_plot.data)},
							{"shapes": [{'line': {'width': 0, 'color' : self.main_plot_highlight_col, 'dash' : 'dash'},
										'fillcolor': self.main_plot_highlight_col,
										#'layer': 'below',
										'opacity': self.main_plot_highlight_alpha,
										'type': 'rect',
										'x0': 0,
										'x1': 1,
										'xref': 'x5 domain',
										'y0': group-(self.id_step/2),
										'y1': 100+(self.id_step/2),
										'yref': 'y5'},
										
										#lower right highlight
										{'line': {'width': 0, 'color' : self.main_plot_highlight_col, 'dash' : 'dash'},
										'fillcolor': self.main_plot_highlight_col,
										#'layer': 'below',
										'opacity': self.main_plot_highlight_alpha,
										'type': 'rect',
										'x0': 0,
										'x1': 1,
										'xref': 'x6 domain',
										'y0': group-(self.id_step/2),
										'y1': 100+(self.id_step/2),
										'yref': 'y6'},
										
										#mid left box
										left_box,
										
										#mid right box
										right_box
										
										]
							}],
						
						label = str(round(group, 2))
				)
				#lower left and lower right data
				step["args"][0]["visible"][0] = True #lower left
				step["args"][0]["visible"][1] = True #lower right
				step["args"][0]["visible"][2] = True #annotation mid left
				for i in step_groups[group]:
					step["args"][0]["visible"][i] = True
					
				steps.append(step)
				
			id_slider = [dict(
				#match default viz
				active=which_viz,
				currentvalue={"prefix": "Percent ID cutoff: ", "suffix": "%"},
				pad={"t": 50},
				bgcolor = self.in_group_depth_col,
				steps=steps,
				#tickfont = {"size" : self.axis_font_size},
				font = {"size" : self.axis_font_size}
			)]
				
				
			overall_plot['layout']['xaxis2'].pop('matches')

			overall_plot['layout']['xaxis2']['showticklabels'] = True
			
			overall_plot['layout']['yaxis3']['showticklabels'] = False
			overall_plot['layout']['yaxis3'].pop('matches')
			
			overall_plot['layout']['xaxis4'].pop('matches')
			overall_plot['layout']['yaxis4'].pop('matches')
			overall_plot['layout']['yaxis4']['showticklabels'] = False
			
			
			#print(overall_plot['layout'])
			
			overall_plot.update_layout(showlegend = False)
			overall_plot.update_layout(margin = dict(t=25))
			overall_plot.update_xaxes(showgrid=False)
			overall_plot.update_yaxes(showgrid=False)
			
			overall_plot.update_xaxes(tickfont = {"size" : self.axis_font_size})
			overall_plot.update_yaxes(tickfont = {"size" : self.axis_font_size})
			
			#overall_plot.update_layout(hovermode="x unified")
			
			overall_plot.update_layout(
				sliders=id_slider,
				plot_bgcolor='#f2f2f2'
			)
			
			
			#print(overall_plot.data[0])
			
			#relayer = list(overall_plot.data)[2:]
			#relayer.append(overall_plot.data[0])
			#relayer.append(overall_plot.data[1])
			#relayer = tuple(relayer)
			#overall_plot.data = relayer
			#relayer = None
			
			#print(overall_plot.data[0])
			#print(overall_plot.data[-2])
			
			overall_plot.write_html(self.plot_name, config = self.html_config)
		
		
	def build(self):
		self.bin_raw()
		self.concatenate()
		self.make_plots()
		

def plot_opts():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
			description='''
	''')
	parser.add_argument('-d', '--database',  dest = 'db', default = None, 
	help =  'Path to the RecruitPlotEasy database you want to make plots from. Required.')
	
	parser.add_argument('-p', '--proteins',  dest = 'prot', action='store_true', 
	help =  "Plot proteins. You had to predict proteins in the database build, or this won't work. Overrides --width")
	
	parser.add_argument('-w', '--width',  dest = 'width', type=int, default = 1000, 
	help =  'Genome bin width. The genome will be divided into bins of [width] size for plotting. Def. 1000')
	
	parser.add_argument('-i', '--id_step',  dest = 'height', type=float, default = 0.5, 
	help =  'Pct. ID bin height. Reads will be binned every [height] pct. ID. Default 0.5.')
	
	parser.add_argument('-o', '--output',  dest = 'outdir', default = "recruitment_plots", 
	help =  'Base directory for outputs. Defaults to creating "recruitment_plots/" in the CWD.')
	
	parser.add_argument('-m', '--mag',  dest = 'mag', default = None, 
	help =  'Plot only this MAG and exit.')
	
	parser.add_argument('--mag_file',  dest = 'mf', default = None, 
	help =  'A file containing the names of MAGs to plot, one per line. Only these MAGs will be plotted.')
	
	args, unknown = parser.parse_known_args()
	
	return parser, args
		
def run_plot():
	parser, opts = plot_opts()
	db = opts.db
	if db is None:
		print("You must supply a database. Quitting.")
		parser.print_help()
		sys.exit()
	else:
		if not os.path.exists(db):
			print("Database", db, "could not be found.")
			print("Have you created one with RecruitPlotEasy build yet?  Quitting.")
			sys.exit()
	
	do_proteins = opts.prot
	
	width = opts.width
	height = opts.height
	
	out = opts.outdir
	
	#Mag files
	mag = opts.mag
	mf = opts.mf
	
	mn = rpdb(db, gen_step = width, id_step = height, do_prot = do_proteins, output_base = out)
	mn.open()
	mn.parse_db()
	
	try:
		mags_to_plot = list(mn.genomes.keys())
	except:
		mags_to_plot = None
	
	if mag is not None:
		#check if the MAG is in the db.
		if mag in mn.genomes:
			print("MAG", mag, "detected. Plotting...")
			mags_to_plot = [mag]
		else:
			mags_to_plot = None
			print("Requested MAG:", mag, "was not found in the database. It cannot be plotted.")
			print("Consider RecruitPlotEasy's describe module to take a look at your databases' contents.")
			
	if mf is not None:
		mags_to_plot = []
		with open(mf) as fh:
			for line in fh:
				mag = line.strip()
				mags_to_plot.append(mag)
	
	if mags_to_plot is not None:
		if len(mags_to_plot) == 0:
			mags_to_plot = None
		else:
			mags_to_plot = ['"'+m+'"' for m in mags_to_plot]
			mags_to_plot = set(mags_to_plot)
	
	if mags_to_plot is None:
		print("No valid MAGs or genomes were detected for plotting.\nRecruitPlotEasy needs at least one to proceed.\n\nQuitting.")
		sys.exit()
	
	if do_proteins:
		if mn.proteins is None:
			print("No proteins have been added to this database. Quitting.")
			mn.close()
			sys.exit()
			
	if do_proteins:
		print("Plotting proteins with pct ID height", height)
	else:
		print("Plotting genome with bin width ", width, 'bases and pct ID height ', height)
	
	
	for sample in mn.samples:
		mn.set_sample(sample)
		for mag in mn.mags_in_sample:
			if mag in mags_to_plot:
				mn.set_mag(mag)
				mn.craft_query()
				mn.load_sample()
		
	mn.close()







