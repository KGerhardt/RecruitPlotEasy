import sys
import os
import re
import numpy as np

#from .agnostic_reader import agnostic_reader
from .agnostic_reader import agnostic_reader
from .cig_parser import cigar_to_pos
from .mdz_parser import mdz_to_match_count

class sp_iter:
	def __init__(self, file):
		self.handle = agnostic_reader(file)

	def parse_sam_entry(self, line):
		if "MD:Z:" not in line:
			return None
			
		segs = line.strip().split("\t")
		query = segs[0]
		
		#Unaligned read.
		if query is None or query == "":
			return None
		
		target = segs[2]
		#Sam is 1-indexed. For python, we want 0.
		start = int(segs[3]) - 1
		cig = segs[5]
		read_len = len(segs[9])
		covered_ranges = cigar_to_pos(cig, start)
		covered_ranges.parse_cig_sam()
		covered_ranges = covered_ranges.get_ranges()
		#Find the MD:Z: tag
		iter = len(segs)-1
		mdz_seg = segs[iter]
		# If it's not the correct field, proceed until it is.
		while not mdz_seg.startswith("MD:Z:"):
			iter -= 1
			mdz_seg = segs[iter]
		#Remove the MD:Z: flag from the start
		mdz_seg = mdz_seg[5:]
		#Get num matches
		parsed_mds = mdz_to_match_count(mdz_seg)
		matches, mismatches = parsed_mds.get_match_mismatch()
		
		total_align = matches+mismatches
		
		pct_id_local = (matches / total_align) * 100
		#The second item in the last tuple of covered range is the last pos covered + 1, adjust it back
		last_pos_aligned = covered_ranges[len(covered_ranges)-1][1] - 1
		#We don't need to pad the last pos covered
		pct_id_global = (matches / read_len) * 100
		
		pct_alignment = (total_align / read_len) * 100
		'''
		as_np_list = []
		for tuple in covered_ranges:
			for index in tuple:
				as_np_list.append(index)
		covered_ranges = np.array(as_np_list, dtype = np.int32)
		'''
		return [query, target, pct_id_local, pct_id_global, pct_alignment, covered_ranges]
		
	def __next__(self):
		#Skip header.
		line = self.handle.readline()
		while line.startswith("@"):
			line = self.handle.readline()
		
		if line:
			line = self.parse_sam_entry(line)
			return line
		else:
			self.handle.close()
			raise StopIteration

class sam_parser:
	def __init__(self, file):
		self.file = file
		self.header = []
		self.handle = None
		self.get_header()
		
	def get_header(self):
		self.handle = agnostic_reader(self.file)
		line = self.handle.readline()
		while line.startswith("@"):
			if line.startswith("@SQ"):
				self.header.append(line)
			line = self.handle.readline()
		self.handle.close()
		self.handle = None
		
	def __iter__(self):
		return sp_iter(self.file)
			
