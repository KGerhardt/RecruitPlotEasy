from .agnostic_reader import agnostic_reader
import sys
import numpy as np

class blast_iter:
	def __init__(self, file):
		self.handle = agnostic_reader(file)
		
	def parse_blast_record(self, blast_record):
		try:
			blast_record = blast_record.strip().split("\t")
			#Outfmt 6 - query, target, %ID, aln length, mismatch, query aln st, query aln end, target aln start (the interesting one), target aln end, e-value, bit score
			#outfmt6 doesn't include seqlen.
			query = blast_record[0]
			target = blast_record[1]
			pct_id_local = float(blast_record[2])
			
			q_aln_1 = int(blast_record[6]) 
			q_aln_2 = int(blast_record[7]) 
			
			ref_aln_1 = int(blast_record[8]) 
			ref_aln_2 = int(blast_record[9])
			
			#Blast specifies read alignment start and end relative to the direction of the read
			#i.e. start < end for reads aligning to complement
			query_start = min([q_aln_1, q_aln_2])
			query_end = max([q_aln_1, q_aln_2])
			
			#Same for ref.
			ref_start = min([ref_aln_1, ref_aln_2])
			ref_end = max([ref_aln_1, ref_aln_2])

			#Padded by 1 or a position is lost
			#alignment_length = ref_end - ref_start + 1
			
			#These is not a calculable value for a standard blastn -outfmt6
			#They require read length, which is NOT included.
			pct_id_global = -1.0
			pct_alignment = -1.0
			
			covered_ranges = [(ref_start, ref_end+1, )]
			'''
			as_np_list = []
			for tuple in covered_ranges:
				for index in tuple:
					as_np_list.append(index)
			covered_ranges = np.array(as_np_list, dtype = np.int32)
			'''
			results = [query, target, pct_id_local, pct_id_global, pct_alignment, covered_ranges]
		except:
			results = None
			
		return results
		
	def __next__(self):
		#Skip header.
		blast_record = self.handle.readline()
		#Handle magicblast header
		while blast_record.startswith("#"):
			blast_record = self.handle.readline()
		
		if blast_record:
			line = self.parse_blast_record(blast_record)
			return line
		else:
			raise StopIteration
		
class blast_parser:
	def __init__(self, file):
		self.file = file
		
	def __iter__(self):
		return blast_iter(self.file)


