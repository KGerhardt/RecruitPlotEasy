from itertools import groupby

class cigar_to_pos:
	def __init__(self, cigar, start):
		self.cigar = cigar
		self.start = start
		#self.consumes_read = set("M", "I", "S", "=", "X")
		self.consumes_ref = set(["M", "D", "N", "=", "X"])
		self.last = self.start
		self.covers = []
		
	def parse_cig_sam(self):
		splits = groupby(self.cigar, lambda c: c.isdigit())
		for g, n in splits:
			num_bases = int("".join(n))
			op = "".join(next(splits)[1])
			if op in self.consumes_ref:
				#Python doesn't include the last base here, so we pad the end of the cover by 1.
				self.covers.append((self.last, self.last + num_bases + 1,))
				self.last += num_bases
				
	def parse_cig_bam(self):
		for tup in self.cigar:
			num_bases = tup[0]
			op = tup[1]
			if op in self.consumes_ref:
				#Python doesn't include the last base here, so we pad the end of the cover by 1.
				self.covers.append((self.last, self.last + num_bases + 1,))
				self.last += num_bases
				
	def get_ranges(self):
		return self.covers