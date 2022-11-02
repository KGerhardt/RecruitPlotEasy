import re
class mdz_to_match_count:
	def __init__(self, mdz):
		self.mdz = mdz
		self.matches = 0
		self.mismatches = 0
		self.parse_mdz()
	
	def parse_mdz(self):
		match_count = re.findall('[0-9]+', self.mdz)
		for num in match_count:
			self.matches += int(num)
			
		self.mismatches = len([i for i in self.mdz if i.isalpha()]) 
	
	def get_match_mismatch(self):
		return self.matches, self.mismatches