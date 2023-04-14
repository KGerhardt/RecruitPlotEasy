import os
class file_checker:
	def __init__(self, file):
		self.file = file
		self.format = "Unknown"
		
		self.original_name = None
		self.sqlname = None
		self.get_sql_name()
	
	def sql_safe(self, string):
		#Sanitize for SQL
		#These are chars safe for sql
		sql_safe = set('_abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789')
		current_chars = set(string)
		#self.sql_name = self.basename
		#Identify SQL-unsafe characters as those outside the permissible set and replace all with underscores.
		for char in current_chars - sql_safe:
			string = string.replace(char, "_")
		
		return string
	
	
	def get_sql_name(self):
		base = os.path.basename(self.file)
		self.original_name = base
		#while base != os.path.splitext(base)[0]:
		##	base = os.path.splitext(base)[0]
		
		self.sqlname = self.sql_safe(self.original_name)
		
	def check_is_fasta(self):
		fh = open(self.file, "r")
		fmt_fine = True
		for i in range(0, 30):
			try:
				line = fh.readline()
			except:
				fmt_fine = False
				break
			if not line.startswith(">"):
				if not set(line) <= set("ATCGatcgNn\n"):
					fmt_fine = False
		fh.close()
		
		if fmt_fine:
			self.format = "fasta"
		
	def check_is_blast(self):
		fh = open(self.file, "r")
		fmt_fine = True
		try:
			line = fh.readline().strip()
		except:
			fmt_fine = False
		else:
			while line.startswith("#"):
				line = fh.readline().strip()
			segment = line.split()
			if len(segment) < 9:
				fmt_fine = False
			else:
				if segment[8].isnumeric() and segment[9].isnumeric():
					fmt_fine = True
				else:
					fmt_fine = False
		fh.close()
		
		if fmt_fine:
			self.format = "blast"
		
	def check_is_sam(self):
		fh = open(self.file, "r")
		fmt_fine = True
		try:
			line = fh.readline()
		except :
			fmt_fine = False
		else:
			#Hm. This might flag on a fastq
			if line.startswith("@"):
				pass
			else:
				fmt_fine = False
		fh.close()
		
		if fmt_fine:
			self.format = "sam"
			
	def check_is_bam(self):
		try:
			fh = open(self.file, "rb")
			fmt_fine = True
			head = fh.read(4)
			if head != b'BAM\1' and head != b"\x1f\x8b\x08\x04":
				fmt_fine = False
			fh.close()
		except:
			fmt_fine = False
		if fmt_fine:
			self.format = "bam"
	
	def run_checks(self):
		self.check_is_fasta()
		self.check_is_blast()
		self.check_is_sam()
		self.check_is_bam()
		