import pyrodigal as pd

class pyrodigal_manager:
	def __init__(self, file = None, aa_out = None, nt_out = None, is_meta = False, full_headers = True, trans_table = 11,
				num_bp_fmt = True, verbose = True, do_compress = "0", compare_against = None):
		#Input NT sequences
		self.file = file
		
		#List of seqs read from input file.
		self.sequences = None
		#Concatenation of up to first 32 million bp in self.sequences - prodigal caps at this point.
		self.training_seq = None
		
		#Predicted genes go here
		self.predicted_genes = None
		#Record the translation table used.
		self.trans_table = trans_table
		
		#This is the pyrodigal manager - this does the gene predicting.
		self.manager = pd.OrfFinder(meta=is_meta)
		self.is_meta = is_meta
		
		#Full prodigal header information includes more than just a protein number.
		#If full_headers is true, protein deflines will match prodigal; else, just protein ID.
		self.full_headers = full_headers
		
		#Prodigal prints info to console. I enhanced the info and made printing default, but also allow them to be totally turned off.
		self.verbose = verbose
		
		#Prodigal formats outputs with 70 bases per line max
		self.num_bp_fmt = num_bp_fmt
		
		#File names for outputs
		self.aa_out = aa_out
		self.nt_out = nt_out
		
		#List of proteins in excess of 100K base pairs (HMMER's limit) and their lengths. This is also fastAAI specific.
		self.excluded_seqs = {}
		
		#Gzip outputs if asked.
		self.compress = do_compress
		
		self.labeled_proteins = None
		
		self.annot = None
		self.source_gen = None
		
		#Normally, we don't need to keep an input sequence after it's had proteins predicted for it - however
		#For FastAAI and MiGA's purposes, comparisons of two translation tables is necessary.
		#Rather than re-importing sequences and reconstructing the training sequences, 
		#keep them for faster repredict with less I/O
		self.compare_to = compare_against
		if self.compare_to is not None:
			self.keep_seqs = True
			self.keep_after_train = True
		else:
			self.keep_seqs = False
			self.keep_after_train = False
	
	#Imports a fasta as binary.
	def import_sequences(self):
		if self.sequences is None:
			self.sequences = {}
			
		#check for zipped and import as needed.
		with open(self.file, 'rb') as test_gz:
			#Gzip magic number
			is_gz = (test_gz.read(2) == b'\x1f\x8b')
		
		if is_gz:
			fh = gzip.open(self.file)
		else:
			fh = open(self.file, "rb")
		
		imp = fh.readlines()
		
		fh.close()
		
		cur_seq = None
		for s in imp:
			s = s.decode().strip()
			#> is 62 in ascii. This is asking if the first character is '>'
			if s.startswith(">"):
				#Skip first cycle, then do for each after
				if cur_seq is not None:
					self.sequences[cur_seq] = ''.join(self.sequences[cur_seq])
					self.sequences[cur_seq] = self.sequences[cur_seq].encode()
					#print(cur_seq, len(self.sequences[cur_seq]))
				cur_seq = s[1:]
				cur_seq = cur_seq.split()[0]
				cur_seq = cur_seq.encode('utf-8')
				self.sequences[cur_seq] = []
			else:
				#Remove the newline character.
				#bases = s[:-1]
				self.sequences[cur_seq].append(s)
		
		#Final set
		self.sequences[cur_seq] = ''.join(self.sequences[cur_seq])
		self.sequences[cur_seq] = self.sequences[cur_seq].encode()
		
		#Now we have the data, go to training.
		if not self.is_meta:
			self.train_manager()
		
	def convert_sequences(self, contents):
		self.sequences = {}
		for seq in contents:
			name = seq.encode()
			nt = contents[seq].encode()
			self.sequences[name] = nt
		
	#Collect up to the first 32 million bases for use in training seq.
	def train_manager(self):
		running_sum = 0
		seqs_added = 0
		if self.training_seq is None:
			self.training_seq = []
			for seq in self.sequences:
				running_sum += len(self.sequences[seq])
				if seqs_added > 0:
					#Prodigal interleaving logic - add this breaker between sequences, starting at sequence 2
					self.training_seq.append(b'TTAATTAATTAA')
					running_sum += 12
					
				seqs_added += 1
					
				#Handle excessive size
				if running_sum >= 32000000:					
					print("Warning:  Sequence is long (max 32000000 for training).")
					print("Training on the first 32000000 bases.")
				
					to_remove = running_sum - 32000000
					
					#Remove excess characters
					cut_seq = self.sequences[seq][:-to_remove]
					#Add the partial seq
					self.training_seq.append(cut_seq)
					
					#Stop the loop and move to training
					break
				
				#add in a full sequence
				self.training_seq.append(self.sequences[seq])

			if seqs_added > 1:
				self.training_seq.append(b'TTAATTAATTAA')
				
			self.training_seq = b''.join(self.training_seq)
		
		if len(self.training_seq) < 20000:
			if self.verbose:
				print("Can't train on 20 thousand or fewer characters. Switching to meta mode.")
			self.manager = pd.OrfFinder(meta=True)
			self.is_meta = True
		else:
			if self.verbose:
				print("")
				#G is 71, C is 67; we're counting G + C and dividing by the total.
				gc = round(((self.training_seq.count(67) + self.training_seq.count(71))/ len(self.training_seq)) * 100, 2)
				print(len(self.training_seq), "bp seq created,", gc, "pct GC")
				
			#Train
			self.manager.train(self.training_seq, translation_table = self.trans_table)
		
		if not self.keep_after_train:
			#Clean up
			self.training_seq = None
		
	def predict_genes(self):
		if self.is_meta:
			if self.verbose:
				print("Finding genes in metagenomic mode")
		else:
			if self.verbose:
				print("Finding genes with translation table", self.trans_table)
				print("")
			
		self.predicted_genes = {}
		for seq in self.sequences:
			
			if self.verbose:
				print("Finding genes in sequence", seq.decode(), "("+str(len(self.sequences[seq]))+ " bp)... ", end = '')
				
			self.predicted_genes[seq] = self.manager.find_genes(self.sequences[seq])
				
			#If we're comparing multiple tables, then we want to keep these for re-prediction.
			if not self.keep_seqs:
				#Clean up
				self.sequences[seq] = None
			
			if self.verbose:
				print("done!")
			
	#Predict genes with an alternative table, compare results, and keep the winner.	
	def compare_alternative_table(self, table):
		if table == self.trans_table:
			print("You're trying to compare table", table, "with itself.")
		else:
			if self.verbose:
				print("Comparing translation table", self.trans_table, "against table", table)
			old_table = self.trans_table
			old_genes = self.predicted_genes
			old_size = 0
			for seq in self.predicted_genes:
				for gene in self.predicted_genes[seq]:
					old_size += (gene.end - gene.begin)
			
			self.trans_table = table
			self.train_manager()
			self.predict_genes()
				
			new_size = 0
			for seq in self.predicted_genes:
				for gene in self.predicted_genes[seq]:
					new_size += (gene.end - gene.begin)
			
			if (old_size / new_size) > 1.1:
				if self.verbose:
					print("Translation table", self.trans_table, "performed better than table", old_table, "and will be used instead.")
			else:
				if self.verbose:
					print("Translation table", self.trans_table, "did not perform significantly better than table", old_table, "and will not be used.")
				self.trans_table = old_table
				self.predicted_genes = old_genes
			
			#cleanup
			old_table = None
			old_genes = None
			old_size = None
			new_size = None
		
	def predict_and_compare(self):
		self.predict_genes()
	
		#Run alt comparisons in gene predict.
		if self.compare_to is not None:
			while len(self.compare_to) > 0:
				try:
					next_table = int(self.compare_to.pop(0))
					
					if len(self.compare_to) == 0:
						#Ready to clean up.
						self.keep_after_train = True
						self.keep_seqs = True
					
					self.compare_alternative_table(next_table)
				except:
					print("Alternative table comparison failed! Skipping.")
		
	#Break lines into size base pairs per line. Prodigal's default for bp is 70, aa is 60.
	def num_bp_line_format(self, string, size = 70):
		#ceiling funciton without the math module
		ceiling = int(round((len(string)/size)+0.5, 0))
		formatted = '\n'.join([string[(i*size):(i+1)*size] for i in range(0, ceiling)])
		return formatted
	
	def translate(self):
		renamed = {}
		self.annot = {}
		self.source_gen = {}
		for seq in self.predicted_genes:
			count = 1
			source = seq.decode()
			seqname = source + "_"
			
			for gene in self.predicted_genes[seq]:
				gene_name = seqname + str(count)
				count += 1
				renamed[gene_name] = gene.translate()
				#print(dir(gene))
				annotation = ";".join(["GC:" + str(round(gene.gc_cont, 4)),
										"Start_Type:" +gene.start_type])
				#print(gene.strand, gene.begin, gene.end, str(gene._gene_data),)
				this_annot = (gene.strand, gene.begin, gene.end, annotation,)
				self.annot[gene_name] = this_annot
				self.source_gen[gene_name] = source
				
		self.predicted_genes = renamed
	
	#Writeouts
	def write_nt(self):
		if self.nt_out is not None:
			if self.verbose:
				print("Writing nucleotide sequences... ")
			if self.compress == '1' or self.compress == '2':
				out_writer = gzip.open(self.nt_out+".gz", "wb")
				
				content = b''
				
				for seq in self.predicted_genes:
					seqname = b">"+ seq + b"_"
					#Gene counter
					count = 1
					for gene in self.predicted_genes[seq]:
						#Full header lines
						if self.full_headers:
							content += b' # '.join([seqname + str(count).encode(), str(gene.begin).encode(), str(gene.end).encode(), str(gene.strand).encode(), gene._gene_data.encode()])
						else:
							#Reduced headers if we don't care.
							content += seqname + str(count).encode()
							
						content += b'\n'
							
						if self.num_bp_fmt:
							#60 bp cap per line
							content += self.num_bp_line_format(gene.sequence(), size = 70).encode()
						else:
							#One-line sequence.
							content += gene.sequence().encode()
							
						content += b'\n'
						count += 1
				
				out_writer.write(content)
				out_writer.close()
			
			if self.compress == '0' or self.compress == '2':
				out_writer = open(self.nt_out, "w")
			
				for seq in self.predicted_genes:
					#Only do this decode once.
					seqname = ">"+ seq.decode() +"_"
					#Gene counter
					count = 1
					
					for gene in self.predicted_genes[seq]:
						#Full header lines
						if self.full_headers:
							#Standard prodigal header
							print(seqname + str(count), gene.begin, gene.end, gene.strand, gene._gene_data, sep = " # ", file = out_writer)
						else:
							#Reduced headers if we don't care.
							print(seqname + str(count), file = out_writer)
							
						if self.num_bp_fmt:
							#60 bp cap per line
							print(self.num_bp_line_format(gene.sequence(), size = 70), file = out_writer)
						else:
							#One-line sequence.
							print(gene.sequence(), file = out_writer)
							
						count += 1
							
				out_writer.close()
		
	def write_aa(self):
		if self.aa_out is not None:
			if self.verbose:
				print("Writing amino acid sequences...")
				
			self.labeled_proteins = {}
			content = ''
			for seq in self.predicted_genes:
				count = 1
				seqname = ">"+ seq.decode() + "_"
				for gene in self.predicted_genes[seq]:
					prot_name = seqname + str(count)
					translation = gene.translate()
					self.labeled_proteins[prot_name[1:]] = translation
					defline = " # ".join([prot_name, str(gene.begin), str(gene.end), str(gene.strand), str(gene._gene_data)])
					content += defline
					content += "\n"
					count += 1
					content += self.num_bp_line_format(translation, size = 60)
					content += "\n"					
				
			if self.compress == '0' or self.compress == '2':
				out_writer = open(self.aa_out, "w")
				out_writer.write(content)
				out_writer.close()
				
			if self.compress == '1' or self.compress == '2':
				content = content.encode()
				out_writer = gzip.open(self.aa_out+".gz", "wb")
				out_writer.write(content)
				out_writer.close()
				
	def run_for_fastaai(self):
		self.verbose = False
		self.import_sequences()
		self.train_manager()
		self.predict_and_compare()
		self.write_aa()
	
#Testing section.	
#import sys
#mn = pyrodigal_manager(file = sys.argv[1])
#mn.import_sequences()
#mn.train_manager()
#mn.predict_genes()
#mn.translate()
