#!/usr/bin/env python3

import sys
import re
import bisect
import sqlite3
import argparse
import platform
import importlib
import os
import zlib
import time
import tempfile
import subprocess
from array import array
from struct import unpack
from urllib.request import pathname2url

if platform.system() != "Windows":
	try:
		lib = importlib.import_module("pysam")
	except:
		print("", end = "")
	else:
		import pysam
	
def get_sys():
	return(platform.system())

CtoPy	   = { 'A':'<c', 'c':'<b', 'C':'<B', 's':'<h', 'S':'<H', 'i':'<i', 'I':'<I', 'f':'<f' }
py4py	   = { 'A':  1 , 'c':  1 , 'C':  1 , 's':  2 , 'S':  2 , 'i':  4 , 'I':  4 , 'f':  4  }
dna_codes   = '=ACMGRSVTWYHKDBN'
cigar_codes = 'MIDNSHP=X'
parse_codes = {
	'sam':					 ' The current alignment in SAM format.',
	'bam':					 ' All the bytes that make up the current alignment ("read"),\n							  still in binary just as it was in the BAM file. Useful\n							  when creating a new BAM file of filtered alignments.',
	'sam_qname':			   ' [1st column in SAM] The QNAME (fragment ID) of the alignment.',
	'bam_qname':			   ' The original bytes before decoding to sam_qname.',
	'sam_flag':				' [2nd column in SAM] The FLAG number of the alignment.',
	'bam_flag':				' The original bytes before decoding to sam_flag.',
	'sam_refID':			   ' The chromosome ID (not the same as the name!).\n							  Chromosome names are stored in the BAM header (file_chromosomes),\n							  so to convert refIDs to chromsome names one needs to do:\n							  "my_bam.file_chromosomes[read.sam_refID]" (or use sam_rname)\n							  But for comparisons, using the refID is much faster that using\n							  the actual chromosome name (for example, when reading through a\n							  sorted BAM file and looking for where last_refID != this_refID)\n							  Note that when negative the alignment is not aligned, and thus one\n							  must not perform my_bam.file_chromosomes[read.sam_refID]\n							  without checking that the value is positive first.',
	'sam_rname':			   ' [3rd column in SAM] The actual chromosome/contig name for the\n							  alignment. Will return "*" if refID is negative.',
	'bam_refID':			   ' The original bytes before decoding to sam_refID.',
	'sam_pos1':				' [4th column in SAM] The 1-based position of the alignment. Note\n							  that in SAM format values less than 1 are converted to "0" for\n							  "no data" and sam_pos1 will also do this.',
	'sam_pos0':				' The 0-based position of the alignment. Note that in SAM all\n							  positions are 1-based, but in BAM they are stored as 0-based.\n							  Unlike sam_pos1, negative values are kept as negative values,\n							  essentially giving one the decoded value as it was stored.',
	'bam_pos':				 ' The original bytes before decoding to sam_pos*.',
	'sam_mapq':				' [5th column in SAM] The Mapping Quality of the current alignment.',
	'bam_mapq':				' The original bytes before decoding to sam_mapq.',
	'sam_cigar_string':		' [6th column in SAM] The CIGAR string, as per the SAM format.\n							  Allowed values are "MIDNSHP=X".',
	'sam_cigar_list':		  ' A list of tuples with 2 values per tuple:\n							  the number of bases, and the CIGAR operation applied to those\n							  bases. Faster to calculate than sam_cigar_string.',
	'bam_cigar':			   ' The original bytes before decoding to sam_cigar_*.',
	'sam_next_refID':		  ' The sam_refID of the alignment\'s mate (if any). Note that as per\n							  sam_refID, this value can be negative and is not the actual\n							  chromosome name (see sam_pnext1).',
	'sam_rnext':			   ' [7th column in SAM] The chromosome name of the alignment\'s mate.\n							  Value is "*" if unmapped. Note that in a SAM file this value\n							  is "=" if it is the same as the sam_rname, however pybam will\n							  only do this if the user prints the whole SAM entry with "sam".',
	'bam_next_refID':		  ' The original bytes before decoding to sam_next_refID.',
	'sam_pnext1':			  ' [8th column in SAM] The 1-based position of the alignment\'s mate.\n							  Note that in SAM format values less than 1 are converted to "0"\n							  for "no data", and sam_pnext1 will also do this.',
	'sam_pnext0':			  ' The 0-based position of the alignment\'s mate. Note that in SAM all\n							  positions are 1-based, but in BAM they are stored as 0-based.\n							  Unlike sam_pnext1, negative values are kept as negative values\n							  here, essentially giving you the value as it was stored in BAM.',
	'bam_pnext':			   ' The original bytes before decoding to sam_pnext0.',
	'sam_tlen':				' [9th column in SAM] The TLEN value.',
	'bam_tlen':				' The original bytes before decoding to sam_tlen.',
	'sam_seq':				 ' [10th column in SAM] The SEQ value (DNA sequence of the alignment).\n							  Allowed values are "ACGTMRSVWYHKDBN and =".',
	'bam_seq':				 ' The original bytes before decoding to sam_seq.',
	'sam_qual':				' [11th column in SAM] The QUAL value (quality scores per DNA base\n							  in SEQ) of the alignment.',
	'bam_qual':				' The original bytes before decoding to sam_qual.',
	'sam_tags_list':		   ' A list of tuples with 3 values per tuple: a two-letter TAG ID, the\n							  type code used to describe the data in the TAG value (see SAM spec.\n							  for details), and the value of the TAG. Note that the BAM format\n							  has type codes like "c" for a number in the range -127 to +127,\n							  and "C" for a number in the range of 0 to 255.\n							  In a SAM file however, all numerical codes appear to just be stored\n							  using "i", which is a number in the range -2147483647 to +2147483647.\n							  sam_tags_list will therefore return the code used in the BAM file,\n							  and not "i" for all numbers.',
	'sam_tags_string':		 ' [12th column a SAM] Returns the TAGs in the same format as would be found \n							  in a SAM file (with all numbers having a signed 32bit code of "i").',
	'bam_tags':				' The original bytes before decoding to sam_tags_*.',
	'sam_bin':				 ' The bin value of the alignment (used for indexing reads).\n							  Please refer to section 5.3 of the SAM spec for how this\n							  value is calculated.',
	'bam_bin':				 ' The original bytes before decoding to sam_bin.',
	'sam_block_size':		  ' The number of bytes the current alignment takes up in the BAM\n							  file minus the four bytes used to store the block_size value\n							  itself. Essentially sam_block_size +4 == bytes needed to store\n							  the current alignment.',
	'bam_block_size':		  ' The original bytes before decoding to sam_block_size.',
	'sam_l_read_name':		 ' The length of the QNAME plus 1 because the QNAME is terminated\n							  with a NUL byte.',
	'bam_l_read_name':		 ' The original bytes before decoding to sam_l_read_name.',
	'sam_l_seq':			   ' The number of bases in the seq. Useful if you just want to know\n							  how many bases are in the SEQ but do not need to know what those\n							  bases are (which requires more decoding effort).',
	'bam_l_seq':			   ' The original bytes before decoding to sam_l_seq.',
	'sam_n_cigar_op':		  ' The number of CIGAR operations in the CIGAR field. Useful if one\n							  wants to know how many CIGAR operations there are, but does not\n							  need to know what they are.',
	'bam_n_cigar_op':		  ' The original bytes before decoding to sam_n_cigar_op.',
	'file_alignments_read':	' A running counter of the number of alignments ("reads"),\n							  processed thus far. Note the BAM format does not store\n							  how many reads are in a file, so the usefulness of this\n							  metric is somewhat limited unless one already knows how\n							  many reads are in the file.',
	'file_binary_header':	  ' From the first byte in the file, until the first byte of\n							  the first read. The original binary header.',
	'file_bytes_read':		 ' A running counter of the bytes read from the file. Note\n							  that as data is read in arbitary chunks, this is literally\n							  the amount of data read from the file/pipe by pybam.',
	'file_chromosome_lengths': ' The binary header of the BAM file includes chromosome names\n							  and chromosome lengths. This is a dictionary of chromosome-name\n							  keys and chromosome-length values.',
	'file_chromosomes':		' A list of chromosomes from the binary header.',
	'file_decompressor':	   ' BAM files are compressed with bgzip. The value here reflects\n							  the decompressor used. "internal" if pybam\'s internal\n							  decompressor is being used, "gzip" or "pigz" if the system\n							  has these binaries installed and pybam can find them.\n							  Any other value reflects a custom decompression command.',
	'file_directory':		  ' The directory the input BAM file can be found in. This will be\n							  correct if the input file is specified via a string or python\n							  file object, however if the input is a pipe such as sys.stdin, \n							  then the current working directory will be used.',
	'file_header':			 ' The ASCII portion of the BAM header. This is the typical header\n							  users of samtools will be familiar with.',
	'file_name':			   ' The file name (base name) of input file if input is a string or\n							  python file object. If input is via stdin this will be "<stdin>"'
}

class read():
	'''
	[ Dynamic Parser Example ]
	for alignment in pybam.read('/my/data.bam'):
		print alignment.sam_seq

	[ Static Parser Example ]
	for seq,mapq in pybam.read('/my/data.bam',['sam_seq','sam_mapq']):
		print seq
		print mapq

	[ Mixed Parser Example ]
	my_bam = pybam.read('/my/data.bam',['sam_seq','sam_mapq'])
	print my_bam._static_parser_code
	for seq,mapq in my_bam:
		if seq.startswith('ACGT') and mapq > 10:
			print my_bam.sam

	[ Custom Decompressor (from file path) Example ]
	my_bam = pybam.read('/my/data.bam.lzma',decompressor='lzma --decompress --stdout /my/data.bam.lzma')

	[ Custom Decompressor (from file object) Example ]
	my_bam = pybam.read(sys.stdin,decompressor='lzma --decompress --stdout') # data given to lzma via stdin
	
	[ Force Internal bgzip Decompressor ]
	my_bam = pybam.read('/my/data.bam',decompressor='internal')

	"print pybam.wat" in the python terminal to see the possible parsable values,
	or visit http://github.com/JohnLonginotto/pybam for the latest info.
	'''

	def __init__(self,f,fields=False, decompressor="internal"):
		
		self.file_bytes_read		 = 0
		self.file_chromosomes		= []
		self.file_alignments_read	= 0
		self.file_chromosome_lengths = {}

		if fields is not False:
			print(fields)
			if type(fields) is not list or len(fields) == 0:
				raise PybamError('\n\nFields for the static parser must be provided as a non-empty list. You gave a ' + str(type(fields)) + '\n')
			else:
				for field in fields:
					if field.startswith('sam') or field.startswith('bam'):
						if field not in list(parse_codes.keys()):
							raise PybamError('\n\nStatic parser field "' + str(field) + '" from fields ' + str(fields) + ' is not known to this version of pybam!\nPrint "pybam.wat" to see available field names with explinations.\n')
					else:
						raise PybamError('\n\nStatic parser field "' + str(field) + '" from fields ' + str(fields) + ' does not start with "sam" or "bam" and thus is not an avaliable field for the static parsing.\nPrint "pybam.wat" in interactive python to see available field names with explinations.\n')

		if decompressor:
			if type(decompressor) is str:
				 if decompressor != 'internal' and '{}' not in decompressor: raise PybamError('\n\nWhen a custom decompressor is used and the input file is a string, the decompressor string must contain at least one occurence of "{}" to be substituted with a filepath by pybam.\n')
			else: raise PybamError('\n\nUser-supplied decompressor must be a string that when run on the command line decompresses a named file (or stdin), to stdout:\ne.g. "lzma --decompress --stdout {}" if pybam is provided a path as input file, where {} is substituted for that path.\nor just "lzma --decompress --stdout" if pybam is provided a file object instead of a file path, as data from that file object will be piped via stdin to the decompression program.\n')

		## First we make a generator that will return chunks of uncompressed data, regardless of how we choose to decompress:
		def generator():
			DEVNULL = open(os.devnull, 'wb')
			# First we need to figure out what sort of file we have - whether it's gzip compressed, uncompressed, or something else entirely!
			if type(f) is str:
				try: 
					self._file = open(f,'rb')
				except: 
					raise PybamError('\n\nCould not open "' + str(self._file.name) + '" for reading!\n')
				try: 
					magic = self._file.read(4)
				except: 
					raise PybamError('\n\nCould not read from "' + str(self._file.name) + '"!\n')
			elif type(f) is file:
				self._file = f
				try: 
					magic = self._file.read(4)
				except: 
					raise PybamError('\n\nCould not read from "' + str(self._file.name) + '"!\n')
			else: 
				raise PybamError('\n\nInput file was not a string or a file object. It was: "' + str(f) + '"\n')

			self.file_name = os.path.basename(os.path.realpath(self._file.name))
			self.file_directory = os.path.dirname(os.path.realpath(self._file.name))
			
			if magic == b'BAM\1':
				# The user has passed us already unzipped BAM data! Job done :)
				data = b'BAM\1' + self._file.read(35536)
				self.file_bytes_read += len(data)
				self.file_decompressor = 'None'
				while data:
					yield data
					data = self._file.read(35536)
					self.file_bytes_read += len(data)
				self._file.close()
				DEVNULL.close()
				return

			elif magic == b"\x1f\x8b\x08\x04":  # The user has passed us compressed gzip/bgzip data, which is typical for a BAM file

				if decompressor is not False and decompressor != 'internal':
					if type(f) is str: 
						self._subprocess = subprocess.Popen(									decompressor.replace('{}',f),	shell=True, stdout=subprocess.PIPE, stderr=DEVNULL)
					else:
						self._subprocess = subprocess.Popen('{ printf "'+magic+'"; cat; } | ' + decompressor, stdin=self._file, shell=True, stdout=subprocess.PIPE, stderr=DEVNULL)
					self.file_decompressor = decompressor
					data = self._subprocess.stdout.read(35536)
					self.file_bytes_read += len(data)
					while data:
						yield data
						data = self._subprocess.stdout.read(35536)
						self.file_bytes_read += len(data)
					self._file.close()
					DEVNULL.close()
					return

				# else look for pigz or gzip:
				else:				
					try:
						self._subprocess = subprocess.Popen(["pigz"],stdin=DEVNULL,stdout=DEVNULL,stderr=DEVNULL)
						if self._subprocess.returncode is None: self._subprocess.kill()
						use = 'pigz'
					except OSError:
						try:
							self._subprocess = subprocess.Popen(["gzip"],stdin=DEVNULL,stdout=DEVNULL,stderr=DEVNULL)
							if self._subprocess.returncode is None: self._subprocess.kill()
							use = 'gzip'
						except OSError: 
							use = 'internal'
					if use != 'internal' and decompressor != 'internal':
						if type(f) is str: self._subprocess = subprocess.Popen([								   use , '--decompress','--stdout',	   f		   ], stdout=subprocess.PIPE, stderr=DEVNULL)
						else:			  self._subprocess = subprocess.Popen('{ printf "'+magic+'"; cat; } | ' + use + ' --decompress  --stdout', stdin=f, shell=True, stdout=subprocess.PIPE, stderr=DEVNULL)
						time.sleep(1)
						if self._subprocess.poll() == None:
							data = self._subprocess.stdout.read(35536)
							self.file_decompressor = use
							self.file_bytes_read += len(data)
							while data:
								yield data
								data = self._subprocess.stdout.read(35536)
								self.file_bytes_read += len(data)
							self._file.close()
							DEVNULL.close()
							return

					# Python's gzip module can't read from a stream that doesn't support seek(), and the zlib module cannot read the bgzip format without a lot of help:
					self.file_decompressor = 'internal'
					raw_data = magic + self._file.read(65536)
					self.file_bytes_read = len(raw_data)
					internal_cache = []
					blocks_left_to_grab = 50
					bs = 0
					checkpoint = 0
					decompress = zlib.decompress
					while raw_data:
						if len(raw_data) - bs < 35536:
							raw_data = raw_data[bs:] + self._file.read(65536)
							self.file_bytes_read += len(raw_data) - bs
							bs = 0
						magic = raw_data[bs:bs+4]
						if not magic: break # a child's heart
						if magic != b"\x1f\x8b\x08\x04": raise PybamError('\n\nThe input file is not in a format I understand. First four bytes: ' + repr(magic) + '\n')
						try:
							more_bs = bs + unpack("<H", raw_data[bs+16:bs+18])[0] +1
							internal_cache.append(decompress(raw_data[bs+18:more_bs-8],-15))
							bs = more_bs
						except: ## zlib doesnt have a nice exception for when things go wrong. just "error"
							header_data = magic + raw_data[bs+4:bs+12]
							header_size = 12
							extra_len = unpack("<H", header_data[-2:])[0]
							while header_size-12 < extra_len:
								header_data += raw_data[bs+12:bs+16]
								subfield_id = header_data[-4:-2]
								subfield_len = unpack("<H", header_data[-2:])[0]
								subfield_data = raw_data[bs+16:bs+16+subfield_len]
								header_data += subfield_data
								header_size += subfield_len + 4
								if subfield_id == 'BC': block_size = unpack("<H", subfield_data)[0]
							raw_data = raw_data[bs+16+subfield_len:bs+16+subfield_len+block_size-extra_len-19]
							crc_data = raw_data[bs+16+subfield_len+block_size-extra_len-19:bs+16+subfield_len+block_size-extra_len-19+8] # I have left the numbers in verbose, because the above try is the optimised code.
							bs = bs+16+subfield_len+block_size-extra_len-19+8
							zipped_data = header_data + raw_data + crc_data
							internal_cache.append(decompress(zipped_data,47)) # 31 works the same as 47.

							# Although the following in the bgzip code from biopython, its not needed if you let zlib decompress the whole zipped_data, header and crc, because it checks anyway (in C land)
							# I've left the manual crc checks in for documentation purposes:
							'''
							expected_crc = crc_data[:4]
							expected_size = unpack("<I", crc_data[4:])[0]
							if len(unzipped_data) != expected_size: print 'ERROR: Failed to unpack due to a Type 1 CRC error. Could the BAM be corrupted?'; exit()
							crc = zlib.crc32(unzipped_data)
							if crc < 0: crc = pack("<i", crc)
							else:	   crc = pack("<I", crc)
							if expected_crc != crc: print 'ERROR: Failed to unpack due to a Type 2 CRC error. Could the BAM be corrupted?'; exit()
							'''
						blocks_left_to_grab -= 1
						if blocks_left_to_grab == 0:
							yield b''.join(internal_cache)
							internal_cache = []
							blocks_left_to_grab = 50
					self._file.close()
					DEVNULL.close()
					if internal_cache != b'':
						yield b''.join(internal_cache)
					return

			elif decompressor is not False and decompressor != 'internal':
				# It wouldn't be safe to just print to the shell four random bytes from the beginning of a file, so instead it's
				# written to a temp file and cat'd. The idea here being that we trust the decompressor string as it was written by 
				# someone with access to python, so it has system access anyway. The file/data, however, should not be trusted.
				magic_file = os.path.join(tempfile.mkdtemp(),'magic')
				with open(magic_file,'wb') as mf: 
					mf.write(magic)
				if type(f) is str: 
					self._subprocess = subprocess.Popen(decompressor.replace('{}',f),	shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				else:
					self._subprocess = subprocess.Popen('{ cat "'+magic_file+'"; cat; } | ' + decompressor, stdin=self._file, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				self.file_decompressor = decompressor
				data = self._subprocess.stdout.read(35536)
				self.file_bytes_read += len(data)
				while data:
					yield data
					data = self._subprocess.stdout.read(35536)
					self.file_bytes_read += len(data)
				self._file.close()
				DEVNULL.close()
				return
			else:
				raise PybamError('\n\nThe input file is not in a format I understand. First four bytes: ' + repr(magic) + '\n')

		## At this point, we know that whatever decompression method was used, a call to self._generator will return some uncompressed data.
		self._generator = generator()
		
		print(self._generator)
				
		## So lets parse the BAM header:
		header_cache = b''
				
		while len(header_cache) < 4:
			header_cache += next(self._generator)
		p_from = 0; p_to = 4
		if header_cache[p_from:p_to] != b'BAM\x01':
			raise PybamError('\n\nInput file ' + self.file_name + ' does not appear to be a BAM file.\n')
			
		## Parse the BAM header:
		p_from = p_to; p_to += 4
		length_of_header = unpack('<i',header_cache[p_from:p_to])[0]
		p_from = p_to; p_to += length_of_header
		while len(header_cache) < p_to: header_cache += str(next(self._generator))
		self.file_header = header_cache[p_from:p_to]
		p_from = p_to; p_to += 4
		while len(header_cache) < p_to: header_cache += str(next(self._generator))
		number_of_reference_sequences = unpack('<i',header_cache[p_from:p_to])[0]
		
		for _ in range(number_of_reference_sequences):
			p_from = p_to; p_to += 4
			while len(header_cache) < p_to: header_cache += next(self._generator)
			l_name = unpack('<l',header_cache[p_from:p_to])[0]
			p_from = p_to; p_to += l_name
			while len(header_cache) < p_to: header_cache += next(self._generator)
			self.file_chromosomes.append(header_cache[p_from:p_to -1].decode('ascii'))
			p_from = p_to; p_to += 4
			while len(header_cache) < p_to: header_cache += next(self._generator)
			self.file_chromosome_lengths[self.file_chromosomes[-1]] = unpack('<l',header_cache[p_from:p_to])[0]

		self.file_bytes_read = p_to
		self.file_binary_header = memoryview(header_cache[:p_to])
		header_cache = header_cache[p_to:]

		
		# A quick check to make sure the header of this BAM file makes sense:
		chromosomes_from_header = []
		for line in str(self.file_header).split('\\n'):
			if line.startswith('@SQ\\tSN:'):
				chromosomes_from_header.append(line.split('\\t')[1][3:])
		if chromosomes_from_header != self.file_chromosomes:
			raise PybamWarn('For some reason the BAM format stores the chromosome names in two locations,\n	   the ASCII text header we all know and love, viewable with samtools view -H, and another special binary header\n	   which is used to translate the chromosome refID (a number) into a chromosome RNAME when you do bam -> sam.\n\nThese two headers should always be the same, but apparently they are not:\nThe ASCII header looks like: ' + self.file_header + '\nWhile the binary header has the following chromosomes: ' + self.file_chromosomes + '\n')
				
		#HL - decoding header
		self.file_header = self.file_header.decode()
	
		
		## Variable parsing:
		def new_entry(header_cache):
			cache = header_cache # we keep a small cache of X bytes of decompressed BAM data, to smoothen out disk access.
			p = 0 # where the next alignment/entry starts in the cache
			while True:
				try:
					while len(cache) < p + 4: cache = cache[p:] + next(self._generator); p = 0 # Grab enough bytes to parse blocksize
					self.sam_block_size  = unpack('<i',cache[p:p+4])[0]
					self.file_alignments_read += 1
					while len(cache) < p + 4 + self.sam_block_size:
						cache = cache[p:] + next(self._generator); p = 0 # Grab enough bytes to parse entry
				except StopIteration: break
				self.bam = cache[p:p + 4 + self.sam_block_size]
				p = p + 4 + self.sam_block_size
				yield self
		self._new_entry = new_entry(header_cache)
		
		

		def compile_parser(self,fields):
			temp_code = ''
			end_of_qname = False
			end_of_cigar = False
			end_of_seq = False
			end_of_qual = False
			dependencies = set(fields)

			if 'bam' in fields:
				fields[fields.index('bam')] = 'self.bam'

			if 'sam_block_size' in fields:
				fields[fields.index('sam_block_size')] = 'self.sam_block_size'

			if 'sam'				  in dependencies:
				dependencies.update(['sam_qname','sam_flag','sam_rname','sam_pos1','sam_mapq','sam_cigar_string','bam_refID','bam_next_refID','sam_rnext','sam_pnext1','sam_tlen','sam_seq','sam_qual','sam_tags_string'])

			if 'sam_tags_string'	  in dependencies:
				dependencies.update(['sam_tags_list'])

			if 'sam_pos1' in dependencies:
				temp_code += "\n		sam_pos1 = (0 if sam_pos0 < 0 else sam_pos0 + 1)"
				dependencies.update(['sam_pos0'])

			if 'sam_pnext1' in dependencies:
				temp_code += "\n		sam_pnext1 = (0 if sam_pnext0 < 0 else sam_pnext0 + 1)"
				dependencies.update(['sam_pnext0'])

			if 'sam_qname' in dependencies or 'bam_qname' in dependencies:
				temp_code += "\n		_end_of_qname = 36 + sam_l_read_name"
				dependencies.update(['sam_l_read_name'])
				end_of_qname = True

			if 'sam_cigar_string' in dependencies or 'sam_cigar_list' in dependencies or 'bam_cigar' in dependencies:
				if end_of_qname:
					pass
				else:
					temp_code += "\n		_end_of_qname = 36 + sam_l_read_name"
				temp_code += "\n		_end_of_cigar = _end_of_qname + (4*sam_n_cigar_op)"
				dependencies.update(['sam_l_read_name','sam_n_cigar_op'])
				end_of_cigar = True

			if 'sam_seq' in dependencies or 'bam_seq' in dependencies:
				if end_of_cigar:
					pass
				elif end_of_qname:
					temp_code += "\n		_end_of_cigar = _end_of_qname + (4*sam_n_cigar_op)"
				else:
					temp_code += "\n		_end_of_cigar = 36 + sam_l_read_name + (4*sam_n_cigar_op)"
				temp_code += "\n		_end_of_seq = _end_of_cigar + (-((-sam_l_seq)//2))"
				dependencies.update(['sam_l_seq','sam_n_cigar_op','sam_l_read_name'])
				end_of_seq = True

			if 'sam_qual' in dependencies or 'bam_qual' in dependencies:
				if end_of_seq:
					pass
				elif end_of_cigar:
					temp_code += "\n		_end_of_seq   = _end_of_cigar + (-((-sam_l_seq)//2))"
				elif end_of_qname:
					temp_code += "\n		_end_of_seq   = _end_of_qname + (4*sam_n_cigar_op) + (-((-sam_l_seq)//2))"
				else:
					temp_code += "\n		_end_of_seq   = 36 + sam_l_read_name + (4*sam_n_cigar_op) + (-((-sam_l_seq)//2))"
				temp_code += "\n		_end_of_qual = _end_of_seq + sam_l_seq"
				dependencies.update(['sam_l_seq','sam_n_cigar_op','sam_l_read_name'])
				end_of_qual = True

			if 'sam_tags_list' in dependencies or 'bam_tags' in dependencies:
				if end_of_qual:
					pass
				elif end_of_seq:
					temp_code += "\n		_end_of_qual = _end_of_seq + sam_l_seq"
				elif end_of_cigar:
					temp_code += "\n		_end_of_qual = _end_of_cigar + (-((-sam_l_seq)//2)) + sam_l_seq"
				elif end_of_qname:
					temp_code += "\n		_end_of_qual = _end_of_qname + (4*sam_n_cigar_op) + (-((-sam_l_seq)//2)) + sam_l_seq"
				else:
					temp_code += "\n		_end_of_qual = 36 + sam_l_read_name + (4*sam_n_cigar_op) + (-((-sam_l_seq)//2)) + sam_l_seq"
				dependencies.update(['sam_l_seq','sam_n_cigar_op','sam_l_read_name'])

			if 'sam_rname'	   in dependencies:
				temp_code += "\n		sam_rname = '*' if sam_refID < 0 else self.file_chromosomes[sam_refID]"
				dependencies.update(['sam_refID'])

			if 'sam_rnext'	   in dependencies:
				temp_code += "\n		sam_rnext = '*' if sam_next_refID < 0 else self.file_chromosomes[sam_next_refID]"
				dependencies.update(['sam_next_refID'])

			## First we figure out what data from the static portion of the BAM entry we'll need:
			tmp = {}
			tmp['code'] = 'def parser(self):\n	from array import array\n	from struct import unpack\n	for _ in self._new_entry:'
			tmp['last_start'] = None
			tmp['name_list']  = []
			tmp['dtype_list'] = []
			def pack_up(name,dtype,length,end,tmp):
				if name in dependencies:
					if tmp['last_start'] is None:
						tmp['last_start'] = end - length
					tmp['name_list'].append(name)
					tmp['dtype_list'].append(dtype)
				elif tmp['last_start'] is not None:
					tmp['code'] += '\n		' + ', '.join(tmp['name_list']) + ' = unpack("<' + ''.join(tmp['dtype_list']) + '",self.bam[' + str(tmp['last_start']) + ':' + str(end-length) + '])'
					if len(tmp['dtype_list']) == 1:
						tmp['code'] += '[0]'
					tmp['last_start'] = None
					tmp['name_list']  = []
					tmp['dtype_list'] = []

			pack_up('sam_refID',	   'i',4, 8,tmp)
			pack_up('sam_pos0',		'i',4,12,tmp)
			pack_up('sam_l_read_name', 'B',1,13,tmp)
			pack_up('sam_mapq',		'B',1,14,tmp)
			pack_up('sam_bin',		 'H',2,16,tmp)
			pack_up('sam_n_cigar_op',  'H',2,18,tmp)
			pack_up('sam_flag',		'H',2,20,tmp)
			pack_up('sam_l_seq',	   'i',4,24,tmp)
			pack_up('sam_next_refID',  'i',4,28,tmp)
			pack_up('sam_pnext0',	  'i',4,32,tmp)
			pack_up('sam_tlen',		'i',4,36,tmp)
			pack_up( None,			None,0,36,tmp) # To add anything not yet added.
			code = tmp['code']
			del tmp

			code += temp_code

			# Fixed-length BAM data (where we just grab the bytes, we dont unpack) can, however, be grabbed individually.
			if 'bam_block_size'  in dependencies: code += "\n		bam_block_size   = self.bam[0			 : 4				]"
			if 'bam_refID'	   in dependencies: code += "\n		bam_refID		= self.bam[4			 : 8				]"
			if 'bam_pos'		 in dependencies: code += "\n		bam_pos		  = self.bam[8			 : 12			   ]"
			if 'bam_l_read_name' in dependencies: code += "\n		bam_l_read_name  = self.bam[12			: 13			   ]"
			if 'bam_mapq'		in dependencies: code += "\n		bam_mapq		 = self.bam[13			: 14			   ]"
			if 'bam_bin'		 in dependencies: code += "\n		bam_bin		  = self.bam[14			: 16			   ]"
			if 'bam_n_cigar_op'  in dependencies: code += "\n		bam_n_cigar_op   = self.bam[16			: 18			   ]"
			if 'bam_flag'		in dependencies: code += "\n		bam_flag		 = self.bam[18			: 20			   ]"
			if 'bam_l_seq'	   in dependencies: code += "\n		bam_l_seq		= self.bam[20			: 24			   ]"
			if 'bam_next_refID'  in dependencies: code += "\n		bam_next_refID   = self.bam[24			: 28			   ]"
			if 'bam_pnext'	   in dependencies: code += "\n		bam_pnext		= self.bam[28			: 32			   ]"
			if 'bam_tlen'		in dependencies: code += "\n		bam_tlen		 = self.bam[32			: 36			   ]"
			if 'bam_qname'	   in dependencies: code += "\n		bam_qname		= self.bam[36			: _end_of_qname	]"
			if 'bam_cigar'	   in dependencies: code += "\n		bam_cigar		= self.bam[_end_of_qname : _end_of_cigar	]"
			if 'bam_seq'		 in dependencies: code += "\n		bam_seq		  = self.bam[_end_of_cigar : _end_of_seq	  ]"
			if 'bam_qual'		in dependencies: code += "\n		bam_qual		 = self.bam[_end_of_seq   : _end_of_qual	 ]"
			if 'bam_tags'		in dependencies: code += "\n		bam_tags		 = self.bam[_end_of_qual  :				  ]"

			if 'sam_qname'	   in dependencies:
				if 'bam_qname'   in dependencies: code += "\n		sam_qname		= bam_qname[:-1]"
				else:							 code += "\n		sam_qname		= self.bam[36			: _end_of_qname -1 ]"

			if 'sam_cigar_list'  in dependencies:
				if 'bam_cigar'   in dependencies: code += "\n		sam_cigar_list = [( cig >> 4  , cigar_codes[cig & 0b1111]) for cig in array('I', bam_cigar) ]"
				else:							 code += "\n		sam_cigar_list = [( cig >> 4  , cigar_codes[cig & 0b1111]) for cig in array('I', self.bam[_end_of_qname : _end_of_cigar]) ]"

			if 'sam_cigar_string'in dependencies:
				if 'bam_cigar'   in dependencies: code += "\n		sam_cigar_string = ''.join([		str(cig >> 4) + cigar_codes[cig & 0b1111] for cig	 in array('I', bam_cigar)])"
				else:							 code += "\n		sam_cigar_string = ''.join([		str(cig >> 4) + cigar_codes[cig & 0b1111] for cig	 in array('I', self.bam[_end_of_qname : _end_of_cigar]) ])"

			if 'sam_seq'		 in dependencies:
				if 'bam_seq'	 in dependencies: code += "\n		sam_seq = ''.join( [ dna_codes[dna >> 4] + dna_codes[dna & 0b1111]   for dna	 in array('B', bam_seq)])[:sam_l_seq]"
				else:							 code += "\n		sam_seq = ''.join( [ dna_codes[dna >> 4] + dna_codes[dna & 0b1111]   for dna	 in array('B', self.bam[_end_of_cigar : _end_of_seq])])[:sam_l_seq]"

			if 'sam_qual'		in dependencies:
				if 'bam_qual'	in dependencies: code += "\n		sam_qual = b''.join( [				   chr(ord(quality) + 33)		for quality in			bam_qual ])"
				else:							 code += "\n		sam_qual = b''.join( [				   chr(ord(quality) + 33)		for quality in			self.bam[_end_of_seq	   : _end_of_qual  ]])"

			if 'sam_tags_list'		in dependencies:
				code += '''
		sam_tags_list = []
		offset = _end_of_qual
		while offset != len(self.bam):
			tag_name = self.bam[offset:offset+2]
			tag_type = self.bam[offset+2]
			if tag_type == 'Z':
				offset_end = self.bam.index('\\0',offset+3)+1
				tag_data = self.bam[offset+3:offset_end-1]
			elif tag_type in CtoPy:
				offset_end = offset+3+py4py[tag_type]
				tag_data = unpack(CtoPy[tag_type],self.bam[offset+3:offset_end])[0]
			elif tag_type == 'B':
				offset_end = offset+8+(unpack('<i',self.bam[offset+4:offset+8])[0]*py4py[self.bam[offset+3]])
				tag_data = array(self.bam[offset+3] , self.bam[offset+8:offset_end] )
			else:
				print 'PYBAM ERROR: I dont know how to parse BAM tags in this format: ',repr(tag_type)
				print '			 This is simply because I never saw this kind of tag during development.'
				print '			 If you could mail the following chunk of text to john at john.uk.com, i will fix this up for everyone :)'
				print repr(tag_type),repr(self.bam[offset+3:end])
				exit()
			sam_tags_list.append((tag_name,tag_type,tag_data))
			offset = offset_end'''

			if 'sam_tags_string'	  in dependencies:
				code += "\n		sam_tags_string = '\t'.join(A + ':' + ('i' if B in 'cCsSI' else B)  + ':' + ((C.typecode + ',' + ','.join(map(str,C))) if type(C)==array else str(C)) for A,B,C in self.sam_tags_list)"

			if 'sam'				  in dependencies:
				code += "\n		sam = sam_qname + '\t' + str(sam_flag) + '\t' + sam_rname + '\t' + str(sam_pos1) + '\t' + str(sam_mapq) + '\t' + ('*' if sam_cigar_string == '' else sam_cigar_string) + '\t' + ('=' if bam_refID == bam_next_refID else sam_rnext) + '\t' + str(sam_pnext1) + '\t' + str(sam_tlen) + '\t' + sam_seq + '\t' + sam_qual + '\t' + sam_tags_string"

			code += '\n		yield ' + ','.join([x for x in fields]) + '\n'

			self._static_parser_code = code # "code" is the static parser's code as a string (a function called "parser")
			exec_dict = {				   # This dictionary stores things the exec'd code needs to know about, and will store the compiled function after exec()
				'unpack':unpack,
				'array':array,
				'dna_codes':dna_codes,
				'CtoPy':CtoPy,
				'py4py':py4py,
				'cigar_codes':cigar_codes
			}
			exec(code, exec_dict)			# exec() compiles "code" to real code, creating the "parser" function and adding it to exec_dict['parser']
			return(exec_dict['parser'])

		if fields:
			static_parser = compile_parser(self,fields)(self)
			def next_read(): return next(self._new_entry)
		else:
			def next_read(): return next(self._new_entry)					
		self.next = next_read
		
	def __next__(self):
		return next(self._new_entry)

	def __iter__(self): return self
	def __str__(self):  return self.sam

	## Methods to pull out raw bam data from entry (so still in its binary encoding). This can be helpful in some scenarios.
	@property
	def bam_block_size(self):   return			   self.bam[						: 4						 ] 
	@property
	def bam_refID(self):		return			   self.bam[ 4					  : 8						 ] 
	@property
	def bam_pos(self):		  return			   self.bam[ 8					  : 12						]
	@property
	def bam_l_read_name(self):  return			   self.bam[ 12					 : 13						] 
	@property
	def bam_mapq(self):		 return			   self.bam[ 13					 : 14						] 
	@property
	def bam_bin(self):		  return			   self.bam[ 14					 : 16						] 
	@property
	def bam_n_cigar_op(self):   return			   self.bam[ 16					 : 18						] 
	@property
	def bam_flag(self):		 return			   self.bam[ 18					 : 20						] 
	@property
	def bam_l_seq(self):		return			   self.bam[ 20					 : 24						] 
	@property
	def bam_next_refID(self):   return			   self.bam[ 24					 : 28						] 
	@property
	def bam_pnext(self):		return			   self.bam[ 28					 : 32						] 
	@property
	def bam_tlen(self):		 return			   self.bam[ 32					 : 36						] 
	@property
	def bam_qname(self):		return			   self.bam[ 36					 : self._end_of_qname		] 
	@property
	def bam_cigar(self):		return			   self.bam[ self._end_of_qname	 : self._end_of_cigar		] 
	@property
	def bam_seq(self):		  return			   self.bam[ self._end_of_cigar	 : self._end_of_seq		  ] 
	@property
	def bam_qual(self):		 return			   self.bam[ self._end_of_seq	   : self._end_of_qual		 ] 
	@property
	def bam_tags(self):		 return			   self.bam[ self._end_of_qual	  :						   ] 

	@property
	def sam_refID(self):		return unpack( '<i', self.bam[ 4					  :  8						] )[0]
	@property
	def sam_pos0(self):		 return unpack( '<i', self.bam[ 8					  : 12						] )[0]
	@property
	def sam_l_read_name(self):  return unpack( '<B', self.bam[ 12					 : 13						] )[0]
	@property
	def sam_mapq(self):		 return unpack( '<B', self.bam[ 13					 : 14						] )[0]
	@property
	def sam_bin(self):		  return unpack( '<H', self.bam[ 14					 : 16						] )[0]
	@property
	def sam_n_cigar_op(self):   return unpack( '<H', self.bam[ 16					 : 18						] )[0]
	@property
	def sam_flag(self):		 return unpack( '<H', self.bam[ 18					 : 20						] )[0]
	@property
	def sam_l_seq(self):		return unpack( '<i', self.bam[ 20					 : 24						] )[0]
	@property
	def sam_next_refID(self):   return unpack( '<i', self.bam[ 24					 : 28						] )[0]
	@property
	def sam_pnext0(self):	   return unpack( '<i', self.bam[ 28					 : 32						] )[0]
	@property
	def sam_tlen(self):		 return unpack( '<i', self.bam[ 32					 : 36						] )[0]
	@property
	def sam_qname(self):		return			   self.bam[ 36					 : self._end_of_qname -1	 ].decode() # -1 to remove trailing NUL byte
	@property
	def sam_cigar_list(self):   return		  [		  (cig >> 4  , cigar_codes[cig & 0b1111] ) for cig	 in array('I', self.bam[self._end_of_qname	 : self._end_of_cigar ])]
	@property
	def sam_cigar_string(self): return ''.join( [	   str(cig >> 4) + cigar_codes[cig & 0b1111]   for cig	 in array('I', self.bam[self._end_of_qname	 : self._end_of_cigar ])])
	@property
	def sam_seq(self):		  return ''.join( [ dna_codes[dna >> 4] +   dna_codes[dna & 0b1111]   for dna	 in array('B', self.bam[self._end_of_cigar	 : self._end_of_seq   ])])[:self.sam_l_seq] # As DNA is 4 bits packed 2-per-byte, there might be a trailing '0000', so we can either
	@property
	def sam_qual(self):  
		return ''.join( [					  chr(quality + 33)	   for quality in			self.bam[self._end_of_seq	   : self._end_of_qual  ]])
	@property
	def sam_tags_list(self):
		result = []
		offset = self._end_of_qual
		while offset != len(self.bam):
			tag_name = self.bam[offset:offset+2].decode()
			tag_type = chr(self.bam[offset+2])			
			if tag_type == 'Z':
				offset_end = self.bam.index(b'\x00',offset+3)+1
				tag_data = self.bam[offset+3:offset_end-1].decode()
			elif tag_type in CtoPy:
				offset_end = offset+3+py4py[tag_type]
				tag_data = unpack(CtoPy[tag_type],self.bam[offset+3:offset_end])[0]
			elif tag_type == 'B':
				offset_end = offset+8+(unpack('<i',self.bam[offset+4:offset+8])[0]*py4py[self.bam[offset+3]])
				tag_data = array(self.bam[offset+3] , self.bam[offset+8:offset_end] )
			else:
				print('PYBAM ERROR: I dont know how to parse BAM tags in this format: ',repr(tag_type))
				print('			 This is simply because I never saw this kind of tag during development.')
				print('			 If you could mail the following chunk of text to john at john.uk.com, ill fix this up :)')
				print(repr(tag_type),repr(self.bam[offset+3:end]))
				exit()
			result.append((tag_name,tag_type,tag_data))
			offset = offset_end
		return result
	@property
	def sam_tags_string(self):
		return '\t'.join(A + ':' + ('i' if B in 'cCsSI' else B)  + ':' + ((C.typecode + ',' + ','.join(map(str,C))) if type(C)==array else str(C)) for A,B,C in self.sam_tags_list)	

	## BONUS methods - methods that mimic how samtools works.
	@property
	def sam_pos1(self):		 return  0  if self.sam_pos0 < 0 else self.sam_pos0 + 1
	@property
	def sam_pnext1(self):	   return  0  if self.sam_pnext0 < 0 else self.sam_pnext0 + 1
	@property
	def sam_rname(self):		return '*' if self.sam_refID	  < 0 else self.file_chromosomes[self.sam_refID	 ]
	@property
	def sam_rnext(self):		return '*' if self.sam_next_refID < 0 else self.file_chromosomes[self.sam_next_refID]
	@property
	def sam(self):		return (
			self.sam_qname													 + '\t' +
			str(self.sam_flag)												 + '\t' +
			self.sam_rname													 + '\t' +
			str(self.sam_pos1)												 + '\t' +
			str(self.sam_mapq)												 + '\t' +
			('*' if self.sam_cigar_string == '' else self.sam_cigar_string)	+ '\t' +
			('=' if self.bam_refID == self.bam_next_refID else self.sam_rnext) + '\t' +
			str(self.sam_pnext1)											   + '\t' +
			str(self.sam_tlen)												 + '\t' +
			self.sam_seq													   + '\t' + 
			self.sam_qual													  + '\t' +
			self.sam_tags_string
		)

	## Internal methods - methods used to calculate where variable-length blocks start/end
	@property
	def _end_of_qname(self):	 return self.sam_l_read_name   + 36						# fixed-length stuff at the beginning takes up 36 bytes.
	@property
	def _end_of_cigar(self):	 return self._end_of_qname	 + (4*self.sam_n_cigar_op)   # 4 bytes per n_cigar_op
	@property
	def _end_of_seq(self):	   return self._end_of_cigar	 + (-((-self.sam_l_seq)//2)) # {blurgh}
	@property
	def _end_of_qual(self):	  return self._end_of_seq	   + self.sam_l_seq			# qual has the same length as seq

	def __del__(self):
		if self._subprocess.returncode is None: self._subprocess.kill()
		self._file.close()

class PybamWarn(Exception): pass
class PybamError(Exception): pass

def parse_to_mags_identical(contig_file_name, out_file_name):
	
	contigs =  open(contig_file_name, 'r')
	
	mags = open(out_file_name, "w")
	
	for line in contigs:
		if line[0] == ">":
			#set to new contig. One final loop of starts ends counts is needed
			current_contig = line[1:].strip().split()[0]
			print(current_contig, current_contig, sep = "\t", file = mags)
		
	contigs.close()
	
	mags.close()
	return("")
	
def tables_in_sqlite_db(conn):

	tabs = conn.execute("SELECT name FROM sqlite_master WHERE type='table';")
	tables = [
		v[0] for v in tabs.fetchall()
		if v[0] != "sqlite_sequence"
	]

	return tables

def sqldb_creation(contigs, mags, sample_reads, map_format, database):
	""" Read information provided by user and creates SQLite3 database
	
	Arguments:
		contigs {str} -- Location of fasta file with contigs of interest.
		mags {str} -- Location of tab separated file with contigs and their corresponding mags.
		sample_reads {list} -- Location of one or more read mapping results file(s).
		format {str} -- Format of read mapping (blast or sam).
		database {str} -- Name (location) of database to create.
	"""
	
	# ===== Database and table creation =====
	# Create or open database
	conn = sqlite3.connect(database)
	cursor = conn.cursor()
	
	#Check if there are tables ahead of time
	tables = tables_in_sqlite_db(cursor)
	
	#Clean out the old DB to begin with; effectively reinitialize.
	for table in tables:
		cursor.execute('DROP TABLE IF EXISTS '+ table)

	# Create lookup table (always creates a new one)
	cursor.execute('DROP TABLE IF EXISTS lookup_table')
	cursor.execute('CREATE TABLE lookup_table \
		(mag_name TEXT, mag_id INTEGER, contig_name TEXT, contig_id INTEGER)')
	# Create sample_info, mag_info, mags_per_sample, and gene_info tables
	cursor.execute('DROP TABLE IF EXISTS mag_info')
	cursor.execute('DROP TABLE IF EXISTS sample_info')
	cursor.execute('DROP TABLE IF EXISTS mags_per_sample')
	cursor.execute('CREATE TABLE mag_info \
		(mag_id INTEGER, contig_id INTEGER, contig_len INTEGER)')
	cursor.execute('CREATE TABLE sample_info \
		(sample_name TEXT, sample_id TEXT, sample_number INTEGER)')
	cursor.execute('CREATE TABLE mags_per_sample \
		(sample_name TEXT, mag_name TEXT)')
	# ========

	# === Extract sample information and save in into DB ===
	# Rename samples provided to avoid illegal names on files
	sampleid_to_sample = {}
	samples_to_db = []
	sample_number = 1
	for sample_name in sample_reads:
		sample_id = "sample_" + str(sample_number)
		samples_to_db.append((sample_name, sample_id, sample_number))
		sampleid_to_sample[sample_id] = sample_name
		sample_number += 1
	# Enter information into table
	cursor.execute("begin")
	cursor.executemany('INSERT INTO sample_info VALUES(?, ?, ?)', samples_to_db)
	cursor.execute('CREATE UNIQUE INDEX sample_index ON sample_info (sample_name)')
	cursor.execute("commit")
	# ========

	# === Extract contig information and MAG correspondence. Save into DB. ===
	
	# Get contig sizes
	contig_sizes = read_contigs(contigs)
	
	# Get contig - MAG information
	contig_mag_corresp = get_mags(mags)
	
	# Initialize variables
	contig_identifiers = []
	mag_ids = {}
	mag_id = 0
	contig_id = 1
	# The dictionary contig_information is important for speed when filling tables
	contig_information = {}
	# Iterate through contig - MAG pairs
	for contig_name, mag_name in contig_mag_corresp.items():
		# Store MAG (and contig) names and ids
		if mag_name in mag_ids:
			contig_identifiers.append((mag_name, mag_ids[mag_name], contig_name, contig_id))
			contig_information[contig_name] = [mag_name, mag_ids[mag_name], contig_id]
			contig_id += 1
		else:
			mag_id += 1
			mag_ids[mag_name] = mag_id
			contig_identifiers.append((mag_name, mag_ids[mag_name], contig_name, contig_id))
			contig_information[contig_name] = [mag_name, mag_ids[mag_name], contig_id]
			contig_id += 1
	cursor.executemany('INSERT INTO lookup_table VALUES(?, ?, ?, ?)', contig_identifiers)
	cursor.execute('CREATE INDEX mag_name_index ON lookup_table (mag_name)')
	conn.commit()
	# ========

	# === Fill contig length table ===
	contig_lengths = []
	for contig, contig_len in contig_sizes.items():
		# Get mag_id and contig_id
		sql_command = 'SELECT mag_id, contig_id from lookup_table WHERE contig_name = ?'
		cursor.execute(sql_command, (contig,))
		mag_contig_id = cursor.fetchone()		
		contig_lengths.append((mag_contig_id[0], mag_contig_id[1], contig_len))
	
	
	cursor.executemany('INSERT INTO mag_info VALUES(?, ?, ?)', contig_lengths)
	cursor.execute('CREATE INDEX mag_id_index ON mag_info (mag_id)')
	conn.commit()
	# ========

	# === Create one table with information per sample ===
	for sample_name in sampleid_to_sample.keys():
		# Drop if they exist
		cursor.execute('DROP TABLE IF EXISTS ' + sample_name)
		# Create tables once again
		cursor.execute('CREATE TABLE ' + sample_name + \
			' (mag_id INTEGER, contig_id INTEGER, identity FLOAT, start INTEGER, stop INTEGER)')
	# === Retrieve information from read mapping and store it ===
	# Read read mapping file for each sample and fill corresponding table
	for sample_name, mapping_file in sampleid_to_sample.items():
		mags_in_sample = []
		contigs_in_sample = save_reads_mapped(mapping_file, sample_name, map_format, cursor, conn)
		cursor.execute('SELECT contig_name, mag_name, mag_id FROM lookup_table')
		all_contigs = cursor.fetchall()
		for element in all_contigs:
			if element[0] in contigs_in_sample:
				if element[1] not in mags_in_sample:
					mags_in_sample.append(element[1])
				else:
					continue
			else:
				continue
		mags_in_sample = [(mapping_file, x) for x in mags_in_sample]
		cursor.executemany('INSERT INTO mags_per_sample VALUES(?, ?)', mags_in_sample)
		conn.commit()
	cursor.close()
	conn.commit()
	conn.close()

#This is now written with bamnostic so that windows supports bam format, too.
def save_reads_mapped(mapping_file, sample_name, map_format, cursor, conn):
	""" This script reads a read mapping file, extracts the contig to which each read maps,
		the percent id, the start and stop, and stores it in a table per sample.
	
	Arguments:
		mapping_file {str} -- Location of read mapping file.
		sample_name {str} -- Name of database sample name (form sample_#)
		map_format {str} -- Format of read mapping results (blast or sam)
		cursor {obj} -- Cursor to execute db instructions
		conn {obj} -- Connection handler to db.
	"""
	assert (map_format == "sam" or map_format == "bam" or map_format == "blast"), "Mapping format not recognized. Must be one of 'sam' 'bam' or 'blast'"
	contig_mag_corresp = {}
	contigs_in_sample = []
	# Retrieve mag information as mag_name, mag_id, contig_name, contig_id
	sql_command = 'SELECT * from lookup_table'
	cursor.execute(sql_command)
	contig_correspondence = cursor.fetchall()
	for contig_mag in contig_correspondence:
		contig_mag_corresp[contig_mag[2]] = [contig_mag[0], contig_mag[1], contig_mag[3]]
	
	# Read mapping files and fill sample tables
	if map_format == "blast":
		with open(mapping_file) as input_reads:
			record_counter = 0
			records = []
			for line in input_reads:
				# Commit changes after 500000 records
				if record_counter == 500000:
					cursor.execute("begin")
					cursor.executemany('INSERT INTO ' + sample_name + ' VALUES(?, ?, ?, ?, ?)', records)
					cursor.execute("commit")
					record_counter = 0
					records = []
				if line.startswith("#"):
					pass
				else:
					segment = line.split("\t")
					contig_ref = segment[1]
					# Exclude reads not associated with MAGs of interest
					if contig_ref not in contig_mag_corresp:
						continue
					else:
						if contig_ref not in contigs_in_sample:
							contigs_in_sample.append(contig_ref)
						pct_id = float(segment[2])
						pos1 = int(segment[8])
						pos2 = int(segment[9])
						start = min(pos1, pos2)
						end = start+(max(pos1, pos2)-min(pos1, pos2))
						mag_id = contig_mag_corresp[contig_ref][1]
						contig_id = contig_mag_corresp[contig_ref][2]
						records.append((mag_id, contig_id, pct_id, start, end))
						record_counter += 1
			# Commit remaining records
			if record_counter > 0:
				cursor.execute("begin")
				cursor.executemany('INSERT INTO ' + sample_name + ' VALUES(?, ?, ?, ?, ?)', records)
				cursor.execute("commit")
			# Create index for faster access
			cursor.execute('CREATE INDEX ' + sample_name + '_index on ' + sample_name + ' (mag_id)')
			
	if map_format == "sam":
		record_counter = 0
		records = []
		with open(mapping_file) as input_reads:
			for line in input_reads:
				if record_counter == 500000:
					cursor.execute("begin")
					cursor.executemany('INSERT INTO ' + sample_name + ' VALUES(?, ?, ?, ?, ?)', records)
					cursor.execute("commit")
					record_counter = 0
					records = []
				if "MD:Z:" not in line:
					continue
				else :
					segment = line.split()
					contig_ref = segment[2]
					# Exclude reads not associated with MAGs of interest
					if contig_ref not in contig_mag_corresp:
						continue
					else:
						if contig_ref not in contigs_in_sample:
							contigs_in_sample.append(contig_ref)
						# Often the MD:Z: field will be the last one in a magicblast output, but not always.
						# Therefore, start from the end and work in.
						iter = len(segment)-1
						mdz_seg = segment[iter]
						# If it's not the correct field, proceed until it is.
						while not mdz_seg.startswith("MD:Z:"):
							iter -= 1
							mdz_seg = segment[iter]
						#Remove the MD:Z: flag from the start
						mdz_seg = mdz_seg[5:]
						match_count = re.findall('[0-9]+', mdz_seg)
						sum=0
						for num in match_count:
							sum+=int(num)
						total_count = len(''.join([i for i in mdz_seg if not i.isdigit()])) + sum
						pct_id = (sum/(total_count))*100
						start = int(segment[3])
						end = start+total_count-1
						# Get mag_id and contig_id
						mag_id = contig_mag_corresp[contig_ref][1]
						contig_id = contig_mag_corresp[contig_ref][2]
						records.append((mag_id, contig_id, pct_id, start, end))
						record_counter += 1
			# Commit remaining records
			if record_counter > 0:
				cursor.execute("begin")
				cursor.executemany('INSERT INTO ' + sample_name + ' VALUES(?, ?, ?, ?, ?)', records)
				cursor.execute("commit")
			# Create index for faster access
			cursor.execute('CREATE INDEX ' + sample_name + '_index on ' + sample_name + ' (mag_id)')
	
	if map_format == "bam" and get_sys() == "Windows":
		record_counter = 0
		records = []
		
		#This iterator has a set of builtin functions that are called to access pos, ref name, MD:Z:
		for entry in read(mapping_file):
			#This line could allow processing to work like in SAM fmt, but is slower.
			#line = entry.to_string()
			if record_counter == 500000:
				cursor.execute("begin")
				cursor.executemany('INSERT INTO ' + sample_name + ' VALUES(?, ?, ?, ?, ?)', records)
				cursor.execute("commit")
				record_counter = 0
				records = []
			#has_tag returns true if the entry has a %ID relevant field
			mdz_seg = ""
			for a,b,c in entry.sam_tags_list:
				if a == "MD":
					mdz_seg = c
			
			if mdz_seg == "":
				continue
			else :
				#The individual read has a reference ID, and the file has a list of names via IDs.
				#The entry.ref_ID gets the ID number, and the .get_ref returns the actual name from the number
				contig_ref = entry.sam_rname
				#print(contig_ref)
				
				# Exclude reads not associated with MAGs of interest
				if contig_ref not in contig_mag_corresp:
					continue
				else:
					if contig_ref not in contigs_in_sample:
						contigs_in_sample.append(contig_ref)
					
					#Returns the MD:Z: segment
					match_count = re.findall('[0-9]+', mdz_seg)
					sum=0
					for num in match_count:
						sum+=int(num)
					total_count = len(''.join([i for i in mdz_seg if not i.isdigit()])) + sum
					pct_id = (sum/(total_count))*100
					
					
					#BAM files, unlike SAM files, are zero indexed. This +1 adjustment ensures SAM/BAM/R consistency
					start = entry.sam_pos1
					
					end = start+total_count-1
					# Get mag_id and contig_id
					mag_id = contig_mag_corresp[contig_ref][1]
					contig_id = contig_mag_corresp[contig_ref][2]
					
					
					records.append((mag_id, contig_id, pct_id, start, end))
					
					
					record_counter += 1
		# Commit remaining records
		if record_counter > 0:
			cursor.execute("begin")
			cursor.executemany('INSERT INTO ' + sample_name + ' VALUES(?, ?, ?, ?, ?)', records)
			cursor.execute("commit")
		# Create index for faster access
		cursor.execute('CREATE INDEX ' + sample_name + '_index on ' + sample_name + ' (mag_id)')
	
	if map_format == "bam" and get_sys() != "Windows":
		record_counter = 0
		records = []
		
		#This reader has some odd properties - most of it exists as a C interface
		#As a result, the entries are NOT the individual lines of the file and cannot be accessed as such
		#Instead, the iterator 'entry' returns a pointer to a location in memory based on the file
		#This iterator has a set of builtin functions that are called to access pos, ref name, MD:Z:
		input_reads = pysam.AlignmentFile(mapping_file, "rb")
		for entry in input_reads:
			#This line could allow processing to work like in SAM fmt, but is slower.
			#line = entry.to_string()
			if record_counter == 500000:
				cursor.execute("begin")
				cursor.executemany('INSERT INTO ' + sample_name + ' VALUES(?, ?, ?, ?, ?)', records)
				cursor.execute("commit")
				record_counter = 0
				records = []
			#has_tag returns true if the entry has a %ID relevant field
			if not entry.has_tag("MD"):
				continue
			else :
				#No longer needed because of pysam accesses
				#segment = line.split()
				
				#The individual read has a reference ID, and the file has a list of names via IDs.
				#The entry.ref_ID gets the ID number, and the .get_ref returns the actual name from the number
				contig_ref = input_reads.get_reference_name(entry.reference_id)
				
								
				# Exclude reads not associated with MAGs of interest
				if contig_ref not in contig_mag_corresp:
					continue
				else:
					if contig_ref not in contigs_in_sample:
						contigs_in_sample.append(contig_ref)
					
					#Returns the MD:Z: segment
					mdz_seg = entry.get_tag("MD")
					match_count = re.findall('[0-9]+', mdz_seg)
					sum=0
					for num in match_count:
						sum+=int(num)
					total_count = len(''.join([i for i in mdz_seg if not i.isdigit()])) + sum
					pct_id = (sum/(total_count))*100
					
					
					#BAM files, unlike SAM files, are zero indexed. This +1 adjustment ensures SAM/BAM/R consistency
					start = entry.reference_start+1
					
					
					end = start+total_count-1
					# Get mag_id and contig_id
					mag_id = contig_mag_corresp[contig_ref][1]
					contig_id = contig_mag_corresp[contig_ref][2]
					
					
					records.append((mag_id, contig_id, pct_id, start, end))
					
					
					record_counter += 1
		# Commit remaining records
		if record_counter > 0:
			cursor.execute("begin")
			cursor.executemany('INSERT INTO ' + sample_name + ' VALUES(?, ?, ?, ?, ?)', records)
			cursor.execute("commit")
		# Create index for faster access
		cursor.execute('CREATE INDEX ' + sample_name + '_index on ' + sample_name + ' (mag_id)')
	
	conn.commit()
	return contigs_in_sample

def add_sample(database, new_mapping_files, map_format):
	contig_mag_corresp = {}
	samples_dict = {}
	last_sample = 0
	conn = sqlite3.connect(database)
	cursor = conn.cursor()
	# Retrieve all sample information
	sql_command = 'SELECT * from sample_info'
	cursor.execute(sql_command)
	sample_information = cursor.fetchall()
	for sample in sample_information:
		samples_dict[sample[0]] = sample[1]
		if sample[2] > last_sample:
			last_sample = sample[2]
	
	# Retrieve contig - mag correspondence
	sql_command = 'SELECT * from lookup_table'
	cursor.execute(sql_command)
	contig_correspondence = cursor.fetchall()
	for contig_mag in contig_correspondence:
		contig_mag_corresp[contig_mag[2]] = [contig_mag[0], contig_mag[1], contig_mag[3]]
	for new_sample in new_mapping_files:
		# Check if new sample exists
		mags_in_sample = []
		if new_sample in samples_dict:
			sample_name = samples_dict[new_sample]
			# If it does, drop the reads that table from that sample and re-build it
			cursor.execute('DROP TABLE IF EXISTS ' + sample_name)
			cursor.execute('CREATE TABLE ' + sample_name + \
				' (mag_id INTEGER, contig_id INTEGER, identity FLOAT, start INTEGER, stop INTEGER)')
			cursor.execute('DELETE FROM mags_per_sample WHERE sample_name = ?', (new_sample,))
			conn.commit()
			contigs_in_sample = save_reads_mapped(new_sample, sample_name, map_format, cursor, conn)
			cursor.execute('SELECT contig_name, mag_name, mag_id FROM lookup_table')
			all_contigs = cursor.fetchall()
			for element in all_contigs:
				if element[0] in contigs_in_sample:
					if element[1] not in mags_in_sample:
						mags_in_sample.append(element[1])
					else:
						continue
				else:
					continue
			mags_in_sample = [(new_sample, x) for x in mags_in_sample]
			cursor.execute("begin")
			cursor.executemany('INSERT INTO mags_per_sample VALUES(?, ?)', mags_in_sample)
			cursor.execute("commit")

		else:
			# Otherwise create the new table and add the read information
			sample_name = "sample_" + str(last_sample + 1)
			last_sample += 1
			cursor.execute('CREATE TABLE ' + sample_name + \
				' (mag_id INTEGER, contig_id INTEGER, identity FLOAT, start INTEGER, stop INTEGER)')
			contigs_in_sample = save_reads_mapped(new_sample, sample_name, map_format, cursor, conn)
			cursor.execute('SELECT contig_name, mag_name, mag_id FROM lookup_table')
			all_contigs = cursor.fetchall()
			for element in all_contigs:
				if element[0] in contigs_in_sample:
					if element[1] not in mags_in_sample:
						mags_in_sample.append(element[1])
					else:
						continue
				else:
					continue
			mags_in_sample = [(new_sample, x) for x in mags_in_sample]
			cursor.execute("begin")
			cursor.executemany('INSERT INTO mags_per_sample VALUES(?, ?)', mags_in_sample)
			cursor.execute("commit")
			# Add to sample_info table
			new_record = (new_sample, sample_name, last_sample)
			cursor.execute('INSERT INTO sample_info VALUES(?, ?, ?)', new_record)
			conn.commit()
		conn.commit()
	conn.close()

def parse_prodigal_genes(prodigal_gff):
	gene_info = []
	with open(prodigal_gff, 'r') as prodigal_genes:
		for line in prodigal_genes:
			if line.startswith("#"):
				continue
			else:
				line = line.strip().split()
				contig = line[0]
				start = min(int(line[3]), int(line[4]))
				end = max(int(line[3]), int(line[4]))
				strand = line[6]
				annotation = line[8]
				gene_id = annotation.split(";")[0].split("_")[1]
				gene_id = contig + "_" + gene_id
				gene_info.append((contig, gene_id, start, end, strand, annotation))
	return gene_info

def add_gene_information(database, gene_info):
	conn = sqlite3.connect(database)
	cursor = conn.cursor()
	cursor.execute('DROP TABLE IF EXISTS gene_info')
	cursor.execute('CREATE TABLE gene_info \
		(contig_name TEXT, gene_name TEXT, gene_start INTEGER, gene_stop INTEGER, strand TEXT, annotation TEXT)')
	cursor.execute("begin")
	cursor.executemany('INSERT INTO gene_info VALUES(?, ?, ?, ?, ?, ?)', gene_info)
	cursor.execute("commit")

def add_gene_annotation(database, annotation):
	annotations = []
	conn = sqlite3.connect(database)
	cursor = conn.cursor()
	cursor.execute('DROP TABLE IF EXISTS gene_annotation')
	cursor.execute('CREATE TABLE gene_annotation \
		(gene_name TEXT, annotation TEXT)')
	with open(annotation) as gene_annot:
		for line in gene_annot:
			line = line.strip().split("\t")
			annotations.append((line[0], line[1]))
	cursor.execute("begin")
	cursor.executemany('INSERT INTO gene_annotation VALUES(?, ?)', annotations)
	cursor.execute("commit")

#Todo - make multiple
def read_contigs(contig_file_name):
	""" Reads a FastA file and returns
		sequence ids and sizes
	
	Arguments:
		contig_file_name {[str]} -- FastA file location
	Returns:
		[dict] -- Dictionary with ids and sizes
	"""
	
	contig_sizes = {}
	contig_length = 0
	contigs =  open(contig_file_name, 'r')
	
	#The ensuing loop commits a contig to the contig lengths dict every time a new contig is observed, i.e. whenever a current sequence has terminated.
	#This works both for single and splitline multi fastas.
	
	#The first line is manually read in so that its count can be gathered before it is committed - this basically skips the first iteration of the loop.
	current_contig = contigs.readline()[1:].strip().split()[0]
	
	for line in contigs:
		if line[0] == ">":
			#Add the contig that had 
			contig_sizes[current_contig] = contig_length
			
			#set to new contig. One final loop of starts ends counts is needed
			current_contig = line[1:].strip().split()[0]
			contig_length = 0
		else :
			contig_length += len(line.strip())
	
	contigs.close()
	
	#The loop never gets to commit on the final iteration, so this statement adds the last contig.
	contig_sizes[current_contig] = contig_length
	
	return contig_sizes

def get_mags(mag_file):
	""" Reads a file with columns:
		Contig_name MAG_name
		and returns the corresponding MAG per contig
	
	Arguments:
		mag_file {[str]} -- MAG correspondence file location
	Returns:
		[dict] -- Dictionary with contigs and corresponding MAG
	"""
	mag_dict = {}
	with open(mag_file, 'r') as mags:
		for line in mags:
			mag_contig = line.split()
			mag_dict[mag_contig[0]] = mag_contig[1]
	return mag_dict

#The purpose of this function is to take an empty recplot matrix object associated with one MAG, query the database for the sample and MAG in question,
#And fill the database with the returned information.
def fill_matrices(database, mag_id, sample_name, matrices, id_breaks, truncation_behavior = "edges", truncation_degree = 75):
	"""[summary]
	
	Arguments:
		database {str} -- Name of database to use (location).
		mag_id {int} -- ID of mag of interest.
		sample_name {str} -- Name of sample to retrieve reads from.
		matrices {dict} -- Empty matrices to fill.
		id_breaks {list} -- List of identity percentages to include.
	
	Returns:
		matrix [dict] -- Dictionary with list of arrays of start and stop positions
						 and filled matrix to plot.
	"""
	# Retrieve sample_id from sample_name provided
	conn = sqlite3.connect(database)
	cursor = conn.cursor()
	sql_command = 'SELECT sample_id from sample_info WHERE sample_name = ?'
	cursor.execute(sql_command, (sample_name,))
	sample_id = cursor.fetchone()[0]
	# Retrieve all read information from mag_name and sample_name provided
	sql_command = 'SELECT * from ' + sample_id + ' WHERE mag_id = ?'
	cursor.execute(sql_command, (mag_id,))
	
	contig_maxes = {}
	for contig in matrices:
		contig_maxes[contig] = max(matrices[contig][1])
	
	#TODO: Consider keeping/removing this in final...
	#read_counter = 0
	
	#TODO: We shouldn't need to fetch all reads. We can iterate on the cursor without doing this.
	#read_information = cursor.fetchall()
	#read_information is (mag_id contig_id perc_id read_start read_stop)
	#for read_mapped in read_information:
	for read_mapped in cursor:
		if read_mapped[1] in matrices:
			#read_counter += 1
			contig_id = read_mapped[1]
			read_start = read_mapped[3]
			read_stop = read_mapped[4]
			
			#Reads in the first (default) 75 bp don't count
			if read_start < truncation_degree:
				read_start = truncation_degree
			if read_start > read_stop:
				continue
			#Reads in the last (default) 75 bp don't count
			if read_stop > contig_maxes[contig_id] - truncation_degree:
				read_stop = contig_maxes[contig_id] - truncation_degree
			if read_start > read_stop:
				continue
				
			
			
			read_len = read_stop - read_start + 1
			read_id_index = bisect.bisect_right(id_breaks, read_mapped[2]) - 1
			read_start_loc = bisect.bisect_left(matrices[contig_id][1], read_start)
			read_stop_loc = bisect.bisect_left(matrices[contig_id][1], read_stop)
			# If the read falls entirely on a bin add all bases to the bin
			
			if read_start_loc == read_stop_loc:
				matrices[contig_id][2][read_start_loc][read_id_index] += read_len
			# On the contrary split bases between two or more bins
			else:
				for j in range(read_start_loc, read_stop_loc + 1):
					overflow = read_stop - matrices[contig_id][1][j]
					if overflow > 0:
						matrices[contig_id][2][j][read_id_index] += (read_len - overflow)
						read_len = overflow
					else :
						matrices[contig_id][2][j][read_id_index] += read_len

	return matrices

#The purpose of this function is to prepare an empty recplot matrix from a set of contig names and lengths associated with one MAG
def prepare_matrices(database, mag_name, width, bin_height, id_lower):
	#Prep percent identity breaks - always starts at 100 and proceeds down by bin_height steps until it cannot do so again without passing id_lower
	# Prep percent identity breaks - always starts at 100 and proceeds 
	# down by bin_height steps until it cannot do so again without passing id_lower
	id_breaks = []
	current_break = 100
	while current_break > id_lower:
		id_breaks.append(current_break)
		current_break -= bin_height
	id_breaks = id_breaks[::-1]
	
	zeroes = []
	for i in id_breaks:
		zeroes.append(0)
	
	# Retrieve mag_id from provided mag_name
	conn = sqlite3.connect(database)
	cursor = conn.cursor()
	sql_command = 'SELECT mag_id from lookup_table WHERE mag_name = ?'
	cursor.execute(sql_command, (mag_name,))
	mag_id = cursor.fetchone()[0]
	# Retrieve all contigs from mag_name and their sizes
	sql_command = 'SELECT contig_id, contig_len from mag_info WHERE mag_id = ?'
	cursor.execute(sql_command, (mag_id,))
	contig_sizes = cursor.fetchall()
	# Create matrices for each contig in the mag_name provided
	matrix = {}
	#begin reading contigs and determining their lengths.
	
	#id_len is a list of contig name, contig_length
	for id_len in contig_sizes:
		
		starts = []
		ends = []
		pct_id_counts = []
		
		contig_length = id_len[1]
		
		num_bins = int(contig_length / width)
		if num_bins < 1:
			num_bins = 1
		
		bin_width = (contig_length / num_bins)-1

		cur_bin_start = 1
		
		for i in range(1, num_bins):
			starts.append(int(cur_bin_start))
			ends.append(int((cur_bin_start+bin_width)))
			pct_id_counts.append(zeroes[:])
			cur_bin_start+=bin_width+1
		
		#Append final bin; guarantees the final bin = contig length
		starts.append(int(cur_bin_start))
		ends.append(contig_length)
		pct_id_counts.append(zeroes[:])
		
		matrix[id_len[0]] = [starts, ends, pct_id_counts]
		

	return(mag_id, matrix, id_breaks)

def prepare_matrices_genes(database, mag_name, bin_height, id_lower):
	#Prep percent identity breaks - always starts at 100 and proceeds down by bin_height steps until it cannot do so again without passing id_lower
	# Prep percent identity breaks - always starts at 100 and proceeds 
	# down by bin_height steps until it cannot do so again without passing id_lower
	id_breaks = []
	current_break = 100
	while current_break > id_lower:
		id_breaks.append(current_break)
		current_break -= bin_height
	id_breaks = id_breaks[::-1]
	
	zeroes = []
	for i in id_breaks:
		zeroes.append(0)
	
	# Retrieve mag_id from provided mag_name
	conn = sqlite3.connect(database)
	cursor = conn.cursor()
	sql_command = 'SELECT mag_id from lookup_table WHERE mag_name = ?'
	cursor.execute(sql_command, (mag_name,))
	mag_id = cursor.fetchone()[0]
	# Retrieve all contigs from mag_name and their sizes
	sql_command = 'SELECT contig_id, contig_len from mag_info WHERE mag_id = ?'
	cursor.execute(sql_command, (mag_id,))
	contig_sizes = cursor.fetchall()		
	
	#Contig ID is used above, but contig names will be needed to translate from the genes table and the contig data
	contig_curs = get_contig_names(database, mag_name)
	contig_names = []
	
	for item in contig_curs:
		contig_names.append(item[0])
			
	#Retrieve all gene info - order listed
	sql_command = 'SELECT contig_name, gene_name, gene_start, gene_stop, strand, annotation from gene_info WHERE contig_name IN (SELECT contig_name from lookup_table WHERE mag_id = ?)'
	cursor.execute(sql_command, (mag_id,))
	gene_sizes = cursor.fetchall()
	
	gene_matrix = {}
	
	#Separates gene data into contigs to allow for access by contig during matrix creation.
	for item in gene_sizes:
		if item[0] not in gene_matrix:
			gene_matrix[item[0]] = [[item[1]], [item[2]], [item[3]], [item[4]], [item[5]]]
		else:
			gene_matrix[item[0]][0].append(item[1])
			gene_matrix[item[0]][1].append(item[2])
			gene_matrix[item[0]][2].append(item[3])
			gene_matrix[item[0]][3].append(item[4])
			gene_matrix[item[0]][4].append(item[5])
	
	#Final data structures initialization
	matrix = {}
	annotation_matrix = {}
	
	#The script has to be able to access names through IDs, which this dict accomplishes. Otherwise it only works for the first genome in the DB
	contig_assoc = {}
	for i in range(0, len(contig_sizes)):
		contig_assoc[contig_sizes[i][0]] = contig_names[i]
	
	
	#Since genes are separated, id_len can be iterated through to give a length of contig for each and all the assoc. gene data can be accessed by name
		
	#id_len is a list of contig name, contig_length
	for id_len in contig_sizes:
			
		#These are the starts and ends of each gene on this contig
		starts = gene_matrix[str(contig_assoc[id_len[0]])][1][:]
		ends = gene_matrix[str(contig_assoc[id_len[0]])][2][:]
		
		gene_names = gene_matrix[str(contig_assoc[id_len[0]])][0][:]
		strand = gene_matrix[str(contig_assoc[id_len[0]])][3][:]
		annotation = gene_matrix[str(contig_assoc[id_len[0]])][4][:]
		
		final_gene_names = []
		final_gene_strands = []
		final_gene_annots = []
		#I shove around the starts/ends of overlapping genes to make sure that they no longer do, but that also means that 
		#the start/stop used for plotting is not guaranteed to be the real start/stop. I save the true ones here.
		final_gene_starts = []
		final_gene_ends = []

		pct_id_counts = []
				
		contig_length = id_len[1]
		
		
		final_starts = [1]
		final_ends = []
		
		#We want to see the genes in the context of the contig, so we have to add bins for all intergenic regions
		#This loop starts with 1 and creates a bin end based on the next gene start, then adds the gene's start and end, then adds a new start based on the gene's end
		for i in range(0, len(starts)):
			#We have to pad out the annotation information with blanks for the intergenic bins - 
			#determining their locations is easiest done alongside the initial creation of the bins rather than awkwardly afterwards
			
			#blanks
			final_gene_names.append("N/A")
			final_gene_strands.append("N/A")
			final_gene_annots.append("N/A")
			final_gene_starts.append("N/A")
			final_gene_ends.append("N/A")
			#actual gene info
			final_gene_names.append(gene_names[i])
			final_gene_strands.append(strand[i])
			final_gene_annots.append(annotation[i])
			final_gene_starts.append(str(starts[i]))
			final_gene_ends.append(str(ends[i]))
			#starts and stops added here
			final_ends.append(starts[i]-1)
			final_starts.append(starts[i])
			final_ends.append(ends[i])
			final_starts.append(ends[i]+1)
		
		#Caps the ends with the end of the contig
		final_ends.append(contig_length)
		
		#Finishes out the annotaion list - the last bin is always declared as non-genic, and it gets removed if that's untrue.
		final_gene_names.append("N/A")
		final_gene_strands.append("N/A")
		final_gene_annots.append("N/A")
		final_gene_starts.append("N/A")
		final_gene_ends.append("N/A")
		
		#If a gene starts at 1, then there should be a 1, 0 pair in the first position of the final starts/ends and the bin is unnecessary. This removes those
		if final_ends[0] < 1:
			del final_starts[0]
			del final_ends[0]
			del final_gene_annots[0]
			del final_gene_names[0]
			del final_gene_strands[0]
			del final_gene_starts[0]
			del final_gene_ends[0]		
			
		#If a gene terminates at the end of the contig, then the final start would be contig length + 1 and the bin is unnecessary. This removes those
		if final_starts[len(final_starts)-1] > final_ends[len(final_starts)-1]:
			del final_starts[len(final_starts)-1]
			del final_ends[len(final_starts)-1]
			del final_gene_annots[len(final_starts)-1]
			del final_gene_names[len(final_starts)-1]
			del final_gene_strands[len(final_starts)-1]
			del final_gene_starts[len(final_starts)-1]
			del final_gene_ends[len(final_starts)-1]

			
		
		starts = []
		ends = []
		gene_names = []
		strand = []
		annotation = []
		annot_start = []
		annot_end = []
		#If genes overlap, two problems are created:
		# (1) Intergenic region will have end >= start
		# (2) The end of the first overlapping gene will be <= the start of the second
		#This removes the intergenic regions and shoves the start/ends back roughly equal distances when they overlap - marginally inaccurate, but necessary to maintain histogram behavior
		#Has to be added to new lists because I cannot modify length of final starts/final ends inside of loop
		for i in range(0, len(final_starts)):
			if final_ends[i] >= final_starts[i]:
				starts.append(final_starts[i])
				ends.append(final_ends[i])
				gene_names.append(final_gene_names[i])
				strand.append(final_gene_strands[i])
				annotation.append(final_gene_annots[i])
				annot_start.append(final_gene_starts[i])
				annot_end.append(final_gene_ends[i])				
			else:
				midpt = int((final_ends[i-1] + final_starts[i+1])/2)
				#last element added was final_ends[i-1], and this is what needs updated
				ends[len(ends)-1] = midpt
				#This will naturally be added next loop, which should always happen because there will always be 1 more index to iterate over
				final_starts[i+1] = midpt + 1
				
		
				

		#Now add zeroes.
		pct_id_counts = []
		for index in starts:
			pct_id_counts.append(zeroes[:])
		
		#Add them to the data item
		matrix[id_len[0]] = [starts, ends, pct_id_counts]
		annotation_matrix[id_len[0]] = [gene_names, annot_start, annot_end, strand, annotation]

		
	return(mag_id, matrix, id_breaks, annotation_matrix)
	
#This function orchestrates calls to prepare_matrices and fill_matrices. 
#Translations between R and python through reticulate are inefficient and may cause errors with data typing.
#Constraining the transfers to arguments passed from R and a return passed from python alleviates these issues.
def extract_MAG_for_R(database, sample, mag_name, width, bin_height, id_lower):
	mag_id, matrix, id_breaks = prepare_matrices(database, mag_name, width, bin_height, id_lower)
	
	matrix = fill_matrices(database, mag_id, sample, matrix, id_breaks)
	
	return(matrix, id_breaks)
	
def extract_genes_MAG_for_R(database, sample, mag_name, bin_height, id_lower):
	mag_id, matrix, id_breaks, gene_table = prepare_matrices_genes(database, mag_name, bin_height, id_lower)
		
	matrix = fill_matrices(database, mag_id, sample, matrix, id_breaks)
		
	return(matrix, id_breaks, gene_table)
	
#This function queries the database and returns the names of all of the samples present within it
def assess_samples(database):
	# Retrieve sample_id from sample_name provided
	conn = sqlite3.connect(database)
	cursor = conn.cursor()
	sql_command = 'SELECT sample_name from sample_info'
	cursor.execute(sql_command)
	
	samples = []
	
	for samp in cursor:
		
		samples.append(samp)
		
	cursor.close()
	
	return(samples)

#This function queries the database and returns the MAGs covered by a given sample.	
def assess_MAGs(database, sample):
	# Retrieve sample_id from sample_name provided
	conn = sqlite3.connect(database)
	cursor = conn.cursor()
	sql_command = 'SELECT mag_name from mags_per_sample WHERE sample_name = ?'
	cursor.execute(sql_command, (sample,))
	
	mags = []
	
	for mag in cursor:
		mags.append(mag)
		
	cursor.close()
	
	return(mags)
	
def get_contig_names(database, mag_name):
	conn = sqlite3.connect(database)
	cursor = conn.cursor()
	sql_command = 'SELECT contig_name from lookup_table WHERE mag_name = ?'
	cursor.execute(sql_command, (mag_name,))
	
	names = []
	
	for contig in cursor:
		names.append(contig)
		
	cursor.close()
	
	return(names)
	
def add_genes_to_db(database, genes_file, gene_format):
	if(gene_format == "prodigal"):
		gene_information = parse_prodigal_genes(genes_file)
	else:
		print("I don't do that yet")
	
	add_gene_information(database, gene_information)
#add_gene_annotation(database, annotation)

def check_presence_of_genes(database):
	conn = sqlite3.connect(database)
	cursor = conn.cursor()
	sql_command = "SELECT name FROM sqlite_master WHERE type='table'"
	cursor.execute(sql_command,)
	
	checker = False
	
	for name in cursor:
		if str(name) == "('gene_info',)":
			checker = True
	
	cursor.close()
	return(checker)
	
#A function to check if a user-supplied file matches the appropriate fmt we're expecting
def detect_file_format(file):
	detected_format = "none"
	
	#I don't know how this could happen, but good to have it in place?
	toomuch = 0
	
	#If the user supplies a mismatch file, then this will run. Otw, we should only run the requested check to avoid excess effort.
	isfasta = detect_fasta(file)
	isbam = detect_bam(file)
	issam = detect_sam(file)
	isblast = detect_blast(file)
	isdb = detect_is_db(file)
	isassoc = detect_is_assoc(file)
	isprodigalgff = detect_is_prodigal(file)
		
	if isfasta:
		detected_format = "fasta"
		toomuch += 1
	if isbam:
		detected_format = "bam"
		toomuch += 1
	if issam:
		detected_format = "sam"
		toomuch += 1
	if isblast:
		detected_format = "blast"
		toomuch += 1
	if(isdb):
		detected_format = "database"
		toomuch += 1
	if(isassoc):
		detected_format = "assoc"
		toomuch += 1
	if(isprodigalgff):
		detected_format = "genes"
		toomuch += 1
	if toomuch > 1:
		detected_format = "none"
		
	return(detected_format)
	
def detect_fasta(file):
	fh = open(file, "r")
	
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
	
	return(fmt_fine)
	
def detect_bam(file):
	fh = open(file, "rb")
	
	fmt_fine = True
	
	head = fh.read(4)
	
	if head != b'BAM\1' and head != b"\x1f\x8b\x08\x04":
		fmt_fine = False
	
	fh.close()
	
	return(fmt_fine)
	
def detect_blast(file):
	fh = open(file, "r")
	
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
	
	return(fmt_fine)
	
def detect_sam(file):
	fh = open(file, "r")
	
	fmt_fine = True
	
	try:
		line = fh.readline()
	except :
		fmt_fine = False
	else:
		if line.startswith("@"):
			pass
		else:
			fmt_fine = False
	
	fh.close()
	
	return(fmt_fine)
	
def detect_is_db(dbname):
	isdb = True
	try:
		dburi = 'file:{}?mode=rw'.format(pathname2url(dbname))
		conn = sqlite3.connect(dburi, uri=True)
		
		try:
			ff = tables_in_sqlite_db(conn)
		except:
			isdb = False
		else:
			pass
		
		conn.close()
	except :
		isdb = False
	return isdb

#todo write these, include in fmt checker
def detect_is_assoc(file):
	fh = open(file, "r")
	
	fmt_fine = True
	
	#Check first line
	try:
		line = fh.readline()
	except :
		fmt_fine = False
	else :
		segment = line.strip().split()
		if len(segment) != 2:
			fmt_fine = False
	
	#Check the next 29 lines, or all until EOF.
	for i in range(1, 30):
		try:
			line = fh.readline()
		except :
			fmt_fine = False
			break
		else:
			segment = line.strip().split()
			if len(segment) == 2:
				pass
			else:
				#Basically an EOF check - skip the fmt overwrite if the file ended.
				if len(line) == 0:
					break
				fmt_fine = False
				break
	
	fh.close()
	
	return(fmt_fine)
	
def detect_is_prodigal(file):
	fh = open(file, "r")
	
	fmt_fine = True
	
	try:
		line = fh.readline()
	except :
		fmt_fine = False
	else:
		if line.startswith("##gff-version"):
			pass
		else:
			fmt_fine = False
	
	fh.close()
	
	return(fmt_fine)
	
