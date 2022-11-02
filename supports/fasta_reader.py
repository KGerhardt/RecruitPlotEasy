from .agnostic_reader import agnostic_reader
import os
def read_fasta(file, is_mag = False):
	cur_seq = ""
	cur_prot = ""
	
	contents = {}
	deflines = {}
	mag_ids = {}
	
	mag_id = os.path.basename(file)
	while mag_id != os.path.splitext(mag_id)[0]:
		mag_id = os.path.splitext(mag_id)[0]
		
	fasta = agnostic_reader(file)
	for line in fasta:
		if line.startswith(">"):
			if len(cur_seq) > 0:
				contents[cur_prot] = cur_seq
				deflines[cur_prot] = defline
				if is_mag:
					mag_ids[cur_prot] = mag_id
				else:
					mag_ids[cur_prot] = cur_prot
				
			cur_seq = ""
			cur_prot = line.strip().split()[0][1:]
			defline = line.strip()[len(cur_prot)+1 :].strip()
			
		else:
			cur_seq += line.strip()
				
	fasta.close()
	
	#Final iter
	if len(cur_seq) > 0:
		contents[cur_prot] = cur_seq
		deflines[cur_prot] = defline
		if is_mag:
			mag_ids[cur_prot] = mag_id
		else:
			mag_ids[cur_prot] = cur_prot
		
	return contents, deflines, mag_ids
	