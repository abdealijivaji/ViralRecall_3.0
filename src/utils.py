import os
from pathlib import Path
from pyfaidx import Fasta
from multiprocessing import cpu_count


def check_directory_permissions(directory_path : Path) -> bool:
   if os.access(directory_path, os.R_OK | os.W_OK):
      return True
   else:
        return False

valid_bases = set('ATCGN')

def filt_fasta(phagesize : int, genome_file: Fasta) -> list[str]:
	filt_contig_list = []
	
	for contig in genome_file.keys():
		if len(genome_file[contig]) >= phagesize :
			filt_contig_list.append(contig)
	
	return filt_contig_list


	
def is_DNA(genome : Fasta) :
	seq_for_check = set(str(genome[0][:1000]).upper())
	if set(seq_for_check).issubset(valid_bases) :
		pass
	else :
		raise Exception

def mp_cpu(cpu : int | None) -> int :
	'''Returns number of parallel processes to spawn in batch mode'''
	if cpu == None :
		return cpu_count()
	else :
		return cpu