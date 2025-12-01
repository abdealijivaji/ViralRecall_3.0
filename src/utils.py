import os
from pathlib import Path
from pyfaidx import Fasta, FastaIndexingError
from multiprocessing import cpu_count
from pyhmmer import hmmpress, plan7




def load_genome(input : Path) -> Fasta :
	try :
		genome_file = Fasta(input)
		is_DNA(genome_file)
	except FastaIndexingError:
		raise FastaIndexingError(f"Input file {input.name} is not in Fasta format. Please check input file")
	except ValueError:
		raise ValueError(f"{input.name} does not look like a valid DNA sequence. Please check input file")
	
	return genome_file

def compress_hmm(hmmpath : Path) -> None :
    outname = hmmpath.with_suffix('')
    
    with plan7.HMMFile(hmmpath) as hf :
        print(f"Compressing {hmmpath.name} into binary format")
        hmmpress(hf, outname)
    print(f"Compressed {hmmpath.name}")

def prep_hmm(hmm_dir : Path) -> None :
	for file in hmm_dir.glob('*.hmm') :
		if file.with_suffix('.h3m').exists() == False :
			compress_hmm(file)
	

def check_directory_permissions(directory_path : Path) -> bool:
   if os.access(directory_path, os.R_OK | os.W_OK):
      return True
   else:
        return False

valid_bases = set('ATCGN')

def filt_fasta(phagesize : int, genome_file: Fasta) -> list[str]:
	filt_contig_list : list[str] = []
	
	for contig in genome_file.keys():
		if len(genome_file[contig]) >= phagesize :
			filt_contig_list.append(contig)
	
	return filt_contig_list
	
def is_DNA(genome : Fasta) -> None :
	seq_for_check = set(str(genome[0][:1000]).upper())
	if not set(seq_for_check).issubset(valid_bases) :
		raise ValueError

def mp_cpu(cpu : int | None) -> int :
	'''Returns number of parallel processes to spawn in batch mode'''
	if cpu == None :
		return cpu_count()
	else :
		return cpu
	
	