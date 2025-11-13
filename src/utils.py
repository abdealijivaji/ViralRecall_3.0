import argparse
from pathlib import Path
from pyfaidx import Fasta

valid_bases = set('ATCGN')

def filt_fasta(phagesize : int, genome_file: Fasta) -> list:
	filt_contig_list = []
	
	for contig in genome_file.keys():
		if len(genome_file[contig]) >= phagesize :
			filt_contig_list.append(contig)
	
	return filt_contig_list


def parse_args(argv=None) :
	
	args_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="ViralRecall v. 3.0: A flexible command-line tool for predicting NCLDV-like regions in genomic data \nFrank O. Aylward, Virginia Tech Department of Biological Sciences <faylward at vt dot edu>", epilog='*******************************************************************\n\n*******************************************************************')
	args_parser.add_argument('-i', '--input', required=False, help='Input FASTA file (ending in .fna)')
	args_parser.add_argument('-p', '--project', required=False, help='project name for outputs')
	args_parser.add_argument('-w', '--window', required=False, default=int(15), help='sliding window size to use for detecting viral regions (default=15 kb)')
	args_parser.add_argument('-m', '--minsize', required=False, default=int(10), help='minimum length of viral regions to report, in kilobases (default=10 kb)')
	args_parser.add_argument('-s', '--minscore', required=False, default=int(10), help='minimum score of viral regions to report, with higher values indicating higher confidence (default=1)')
	args_parser.add_argument('-e', '--evalue', required=False, default=str(1e-10), help='e-value that is passed to HMMER3 for the VOG hmmsearch (default=1e-10)')

	# args_parser.add_argument('-g', '--minhit', required=False, default=int(4), help='minimum number of viral hits that each viral region must have to be reported (default=4)')
	# args_parser.add_argument('-fl', '--flanking', required=False, default=int(0), help='length of flanking regions upstream and downstream of the viral region to output in the final .fna files (default=0)')
	# args_parser.add_argument('-t', '--cpus', required=False, default=str(1), help='number of cpus to use for the HMMER3 search')
	# args_parser.add_argument('-b', '--batch', type=bool, default=False, const=True, nargs='?', help='Batch mode: implies the input is a folder of .fna files that each will be run iteratively')
	# args_parser.add_argument('-r', '--redo', type=bool, default=False, const=True, nargs='?', help='run without re-launching prodigal and HMMER3 (for quickly re-calculating outputs with different parameters if you have already run once)')
	# args_parser.add_argument('-c', '--contiglevel', type=bool, default=False, const=True, nargs='?', help='calculate contig/replicon level statistics instead of looking at viral regions (good for screening contigs)')
	# args_parser.add_argument('-f', '--figplot', type=bool, default=False, const=True, nargs='?', help='Specify this flag if you would like a plot of the viral-like regions with the output')
	args_parser.add_argument('-v', '--version', action='version', version='ViralRecall v. 2.1')
	args_parser = args_parser.parse_args()

	return args_parser
	
def is_DNA(genome : Fasta) -> bool :
	seq_for_check = set(str(genome[0][:1000]).upper())
	if set(seq_for_check).issubset(valid_bases) :
		return True
	else :
		return False