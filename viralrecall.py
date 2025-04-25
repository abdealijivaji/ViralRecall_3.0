#!/usr/bin/env python
import sys, os, re, shlex, subprocess, pandas, numpy, itertools, argparse, time
from collections import defaultdict
from Bio import SeqIO
from operator import itemgetter
from itertools import islice
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
import math
from pathlib import Path
import pyrodigal_gv

def is_fasta(input):
	"""
	Checks if a file is in FASTA format.
	"""
	seqdict = SeqIO.to_dict(SeqIO.parse(input, "fasta"))
	print(input[1])
	if len(seqdict) < 1:
		#print(input +" does not appear to be in FASTA format! Quitting")
		quit()
# def check_fasta(input):
# 	file = list(SeqIO.parse(input, "fasta"))
	

def filt_contigs(input, phagesize) :
	seq_file = SeqIO.parse(input, "fasta")
	contig_len = int(phagesize)
	filt_seqs = (record for record in seq_file if len(record.seq) > contig_len)
	print(filt_seqs)
	return filt_seqs



def run_program(input, project, database, window, phagesize, minscore, minhit, evalue, cpus, plotflag, redo, flanking, batch, summary_file, contiglevel):
	with open(input) as handle:
		is_fasta(handle)
		is_fasta(filt_contigs(handle))
		


def main(argv=None):

	args_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="ViralRecall v. 2.0: A flexible command-line tool for predicting NCLDV-like regions in genomic data \nFrank O. Aylward, Virginia Tech Department of Biological Sciences <faylward at vt dot edu>", epilog='*******************************************************************\n\n*******************************************************************')
	args_parser.add_argument('-i', '--input', required=True, help='Input FASTA file (ending in .fna)')
	args_parser.add_argument('-p', '--project', required=True, help='project name for outputs')
	args_parser.add_argument('-db', '--database', required=False, default="GVOG", help='Viral HMM database to use. Options are "general" for the general VOG db, "GVOG" for the GVOG db, and "marker" for searching only a set of 10 conserved NCLDV markers (good for screening large datasets). See README for details')
	args_parser.add_argument('-w', '--window', required=False, default=int(15), help='sliding window size to use for detecting viral regions (default=15)')
	args_parser.add_argument('-m', '--minsize', required=False, default=int(10), help='minimum length of viral regions to report, in kilobases (default=10)')
	args_parser.add_argument('-s', '--minscore', required=False, default=int(1), help='minimum score of viral regions to report, with higher values indicating higher confidence (default=1)')
	args_parser.add_argument('-g', '--minhit', required=False, default=int(4), help='minimum number of viral hits that each viral region must have to be reported (default=4)')
	args_parser.add_argument('-e', '--evalue', required=False, default=str(1e-10), help='e-value that is passed to HMMER3 for the VOG hmmsearch (default=1e-10)')
	args_parser.add_argument('-fl', '--flanking', required=False, default=int(0), help='length of flanking regions upstream and downstream of the viral region to output in the final .fna files (default=0)')
	args_parser.add_argument('-t', '--cpus', required=False, default=str(1), help='number of cpus to use for the HMMER3 search')
	args_parser.add_argument('-b', '--batch', type=bool, default=False, const=True, nargs='?', help='Batch mode: implies the input is a folder of .fna files that each will be run iteratively')
	args_parser.add_argument('-r', '--redo', type=bool, default=False, const=True, nargs='?', help='run without re-launching prodigal and HMMER3 (for quickly re-calculating outputs with different parameters if you have already run once)')
	args_parser.add_argument('-c', '--contiglevel', type=bool, default=False, const=True, nargs='?', help='calculate contig/replicon level statistics instead of looking at viral regions (good for screening contigs)')
	args_parser.add_argument('-f', '--figplot', type=bool, default=False, const=True, nargs='?', help='Specify this flag if you would like a plot of the viral-like regions with the output')
	args_parser.add_argument('-v', '--version', action='version', version='ViralRecall v. 2.1')
	args_parser = args_parser.parse_args()

	# set up object names for input/output/database folders
	input = args_parser.input
	project = args_parser.project
	database = args_parser.database
	window = int(args_parser.window)
	phagesize = int(args_parser.minsize)*1000
	minscore = int(args_parser.minscore)
	minhit = int(args_parser.minhit)
	evalue = str(args_parser.evalue)
	cpus = args_parser.cpus
	plotflag = args_parser.figplot
	redo = args_parser.redo
	contiglevel = args_parser.contiglevel
	flanking = args_parser.flanking
	batch = args_parser.batch
	
	# Base dir to establish path of viralrecall.py file
	base_dir = Path(__file__).parent.resolve()
	#print(base_dir)
	project = project.rstrip("/")
	if batch:
		if os.path.isdir(project):
			pass
		else:
			os.mkdir(project)
		
		summary_file = open(os.path.join(project, "batch_summary.txt"), "w")
		summary_file.write("genome\tcontigs_tested\n")

		# if os.path.isdir(project):
		# 	pass
		# else:
		# 	os.mkdir(project)
			
		file_list = os.listdir(input)
		for i in file_list:
			# Remove suffix before creating directory
			dir_name = Path(i)
			dir_name = dir_name.with_suffix('')
			print(dir_name)
			newproject = os.path.join(project, dir_name)
			if os.path.isdir(newproject):
				pass
			else:
				os.mkdir(newproject)
			#newproject = os.path.splitext(newproject)[0]
			newinput = os.path.join(input, i)
			print("Running viralrecall on "+ i + " and output will be deposited in "+ newproject)
			#run_program(newinput, newproject, database, window, phagesize, minscore, minhit, evalue, cpus, plotflag, redo, flanking, batch, summary_file, contiglevel)
	else:
		#summary_file = 1
		# creates folder wherever you want now
		if os.path.isdir(project):
			pass
		else:
			os.mkdir(project)
		summary_file = open(os.path.join(project, "batch_summary.txt"), "w")
		summary_file.write("genome\tcontigs_tested\n")
		
		run_program(input, project, database, window, phagesize, minscore, minhit, evalue, cpus, plotflag, redo, flanking, batch, summary_file, contiglevel)

	return 0

if __name__ == '__main__':
	status = main()
	sys.exit(status)

