#!/usr/bin/env python
import sys, os, argparse, time, warnings
import pandas as pd
from collections import defaultdict, namedtuple
import Bio
from Bio import SeqIO
from pathlib import Path
import pyrodigal_gv
from pyhmmer import easel, plan7, hmmer
import multiprocessing.pool as mp
from pyfaidx import Fasta

warnings.simplefilter('ignore', Bio.BiopythonDeprecationWarning)

valid_bases = set('ATCGN')

orf_finder = pyrodigal_gv.ViralGeneFinder(meta=True)

def predict_proteins(input, contigs, outfile):
	#threads = int(cpus)
	proteins = []
	header = []
	Header = namedtuple("Header", ["contig", "query", "pstart", "pend", "pstrand", "gen_code"])

	with open(outfile, "w") as fout :
		for seqrecord in contigs:
			sequence = bytes(str(input[seqrecord]), 'UTF-8')
			genes = orf_finder.find_genes(sequence)
			genes.write_translations(fout, sequence_id=seqrecord)
			for n, gene in enumerate(genes, start= 1):
				aa = gene.translate()
				prot_id = seqrecord + "_" + str(n)
				head = Header(seqrecord, prot_id, gene.begin, gene.end, gene.strand, gene.translation_table)
				header.append(head)
				proteins.append((prot_id, aa))
	#print(proteins)
	return proteins, header
					

alphabet = easel.Alphabet.amino()

def search_with_pyhmmer(proteins, hmm_path):
	results = []
	Result = namedtuple("Result", ["contig", "query", "HMM_hit", "bitscore", "evalue"])	
	seqs = []
	
	for i, (prot_id, aa) in enumerate(proteins):
		seqs.append(easel.TextSequence(name = bytes(prot_id,  'UTF-8'), sequence= aa).digitize(alphabet))
	
	digseqs = easel.DigitalSequenceBlock(alphabet, seqs)
	
	
	with plan7.HMMFile(hmm_path) as hmm_file:
		# Convert protein predictions into pyhmmer Sequence objects
		hits = hmmer.hmmsearch(hmm_file, digseqs)

		for hitlist in hits:
			for hit in hitlist:
				# hit.description = bytes(hit.name.decode().split(None, maxsplit=1)[1], 'UTF-8')
				Contig =  hit.name.decode().rsplit("_", maxsplit=1)[0]
				eval = "%.3g" % hit.evalue
				results.append(Result(
					Contig ,
					hit.name.decode() ,
					hitlist.query.name.decode(),
					round(hit.score, 2),
					eval,
					# hit.description.decode()
				))
	return results

def get_region(hits, description, contig) :
	''' To get coordinates of putative viral regions '''
	viral_contig_reg = defaultdict(list, { k:[] for k in contig})

	start = []
	end = []
	for j in description:
	
		for i in hits:
			if j.query == i.query :
				start.append(j.pstart)
				end.append(j.pend)
	
	
	
			

def count_hits(hits, contigs):
	''' Rught now does listing of markers for each contig '''
	query2hits = defaultdict(list)
	for i in contigs:
		for j in hits:
			if i == j.contig:
				marker_score = j.HMM_hit + ": " + str(j.bitscore)
				query2hits[i].append(marker_score)

	return query2hits
		
	

	

def run_program(input, out_base, database, window, phagesize, minscore, minhit, evalue, cpus, plotflag, redo, flanking, batch, summary_file, contiglevel):
	
	# genome_file = load_sequences(input)
	#is_DNA(genome_file)
	#filt_contig_list = filt_contigs(genome_file, phagesize)

	# Using pyfaidx to parse fasta and looping 
	genome_file = Fasta(input)
	filt_contig_list = []
	
	for contig in genome_file.keys():
		if len(genome_file[contig]) >= phagesize :
			filt_contig_list.append(contig)
	
	prot_out = out_base + ".faa"  

	proteins, description = predict_proteins(genome_file, filt_contig_list, prot_out)
	desc_df = pd.DataFrame(description)
	#print(desc_df)

	hmm_dir = database
	hmm_results = search_with_pyhmmer(proteins, hmm_dir)
	
	hmmout = out_base + ".hmmout"
	hmm_df = pd.DataFrame(hmm_results) #namedtuples perfectly compatible with pandas dataframe fuck
	#print(hmm_df)
	hmm_df.to_csv(hmmout, index=False, sep= "\t")
	df = pd.merge(desc_df, hmm_df, on = ['contig', 'query'], how = 'left')
	df.fillna({'HMM_hit': "no_hit"}, inplace = True) # add column names to dictionary here to replace NA with no_hit for colums with strings
	df.fillna(float(0), inplace = True)

	#print(df)
	#df.to_csv(out_base, sep = '\t', index = False)

	# Now to calculate score on a rolling window

	print(df)

	#get_region(hmm_results, description, filt_contig_list)
	contig_hits = count_hits(hmm_results, filt_contig_list)
	# print(hmm_df)
	# with open(hmmout, "w") as out:
	# 	for res in hmm_results:
	# 		line = [res.query, res.HMM_hit, str(res.bitscore), str(res.evalue), res.description]
	# 		out.writelines(str(line) + '\n')



def main(argv=None):

	args_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="ViralRecall v. 2.0: A flexible command-line tool for predicting NCLDV-like regions in genomic data \nFrank O. Aylward, Virginia Tech Department of Biological Sciences <faylward at vt dot edu>", epilog='*******************************************************************\n\n*******************************************************************')
	args_parser.add_argument('-i', '--input', required=False, help='Input FASTA file (ending in .fna)')
	args_parser.add_argument('-p', '--project', required=False, help='project name for outputs')
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
	input = "/home/abdeali/viralR_test_input/cat_genomes.fna" # args_parser.input # 
	project = "/home/abdeali/test_out/" #args_parser.project
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
	
	database = "/home/abdeali/packages/ViralRecall_3.0/hmm/merged_GVOGs.hmm"

	# path of viralrecall.py file
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

		file_list = os.listdir(input)
		for i in file_list:
			# Remove suffix before creating directory
			dir_name = Path(i)
			dir_name = dir_name.with_suffix('')
			new_project = os.path.join(project, dir_name)
			out_base = os.path.join(new_project, dir_name)
			if os.path.isdir(new_project):
				pass
			else:
				os.mkdir(new_project)
			#newproject = os.path.splitext(newproject)[0]
			newinput = os.path.join(input, i)
			print("Running viralrecall on "+ i + " and output will be deposited in "+ new_project)
			run_program(newinput, out_base, database, window, phagesize, minscore, minhit, evalue, cpus, plotflag, redo, flanking, batch, summary_file, contiglevel)
	else:
		#summary_file = 1
		# creates folder wherever you want now
		if os.path.isdir(project):
			pass
		else:
			os.mkdir(project)
		out_base = os.path.join(project, os.path.basename(project)) 
		summary_file = open(os.path.join(project, "batch_summary.txt"), "w")
		summary_file.write("genome\tcontigs_tested\n")

		run_program(input, out_base, database, window, phagesize, minscore, minhit, evalue, cpus, plotflag, redo, flanking, batch, summary_file, contiglevel)

	return 0

if __name__ == '__main__':
	status = main()
	sys.exit(status)

