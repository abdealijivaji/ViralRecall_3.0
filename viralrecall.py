#!/usr/bin/env python
import sys, os, pandas, argparse, time, warnings
from collections import defaultdict, namedtuple
import Bio
from Bio import SeqIO
from pathlib import Path
import pyrodigal_gv
import pyhmmer
from pyhmmer.easel import Alphabet, SequenceFile
from pyhmmer.plan7 import HMMFile, Pipeline
import multiprocessing.pool

warnings.simplefilter('ignore', Bio.BiopythonDeprecationWarning)

def load_sequences(input) :
	''' Load nucleotide FASTA sequences '''
	genome_file : list = list(SeqIO.parse(input, "fasta"))
	if not genome_file:
		print(f"{input} does not appear to be in FASTA format!")
	else:
		return genome_file


valid_bases = set('ATCGN')
def is_DNA(sequence) -> None:
	''' Checks if input is protein or DNA '''
	sequence = sequence[0]
	sequence = sequence.upper()
	
	if not set(sequence).issubset(valid_bases):
		print("Input Sequence contains non-DNA letters. Are you sure input is DNA sequence?")

def filt_contigs(input, phagesize) -> list :
	''' Remove contigs less than a certain size to save on compute time '''
	seq_file = input
	contig_len = int(phagesize)
	filt_seqs = [record.id for record in seq_file if len(record.seq) > contig_len]	
	#print(filt_seqs)
	if len(filt_seqs) < 1:
		print("genome file contains no contigs larger than {} kb.\nModify minimum contig length by -m flag".format(int(contig_len/1000)))
	return filt_seqs

def predict_proteins(input, contigs, outfile, cpus):
	record = input
	contig_list = contigs
	threads = int(cpus)
	prot_out = outfile + ".faa"  
	filt_seqs = [record for record in record if record.id in contig_list]
	#print(filt_seqs)
	orf_finder = pyrodigal_gv.ViralGeneFinder(meta=True)
	#with open(outfile, "w") as fout :
	proteins = []
	header = []
	Header = namedtuple("Header", ["contig", "query", "pstart", "pend", "pstrand", "gen_code"])

	with open(prot_out, "w") as fout :
		for seqrecord in filt_seqs:
			genes = orf_finder.find_genes(bytes(seqrecord.seq))
			genes.write_translations(fout, sequence_id=seqrecord.id)
			for n, gene in enumerate(genes, start= 1):
				aa = gene.translate()
				prot_id = seqrecord.id + "_" + str(n)
				head = (Header(seqrecord.id, prot_id, gene.begin, gene.end, gene.strand, gene.translation_table)
							# f"{id} # {gene.begin} # {gene.end} # "
							# + f"{gene.strand} # ID={n}; "
							# + f"partial={int(gene.partial_begin)}{int(gene.partial_end)}; "
							# + f"start_type={gene.start_type} ; rbs_motif={gene.rbs_motif}; "
							# + f"rbs_spacer={gene.rbs_spacer}; "
							# + f"genetic_code={gene.translation_table}; "
							# + f"gc_cont={gene.gc_cont:.3f}"
						)
				header.append(head)
				proteins.append((prot_id, aa))
	#print(proteins)
	return proteins, prot_out, header
					


def search_with_pyhmmer(proteins, project, hmm_path):
	alphabet = Alphabet.amino()
	pipeline = Pipeline(alphabet)
	results = []
	Result = namedtuple("Result", ["contig", "query", "HMM_hit", "bitscore", "evalue"])	
	#seqs: list[pyhmmer.easel.DigitalSequence] = [ pyhmmer.easel.TextSequence(sequence=prot.value(), name=bytes(prot.key(), 'UTF-8')).digitize(alphabet) for prot in proteins ]
	seqs = []
	#outfile = project + "/" + project + ".hmmout"
	#digseqs = []
	for i, (prot_id, aa) in enumerate(proteins):
		seqs.append(pyhmmer.easel.TextSequence(name = bytes(prot_id,  'UTF-8'), sequence= aa).digitize(alphabet))
		#digseqs = pyhmmer.easel.DigitalSequence(alphabet, name = bytes(id,  'UTF-8'), sequence= bytes(aa, 'UTF-8'))
	#digseqs = seqs.digitize(alphabet)
	digseqs = pyhmmer.easel.DigitalSequenceBlock(alphabet, seqs)
	#print(pyhmmer.easel.DigitalSequence(alphabet, seqs))
	
	
	with HMMFile(hmm_path) as hmm_file:
		# Convert protein predictions into pyhmmer Sequence objects
		hits = pyhmmer.hmmer.hmmsearch(hmm_file, digseqs)

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
	# with open(input) as handle:
	# 	is_fasta(handle)
	# 	filt_contigs(handle, phagesize)
	
	genome_file = load_sequences(input)

	is_DNA(genome_file)
	
	filt_contig_list = filt_contigs(genome_file, phagesize)
	
	proteins, outfile, description = predict_proteins(genome_file, filt_contig_list, out_base, cpus)
	

	hmm_dir = database
	hmm_results = search_with_pyhmmer(proteins, out_base, hmm_dir)
	

	hmmout = out_base + ".hmmout"
	hmm_df = pandas.DataFrame(hmm_results) #namedtuples perfectly compatible with pandas dataframe fuck
	hmm_df.to_csv(hmmout, index=False, sep= "\t")

	


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
	input = "/home/abdeali/viralR_test_input/" # args_parser.input # 
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
	batch = True #  args_parser.batch
	
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

