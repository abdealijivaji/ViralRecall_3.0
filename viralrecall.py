#!/usr/bin/env python
from operator import itemgetter
import sys, os, argparse, itertools
import pandas as pd
from collections import defaultdict, namedtuple
from pathlib import Path

import multiprocessing.pool as mp
from pyfaidx import Fasta
from src.proteins import *
from src.utils import *

#warnings.simplefilter('ignore', Bio.BiopythonDeprecationWarning)

valid_bases = set('ATCGN')

# orf_finder = pyrodigal_gv.ViralGeneFinder(meta=True)


def get_region(df) :
	''' To get coordinates of putative viral regions '''
	rollscore = df['rollscore']
	vreg_index = [indx for indx in rollscore.index if rollscore[indx] > 10 ]
	
	return vreg_index
		

def count_hits(hits, contigs):
	''' Right now does listing of markers for each contig '''
	query2hits = defaultdict(list)
	for i in contigs:
		for j in hits:
			if i == j.contig:
				marker_score = j.HMM_hit + ": " + str(j.bitscore)
				query2hits[i].append(marker_score)

	return query2hits
		
	
def sliding_window_mean(df : pd.DataFrame, window_size) -> pd.Series:
	"""
	Calculate a rolling mean of the 'bitscore' column in the DataFrame df based on protein start positions.
	To allow variable window sizes, the 'pstart' column is converted to datetime and used as the index for rolling calculations.
	Offset is the window size in seconds
	"""
	means = pd.Series(0.0, index=df.index, dtype=float)  
	for name, grp in df.groupby('contig'):
		df = grp.loc[:,['pstart','bitscore']] #.sort_values(col_B).reset_index(drop=True)
		df['pstart'] = pd.to_datetime(grp['pstart'], unit='s', origin='unix', utc = True)  # convert starts to datetime index
		rollscr = df.rolling(window=window_size, min_periods=3, center=True, on='pstart').mean()
		means = means.combine(rollscr['bitscore'], max )
	
	return means
	

def run_program(input : str, 
				out_base : str,
				database : str, 
				window: int, 
				phagesize: int, 
				minscore: int, 
				evalue: float) : # , minhit, cpus,  plotflag, redo, flanking, batch, summary_file, contiglevel
	
	''' Main function to run ViralRecall '''

	infile_name = os.path.basename(input)
	out_dir = os.path.dirname(out_base)
	print("Running viralrecall on "+ infile_name + " and output will be deposited in "+ out_dir)
	
	# Using pyfaidx to parse fasta and looping 
	genome_file = Fasta(input)
	
	if is_DNA == False:
		print("{} does not look like a valid DNA sequence. Please check the input file.".format(infile_name))

	filt_contig_list = filt_fasta(phagesize, genome_file)
	
	if not filt_contig_list:
		print("No contigs longer than the minimum phage size of {phagesize} bp were found in the input file.\n" \
		" Not proceeding with the genome. Change phagesize parameter to a smaller value if needed.")
		return
	if os.path.isdir(out_dir) == False:
		os.mkdir(out_dir)

	vir_summary = []
	Summary = namedtuple("Summary", ["file", "contig", "contig_length", "num_viral_region", "vstart", "vend", "vir_length", "num_prots", "num_viral_hits", "score"])
	proteins, description = predict_proteins(genome_file, filt_contig_list, out_base)
	desc_df = pd.DataFrame(description)

	hmm_dir = database
	gvog_hmm = os.path.join(hmm_dir, "gvog.complete.hmm")
	hmm_results = search_with_pyhmmer(proteins, gvog_hmm, out_base, evalue)
	
	hmmout = out_base + ".hmmout"
	hmm_df = pd.DataFrame(hmm_results) #namedtuples perfectly compatible with pandas dataframe
	hmm_df.to_csv(hmmout, index=False, sep= "\t")

	# keep only top hits for each protein
	best_hits = hmm_df.sort_values(by = ["query", "bitscore" ], ascending=False).groupby(['contig', 'query']).nth(0).reset_index(drop=True)
	
	df = pd.merge(desc_df, best_hits, on = ['contig', 'query'], how = 'left')
	df.fillna({'HMM_hit': "no_hit"}, inplace = True) # add column names to dictionary here to replace NA with no_hit for colums with strings
	df.fillna(float(0), inplace = True)

	# Now to calculate score on a rolling window
	# converting starts to timestamp index 
	
	offset = pd.offsets.Second(window)
	# Use offset as the window argument for rolling
	df['rollscore'] = sliding_window_mean(df, offset)
	df.to_csv(out_base + ".tsv", index=False, sep="\t")

	# To extract viral regions, we need to extract regions with scores > threshold (minscore)
	viral_indices = { i : [] for i in filt_contig_list }  # dictionary to hold indices of viral regions for each contig
	for name, grp in df.groupby('contig'):
		above_threshold = grp.loc[grp['rollscore'] > minscore, :]
		if above_threshold.dropna().empty == False:
			thresh_index = above_threshold.index
			grp_index = []
			for key, group in itertools.groupby(enumerate(thresh_index), key =lambda x: x[0] - x[1]):
				indices = [*map(itemgetter(1), group)]
				grp_index.append(indices)
			strt_idx = grp_index[0][0]  # start of the first group
			vstart = above_threshold['pstart'].astype('int64')[strt_idx]  # start position of the first group
			if len(grp_index) == 1:
				
				end_idx = grp_index[0][-1] 
				vend = above_threshold['pend'].astype('int64')[end_idx]  # end position of the first group
				
				if vend - vstart >= phagesize:
					viral_indices[name].append([strt_idx, end_idx])
				
			else:
			# If difference between groups is less than 5 proteins and less than window size (default 15kb)
			# We don't update vstart, 
				for i in range(1, len(grp_index)):
					nxt_strt_idx = grp_index[i][0]
					end_idx = grp_index[i-1][-1]
					prot_diff = nxt_strt_idx - end_idx
					nxt_vstart = above_threshold['pstart'].astype('int64')[nxt_strt_idx]
					vend = above_threshold['pend'].astype('int64')[end_idx]
					bp_diff = nxt_vstart - vend
					
					if prot_diff > 5 and bp_diff > window :
						unq_hit = grp['HMM_hit'].iloc[strt_idx:end_idx].unique()			
						if vend - vstart >= phagesize and len(unq_hit) > 3:
							viral_indices[name].append([strt_idx, end_idx])
						strt_idx = nxt_strt_idx
						vstart = nxt_vstart
		
	viral_df = pd.DataFrame()
	for key, value in viral_indices.items():
		if len(value) >= 1:
			for idx, coords in enumerate(value, start=1):
				vregion_df = df.iloc[coords[0]:coords[1]+1]
				vregion_df.insert(0, 'Viral_region_number', f"vregion_{idx}", allow_duplicates=True)
				viral_df = pd.concat([viral_df, vregion_df], ignore_index=True)
				num_hits= len(vregion_df.loc[vregion_df['HMM_hit'] != "no_hit"])
				vstart , vend = vregion_df['pstart'].astype('int64')[coords[0]] , vregion_df['pend'].astype('int64')[coords[1]]
				vreg_head : str = genome_file[key][vstart:vend].fancy_name # type: ignore
				vreg_seq : str = genome_file[key][vstart:vend].seq # type: ignore
				# print(type(vreg_head))
				vreg_file = os.path.join(out_dir, f"{key}_vregion_{idx}.fna")
				with open(vreg_file, "w") as out:
					out.write(">" + vreg_head + "\n")
					out.write(vreg_seq)
				vir_summary.append(Summary(file=infile_name,
							   contig=key,
							   contig_length= len(genome_file[key][:].seq), # type: ignore
							   num_viral_region= f"vregion_{idx}",
							   vstart= vstart,
							   vend= vend,
							   vir_length= len(vreg_seq),
							   num_prots= len(vregion_df),
							   num_viral_hits= num_hits ,
							   score= vregion_df["bitscore"].mean()))
		else:
			vir_summary.append(Summary(file=infile_name,
							   contig=key,
							   contig_length= len(genome_file[key][:].seq), # type: ignore
							   num_viral_region= "no viral regions detected",
							   vstart= "NA",
							   vend= "NA",
							   vir_length= "NA",
							   num_prots= "NA",
							   num_viral_hits= "NA",
							   score= str(0)))
	summ_file = out_base + "_summary.tsv"
	summ_df = pd.DataFrame(vir_summary)
	summ_df.to_csv(summ_file, index=False, sep= "\t")
		
	vannot = out_base + "_viralregions.annot.tsv"
	viral_df.to_csv(vannot, index=False, sep= "\t")
	print(f"{infile_name} finished")
	return


def main(argv=None):

	args_list = parse_args()

	# set up object names for input/output/database folders
	input =  args_list.input # "/home/abdeali/viralR_test_input/Chlamy_punui_contig.fna" #
	project = args_list.project # "/home/abdeali/viralR_test_output/Chlamy_punui" # 
	# database = args_parser.database
	window = int(args_list.window)*1000 # convert to bp
	phagesize = int(args_list.minsize)*1000
	minscore = int(args_list.minscore)
	# minhit = int(args_parser.minhit)
	evalue = float(args_list.evalue)
	# cpus = args_parser.cpus
	# plotflag = args_parser.figplot
	# redo = args_parser.redo
	# contiglevel = args_parser.contiglevel
	# flanking = args_parser.flanking
	# batch = args_parser.batch
	
	input = os.path.expanduser(input) 
	project = os.path.expanduser(project) 
	
	database = "/home/abdeali/GVOGs/"
	
	base_dir = os.path.dirname(__file__) # path of viralrecall.py file
	database = os.path.join(base_dir, database)
	
	project = project.rstrip("/")

	existence = os.path.exists(input)
	indir = os.path.isdir(input)

	if indir and existence:
		if os.path.isdir(project):
			pass
		else:
			os.mkdir(project)
		
		files = os.scandir(input)
		file_list = [i.name for i in files if i.name.endswith('.fna') or i.name.endswith('.fasta') or i.name.endswith('.fa')]
		arg_list : list[tuple]= []
		for i in file_list:
			# Remove suffix before creating directory
			dir_name =  Path(i).with_suffix('')
			new_project = os.path.join(project, dir_name)
			out_base = os.path.join(new_project, dir_name)
			newinput = os.path.join(input, i)
			arg_list.append((newinput, out_base, database, window, phagesize, minscore, evalue)) 
		
		with mp.Pool() as pool:
			pool.starmap(run_program, arg_list)
			

	elif existence and not indir:
		
		out_base = os.path.join(project, os.path.basename(project)) 
		run_program(input, out_base, database, window, phagesize, minscore, evalue) # , minhit, cpus, plotflag, redo, flanking, batch, summary_file, contiglevel
	
	else:
		print("Input is not a valid directory or file. Please check the input path.")
	

if __name__ == '__main__':
	status = main()
	sys.exit(status)

