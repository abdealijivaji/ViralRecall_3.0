#!/usr/bin/env python

import argparse
import pandas as pd
from collections import namedtuple
from pathlib import Path
import multiprocessing.pool as mp
# from pyfaidx import Fasta, FastaIndexingError
from src.proteins import *
from src.utils import *
from src.vreg_annot import *
import src.log as log

__version__ = 3.0


def parse_args(argv=None) :
	
	args_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="ViralRecall v. 3.0: A flexible command-line tool for predicting NCLDV-like regions in genomic data \nFrank O. Aylward, Virginia Tech Department of Biological Sciences <faylward at vt dot edu>", epilog='*******************************************************************\n\n*******************************************************************')
	args_parser.add_argument('-i', '--input', required=False, help='Input Fasta file or directory of fasta files (ending in .fna, .fasta, or .fa)')
	args_parser.add_argument('-o', '--outdir', required=False, help='Output directory name')
	args_parser.add_argument('-w', '--window', required=False, default=int(15), help='sliding window size to use for detecting viral regions (default=15 kb)')
	args_parser.add_argument('-m', '--minsize', required=False, default=int(10), help='minimum length of viral regions to report, in kilobases (default=10 kb)')
	args_parser.add_argument('-s', '--minscore', required=False, default=int(10), help='minimum score of viral regions to report, with higher values indicating higher confidence (default=10)')
	args_parser.add_argument('-e', '--evalue', required=False, default=str(1e-10), help='e-value that is passed to pyHmmer for hmmsearch (default=1e-10)')
	args_parser.add_argument('-g', '--minhit', required=False, default=int(4), help='minimum number of unique viral hits that each viral region must have to be reported (default=4)')
	# args_parser.add_argument('-fl', '--flanking', required=False, default=int(0), help='length of flanking regions upstream and downstream of the viral region to output in the final .fna files (default=0)')
	args_parser.add_argument('-c', '--cpus', required=False, default=None, help='number of cpus to use for the HMMER3 search')
	# args_parser.add_argument('-r', '--redo', type=bool, default=False, const=True, nargs='?', help='run without re-launching prodigal and HMMER3 (for quickly re-calculating outputs with different parameters if you have already run once)')
	# args_parser.add_argument('-c', '--contiglevel', type=bool, default=False, const=True, nargs='?', help='calculate contig/replicon level statistics instead of looking at viral regions (good for screening contigs)')
	# args_parser.add_argument('-f', '--figplot', type=bool, default=False, const=True, nargs='?', help='Specify this flag if you would like a plot of the viral-like regions with the output')
	args_parser.add_argument('-v', '--version', action='version', version='ViralRecall v. 3.0')
	args_parser = args_parser.parse_args()

	return args_parser
	

def run_program(input : Path, 
				out_base : Path,
				database : Path, 
				window: int, 
				phagesize: int, 
				minscore: int, 
				evalue: float,
				minhit: int) : # , minhit, cpus,  plotflag, redo, flanking, batch, summary_file, contiglevel
	
	''' Main function to run ViralRecall '''

	out_base.parent.mkdir(exist_ok=True)
	logger = log.setup_logger(out_base.parent)
	try :
		genome_file = load_genome(input)
	except Exception as e:
		logger.error(e)
	
	logger.info(f"Running viralrecall on {input.name}")
	

	filt_contig_list = filt_fasta(phagesize, genome_file)
	
	if not filt_contig_list:
		print(f"No contigs longer than the minimum phage size of {phagesize} bp were found in {input.name}.\n" \
		" Not proceeding with the genome. Change phagesize parameter to a smaller value if needed.")
		return

	vir_summary = []
	Summary = namedtuple("Summary", ["file", "contig", "contig_length", "num_viral_region", "vstart", "vend", "vir_length", "num_prots", "num_viral_hits", "score"])
	proteins, description = predict_proteins(genome_file, filt_contig_list, out_base)
	desc_df = pd.DataFrame(description)

	gvog_hmm = database / "gvog_mirus_cat.h3m"
	hits = search_with_pyhmmer(proteins, gvog_hmm, evalue)
	hmm_results = parse_hmmer(hits, out_base)
	
	hmmout = out_base.with_suffix(".hmmout")
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
	df = merge_annot(df)
	df.to_csv(out_base.with_suffix(".tsv"), index=False, sep="\t")

	# To extract viral regions, we need to extract regions with scores > threshold (minscore)
	
	viral_indices = extract_reg(window, phagesize, minscore, filt_contig_list, df, minhit)
		
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
				vreg_file = out_base.parent / f"{key}_vregion_{idx}.fna"
				with open(vreg_file, "w") as out:
					out.write(">" + vreg_head + "\n")
					out.write(vreg_seq)
				vir_summary.append(Summary(file=input.name,
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
			vir_summary.append(Summary(file=input.name,
							   contig=key,
							   contig_length= len(genome_file[key][:].seq), # type: ignore
							   num_viral_region= "no viral regions detected",
							   vstart= "NA",
							   vend= "NA",
							   vir_length= "NA",
							   num_prots= "NA",
							   num_viral_hits= "NA",
							   score= str(0)))
	summ_file = Path(str(out_base) + "_summary.tsv")
	summ_df = pd.DataFrame(vir_summary)
	summ_df.to_csv(summ_file, index=False, sep= "\t")
		
	vannot = Path(str(out_base) + "_viralregions.annot.tsv")
	viral_df.to_csv(vannot, index=False, sep= "\t")
	logger.info(f"Finished running Viralrecall on {input.name} and found {len(viral_indices)} viral regions" )
	return




def main(argv=None):

	args_list = parse_args()

	# set up object names for input/output/database folders
	input =   "~/viralR_test_input/Chlamy_punui_contig.fna" # args_list.input #
	project =  "~/viralR_test_output/Chlamy_punui" # args_list.project #
	# database = args_parser.database
	window = int(args_list.window)*1000 # convert to bp
	phagesize = int(args_list.minsize)*1000
	minscore = int(args_list.minscore)
	minhit = int(args_list.minhit)
	evalue = float(args_list.evalue)
	cpus = args_list.cpus
	# plotflag = args_parser.figplot
	# redo = args_parser.redo
	# contiglevel = args_parser.contiglevel
	# flanking = args_parser.flanking
	# batch = args_parser.batch
	
	input = Path(input).expanduser() 
	project = Path(project).expanduser() 
	
	 # path of viralrecall.py file
	database = Path(__file__).parent / "hmm"
	prep_hmm(database)
	# project = project.rstrip("/")

	existence = input.exists()
	indir = input.is_dir()

	cpus = mp_cpu(cpus)
	if indir and existence:
		if check_directory_permissions(input) == False :
			print(f"Insufficient read/write permissions for {input} \nQuitting")
			return
		
		project.mkdir(exist_ok=True)
		
		try :
			file_list = [i.name for i in input.iterdir() if i.name.endswith(('.fna', '.fasta', '.fa'))]
		except :
			raise FileNotFoundError(f"Input Directory : {input}, does not contain genome fasta files. Quitting")
		arg_list : list[tuple]= []
		for i in file_list:
			# Remove suffix before creating directory
			dir_name =  Path(i).stem
			new_project = project / dir_name
			out_base = new_project / dir_name
			newinput = input / i
			arg_list.append((newinput, out_base, database, window, phagesize, minscore, evalue, minhit)) 
		
		with mp.Pool(cpus) as pool:
			pool.starmap(run_program, arg_list)
			

	elif existence and not indir:
		
		out_base = project / input.stem
		run_program(input, out_base, database, window, phagesize, minscore, evalue, minhit) # , minhit, cpus, plotflag, redo, flanking, batch, summary_file, contiglevel
	
	else:
		print("Input is not a valid directory or file. Please check the input path.")
	

if __name__ == '__main__':
	main()
	

