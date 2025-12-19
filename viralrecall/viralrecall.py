#!/usr/bin/env python

import argparse, time
import pandas as pd
from collections import namedtuple
from pathlib import Path
import multiprocessing.pool as mp
from .proteins import predict_proteins, search_with_pyhmmer, parse_hmmer
from .utils import load_genome, prep_hmm, filt_fasta, check_directory_permissions, mp_cpu, find_db
from .vreg_annot import sliding_window_mean, merge_annot, extract_reg, str_hits
from .log import setup_logger
from .plot_vreg import plot_vreg

__version__ = 3.0


def parse_args(argv=None) :
	
	args_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, prog='viralrecall', add_help=False,
									   description=f"Viralrecall: A command-line tool for predicting NCLDV-like regions in genomic data \nAbdeali (Ali) M. Jivaji, Virginia Tech Department of Biological Sciences <abdeali@vt.edu>", 
									   epilog='*******************************************************************\n\n*******************************************************************')
	required = args_parser.add_argument_group('required arguments')
	required.add_argument('-i', '--input', required=True, 
						  help='Input Fasta file or directory of fasta files (ending in .fna, .fasta, or .fa)')
	required.add_argument('-o', '--outdir', required=True, 
						  help='Output directory name')
	required.add_argument('-d', '--database', required=True, 
						  help='Database directory name, e.g. ~/hmm\n(download the database using viralrecall_database)')
	optional = args_parser.add_argument_group('optional arguments')
	optional.add_argument('-w', '--window', required=False, 
						  default=int(15), help='sliding window size to use for detecting viral regions (default=15 kb)')
	optional.add_argument('-m', '--minsize', required=False, 
						  default=int(10), help='minimum length of viral regions to report, in kilobases (default=10 kb)')
	optional.add_argument('-s', '--minscore', required=False, 
						  default=int(10), help='minimum score of viral regions to report, with higher values indicating higher confidence (default=10)')
	optional.add_argument('-e', '--evalue', required=False, 
						  default=str(1e-10), help='e-value that is passed to pyHmmer for hmmsearch (default=1e-10)')
	optional.add_argument('-g', '--minhit', required=False, 
						  default=int(4), help='minimum number of unique viral hits that each viral region must have to be reported (default=4)')
	optional.add_argument('-c', '--cpus', required=False, 
						  default=None, help='number of cores to use in batch mode (default=all available cores)')
	optional.add_argument('-h', '--help', action='help', 
						  help='show this help message and exit')
	optional.add_argument('-v', '--version', action='version', 
						  version=f'ViralRecall v. {__version__}')
	args_parser = args_parser.parse_args()

	return args_parser
	

def run_program(input : Path, 
				out_base : Path,
				database : Path, 
				window: int, 
				phagesize: int, 
				minscore: int, 
				evalue: float,
				minhit: int) -> None : 
	
	''' Main function to run ViralRecall '''

	out_base.parent.mkdir(exist_ok=True)
	logger = setup_logger(out_base.parent, str(out_base.name))
	try :
		genome_file = load_genome(input)
	except Exception as e:
		logger.error(str(e))
		return
	
	logger.info(f"Running viralrecall on {input.name} and the output is being deposited in {out_base.parent}")
	

	filt_contig_list = filt_fasta(phagesize, genome_file)
	
	if not filt_contig_list:
		logger.error(f"No contigs longer than the minimum phage size of {phagesize} bp were found in {input.name}.\n" \
		" Not proceeding with the genome. Change phagesize parameter to a smaller value if needed.")
		return

	
	proteins, description = predict_proteins(genome_file, filt_contig_list, out_base)
	
	desc_df = pd.DataFrame(description)

	gvog_hmm = database / "gvog_mirus_cat.h3m"
	ncldv_hmm = database / "NCLDV_markers.h3m"	
	
	hits = search_with_pyhmmer(proteins, gvog_hmm, evalue)
	hmm_results = parse_hmmer(hits, out_base)

	hmmout = out_base.with_suffix(".hmmout")
	hmm_df = pd.DataFrame(hmm_results) #namedtuples perfectly compatible with pandas dataframe
	
	ncldv_hits = search_with_pyhmmer(proteins, ncldv_hmm, evalue)
	ncldv_hmm_results = parse_hmmer(ncldv_hits, out_base)
	ncldv_hmm_df = pd.DataFrame(ncldv_hmm_results, columns= ["contig", "query", "HMM_hit", "bitscore", "evalue"])
	
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
	df = merge_annot(df, database)

	vir_summary = []
	Summary = namedtuple("Summary", ["file", "contig", "contig_length", "num_viral_region", 
								  "vstart", "vend", "vir_length", "num_prots", "num_viral_hits", 
								  "score", "NCLDV_markers", "Mirusvirus_hits"])
	
	viral_indices = extract_reg(window, phagesize, minscore, filt_contig_list, df, minhit)
	viral_coords = {}
	viral_df = pd.DataFrame()
	count = 0
	for key, value in viral_indices.items():		
		if len(value) >= 1:
			starts = []
			ends = []
			for idx, coords in enumerate(value, start=1):
				vregion_df = df.iloc[coords[0]:coords[1]+1]
				vregion_df.insert(0, 'Viral_region_number', f"vregion_{idx}", allow_duplicates=True)
				viral_df = pd.concat([viral_df, vregion_df], ignore_index=True)
				num_hits= len(vregion_df.loc[vregion_df['HMM_hit'] != "no_hit"])
				vstart  = vregion_df['pstart'].astype('int64')[coords[0]]
				vend = vregion_df['pend'].astype('int64')[coords[1]]
				starts.append(vstart)
				ends.append(vend)

				vreg_head : str = genome_file[key][vstart:vend].fancy_name # type: ignore
				vreg_seq : str = genome_file[key][vstart:vend].seq # type: ignore
				vreg_nuc_file = out_base.parent / f"{key}_vregion_{idx}.fna"
				
				vprots = ncldv_hmm_df.loc[ncldv_hmm_df['query'].isin(vregion_df['query'])]['HMM_hit']
				NCLDV_hits = str_hits(vprots, "NCLDV")
				mirus_hits = str_hits(vregion_df['HMM_hit'], "Mirus")

				with open(vreg_nuc_file, "w") as out:
					out.write(">" + vreg_head + "\n")
					out.write(vreg_seq)
				vir_summary.append(Summary(file=input.name,
							   contig=key,
							   contig_length= len(genome_file[key]), # type: ignore
							   num_viral_region= f"vregion_{idx}",
							   vstart= vstart,
							   vend= vend,
							   vir_length= len(vreg_seq),
							   num_prots= len(vregion_df),
							   num_viral_hits= num_hits ,
							   score= vregion_df["bitscore"].mean(),
							   NCLDV_markers= NCLDV_hits,
							   Mirusvirus_hits= mirus_hits)
							   )
				count += 1
			viral_coords[key] = [starts, ends]
		else:
			vir_summary.append(Summary(file=input.name,
							   contig=key,
							   contig_length= len(genome_file[key]), # type: ignore
							   num_viral_region= "no viral regions detected",
							   vstart= "NA",
							   vend= "NA",
							   vir_length= "NA",
							   num_prots= "NA",
							   num_viral_hits= "NA",
							   score= str(0),
							   NCLDV_markers= str_hits(ncldv_hmm_df[ncldv_hmm_df['query'] == key]['HMM_hit'], "NCLDV"),
							   Mirusvirus_hits= str_hits(df[df["query"] == key ]['HMM_hit'], "Mirus")
)
							   )
			
	summ_file = Path(str(out_base) + "_summary.tsv")
	summ_df = pd.DataFrame(vir_summary)
	vannot = Path(str(out_base) + "_viralregions.annot.tsv")
	ncldv_mark_file = Path(str(out_base) + "_NCLDV_markers.tsv")


	summ_df.to_csv(summ_file, index=False, sep= "\t")
	viral_df.to_csv(vannot, index=False, sep= "\t")
	hmm_df.to_csv(hmmout, index=False, sep= "\t")
	df.to_csv(out_base.with_suffix(".tsv"), index=False, sep="\t")
	ncldv_hmm_df.to_csv(ncldv_mark_file, index=False, sep= "\t")


	logger.info(f"Creating plots")
	plot_vreg(df, minscore, viral_coords, out_base)

	logger.info(f"Finished running Viralrecall on {input.name} and found {count} viral region/s" )
	return




def main(argv=None):

	args_list = parse_args()

	# set up object names for input/output/database folders
	input =    Path(args_list.input).expanduser() 
	project = Path(args_list.outdir).expanduser()
	db_dir = Path(args_list.database).expanduser() 
	window = int(args_list.window)*1000 # convert to bp
	phagesize = int(args_list.minsize)*1000
	minscore = int(args_list.minscore)
	minhit = int(args_list.minhit)
	evalue = float(args_list.evalue)
	cpus = args_list.cpus

	
	
	database = find_db(db_dir)
	prep_hmm(database)

	existence = input.exists()
	indir = input.is_dir()

	if indir and existence:
		if check_directory_permissions(input) == False :
			print(f"Insufficient read/write permissions for {input} \nQuitting")
			return
		
		project.mkdir(exist_ok=True)
		try :
			file_list = [i.name for i in input.iterdir() if i.name.endswith(('.fna', '.fasta', '.fa'))]
		except FileNotFoundError as f:
			raise FileNotFoundError(f"Input Directory : {input}, does not contain fasta files. Quitting")
		
		cpus : int = mp_cpu(cpus)
		print(f"Running ViralRecall in batch mode on {len(file_list)} files found in {input} and parallelizing across {cpus} cores")
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
			print(f"Finished running ViralRecall on all files in {input} and the output has been deposited in {project}")

	elif existence and not indir:
		
		out_base = project / input.stem
		run_program(input, out_base, database, window, phagesize, minscore, evalue, minhit)

	else:
		print("Input is not a valid directory or file. Please check the input path.")
	

if __name__ == '__main__':
	main()
	

