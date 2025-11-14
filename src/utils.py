from pathlib import Path
from pyfaidx import Fasta
import numpy as np
import pandas as pd

valid_bases = set('ATCGN')

def filt_fasta(phagesize : int, genome_file: Fasta) -> list:
	filt_contig_list = []
	
	for contig in genome_file.keys():
		if len(genome_file[contig]) >= phagesize :
			filt_contig_list.append(contig)
	
	return filt_contig_list


	
def is_DNA(genome : Fasta) -> bool :
	seq_for_check = set(str(genome[0][:1000]).upper())
	if set(seq_for_check).issubset(valid_bases) :
		return True
	else :
		return False
	
def above_threshold_ind(rollscore : pd.Series, minscore) :
	bool_array : list = rollscore > minscore
	return bool_array

def contiguous_true_ranges_numpy(data):
	data = np.array(data)
	edges = np.diff(data.astype(int))
	starts = np.where(edges == 1)[0] + 1
	ends = np.where(edges == -1)[0]
	if data[0]:
		starts = np.insert(starts, 0, 0)
	if data[-1]:
		ends = np.append(ends, len(data) - 1)
	return list(zip(starts, ends))


def extract_reg(window, phagesize, minscore, filt_contig_list, df):
	viral_indices = { i : [] for i in filt_contig_list }  # dictionary to hold indices of viral regions for each contig
	for name, grp in df.groupby('contig'):
		above_threshold = above_threshold_ind(grp['rollscore'], minscore)
		
		viral_ranges = contiguous_true_ranges_numpy(above_threshold)

		strt_idx = viral_ranges[0][0]  # start of the first group
		vstart = grp['pstart'].astype('int64')[strt_idx]  # start position of the first group
		
		if len(viral_ranges) == 1:
			end_idx = viral_ranges[0][-1] 
			vend = grp['pend'].astype('int64')[end_idx]  # end position of the first group
				
			if vend - vstart >= phagesize:
				viral_indices[name].append([strt_idx, end_idx])
				
		else:
			# If difference between groups is less than 5 proteins and less than window size (default 15kb)
			# We don't update vstart, 
			for i in range(1, len(viral_ranges)):
				nxt_strt_idx = viral_ranges[i][0]
				end_idx = viral_ranges[i-1][-1]
				prot_diff = nxt_strt_idx - end_idx
				nxt_vstart = grp['pstart'].astype('int64')[nxt_strt_idx]
				vend = grp['pend'].astype('int64')[end_idx]
				bp_diff = nxt_vstart - vend
						
			if prot_diff > 5 and bp_diff > window :
				unq_hit = grp['HMM_hit'].iloc[strt_idx:end_idx].unique()			
				if vend - vstart >= phagesize and len(unq_hit) > 3:
					viral_indices[name].append([strt_idx, end_idx])
				strt_idx = nxt_strt_idx
				vstart = nxt_vstart
	return viral_indices