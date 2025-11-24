import os
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd


	
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



def merge_annot(df: pd.DataFrame) -> pd.DataFrame :
	full_tab = pd.read_table("hmm/gvog_annotation.tsv", delimiter="\t")
	annot = full_tab.loc[:, ["GVOG", "NCVOG_descs"]]
	annot = annot.rename(columns={'GVOG':'HMM_hit'})
	df = df.merge(annot, on = 'HMM_hit', how='left')
	return df

def count_hits(hits, contigs):
	''' Right now does listing of markers for each contig '''
	query2hits = defaultdict(list)
	for i in contigs:
		for j in hits:
			if i == j.contig:
				marker_score = j.HMM_hit + ": " + str(j.bitscore)
				query2hits[i].append(marker_score)

	return query2hits
		