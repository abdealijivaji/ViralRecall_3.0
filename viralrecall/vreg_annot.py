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

def above_threshold_ind(rollscore : pd.Series, minscore : int) -> pd.Series: 
	bool_array : pd.Series = rollscore > minscore
	return bool_array

def contiguous_true_ranges(bool_array : pd.Series) -> list[tuple[int, int]]:
	s = bool_array
	s1 = s != s.shift()
	starts = s & s1
	ends = s & (~s).shift(-1, fill_value=True)
	return list(zip(starts[starts].index, ends[ends].index))


def extract_reg(window : int, phagesize : int, 
				minscore : int, filt_contig_list: list[str], 
				df : pd.DataFrame , minhit : int) -> dict:
	viral_indices = { i : [] for i in filt_contig_list }  # dictionary to hold indices of viral regions for each contig
	for name, grp in df.groupby('contig'):
		name = str(name)
		above_threshold = above_threshold_ind(grp['rollscore'], minscore)
		
		viral_ranges = contiguous_true_ranges(above_threshold)
		if len(viral_ranges) > 0:
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
						if vend - vstart >= phagesize and len(unq_hit) > minhit:
							viral_indices[name].append([strt_idx, end_idx])
						strt_idx = nxt_strt_idx
						vstart = nxt_vstart
					
				unq_hit = grp['HMM_hit'].iloc[strt_idx:end_idx].unique()
				if vend - vstart >= phagesize and len(unq_hit) > minhit:
						viral_indices[name].append([strt_idx, end_idx])
	return viral_indices



def merge_annot(df: pd.DataFrame, db_dir: Path) -> pd.DataFrame :
	
	annot = pd.read_table(Path(db_dir / "gvog_annotation.tsv"), delimiter="\t")
	return df.merge(annot, on = 'HMM_hit', how='left')

def str_hits(HMM_series: pd.Series, mode: str) -> str:
	NCLDV = {'A32', 'D5', 'mcp', 'mRNAc', 'PolB', 'RNAPL', 'RNAPS', 'RNR', 'SFII', 'VLTF3'}
	Mirus = {'Mirus_JellyRoll' ,  'Mirus_MCP' ,  'Mirus_Portal' ,  'Mirus_Terminase_ATPase' ,  'Mirus_Terminase_merged' ,  'Mirus_Triplex1' ,  'Mirus_Triplex2'}
	
	hmm_hits = set(HMM_series.to_list())
	
	markers = []
	if mode == 'NCLDV':
		markers = list(hmm_hits.intersection(NCLDV))
	elif mode == 'Mirus':
		markers = list(Mirus.intersection(hmm_hits))

	if not markers:
		return "NA"
	else:
		return ";".join(markers)