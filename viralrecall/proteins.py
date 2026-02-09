import pyrodigal_gv
from pathlib import Path
from pyhmmer import easel, plan7, hmmer
from collections import namedtuple
from pyfaidx import Fasta
from typing import Any

amino = easel.Alphabet.amino()

def predict_proteins(input : Fasta, 
					 contigs : list[str], 
					 outbase: Path) -> tuple[Any, list[tuple]]:
	''' 
	Predict proteins using pyrodigal-gv and returns namedtuples of proteins and their headers 
	Also writes predictions in CDS, AA, and gff3 formats.
	'''
	orf_finder = pyrodigal_gv.ViralGeneFinder(meta=True)
	prot_out = outbase.with_suffix(".faa") 
	cds_out = outbase.with_suffix(".cds.fasta")
	gff_out = outbase.with_suffix(".gff")
	gbk_out = outbase.with_suffix(".gbk")

	prot, cds, gff, gbk = open(prot_out, "w") , open(cds_out, "w") , open(gff_out, "w"), open(gbk_out, "w")
	
	proteins : easel.TextSequenceBlock = easel.TextSequenceBlock()
	header : list[tuple] = []
	Header = namedtuple("Header", ["contig", "query", "pstart", "pend", "pstrand", "gen_code"])
	
	with prot, cds, gff, gbk:
		for seqrecord in contigs:
			sequence = bytes(str(input[seqrecord]), 'UTF-8') # loading the sequence as bytes
			genes = orf_finder.find_genes(sequence)
			genes.write_translations(prot, sequence_id=seqrecord)
			genes.write_genes(cds, sequence_id=seqrecord)
			genes.write_gff(gff, sequence_id=seqrecord)
			genes.write_genbank(gbk, sequence_id=seqrecord)
			
			for n, gene in enumerate(genes, start= 1):
				aa = gene.translate()
				prot_id = seqrecord + "_" + str(n)
				head = Header(seqrecord, prot_id, gene.begin, gene.end, gene.strand, gene.translation_table)
				header.append(head)
	
				if len(aa) < 100000: # HMMER can't take proteins longer than 100k aa
					proteins.append(easel.TextSequence(name = bytes(prot_id,  'UTF-8'), sequence= aa))
	return proteins, header

def subset_proteins(proteins: easel.TextSequenceBlock, 
					prot_ids) :
	''' Returns a subset of proteins from the TextSequenceBlock based on a list of protein IDs '''

	for prot_id in prot_ids :
		sub_prot = easel.TextSequence(name= proteins.indexed[prot_id].name,
							sequence= proteins.indexed[prot_id].sequence)
		yield sub_prot


tbl_head = b"#                                                               --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----\ntarget_name        accession  query_name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc target_description\n"

def search_with_pyhmmer(proteins: easel.TextSequenceBlock, 
						hmm_path: Path, 
						evalue: float) :
	
	
	digseqs = easel.DigitalSequenceBlock(amino, proteins.digitize(amino))
	
	with plan7.HMMFile(hmm_path) as hmm_file:
		hits = list(hmmer.hmmsearch(hmm_file, digseqs, E=evalue))
	return hits	
	
		
def parse_hmmer(hits, out_base: Path) -> list[tuple]:	
	results : list[tuple] = []
	Result = namedtuple("Result", ["contig", "query", "HMM_hit", "bitscore", "evalue"])	
	# hmmout = out_base.with_suffix(".tblout")

	# with open(hmmout, 'wb') as outfile:
	# 	hits = list(hits)
	for hitlist in hits:
		# hitlist.write(outfile, header=True) 
		for hit in hitlist:
			if hit.included:
				hmm_hit = hit.name.decode() if hasattr(hit.name, "decode") else hit.name
				Contig =  hmm_hit.rsplit("_", maxsplit=1)[0]
				hitlist_name = hitlist.query.name.decode() if hasattr(hitlist.query.name, "decode") else hitlist.query.name
				eval = "%.3g" % hit.evalue
				results.append(Result(
					Contig ,
					hmm_hit,
					hitlist_name,
					round(hit.score, 2),
					eval
				))
	return results

