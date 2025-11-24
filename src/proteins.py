import pyrodigal_gv
from pathlib import Path
from pyhmmer import easel, plan7, hmmer
from collections import defaultdict, namedtuple
from pyfaidx import Fasta



amino = easel.Alphabet.amino()

def predict_proteins(input : Fasta, 
					 contigs : list[str], 
					 outbase: Path) -> tuple[list, list]:
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
	
	proteins : list = []
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
					proteins.append(easel.TextSequence(name = bytes(prot_id,  'UTF-8'), sequence= aa).digitize(amino))
	return proteins, header

tbl_head = b"#                                                               --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----\ntarget_name        accession  query_name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc target_description\n"

def search_with_pyhmmer(proteins: list, 
						hmm_path: Path, 
						out_base: Path, 
						evalue: float) -> list[tuple]:
	
	results : list[tuple] = []
	Result = namedtuple("Result", ["contig", "query", "HMM_hit", "bitscore", "evalue"])	
	
	digseqs = easel.DigitalSequenceBlock(amino, proteins)
	
	hmmout = out_base.with_suffix(".tblout")
	
	with plan7.HMMFile(hmm_path) as hmm_file:
		# Convert protein predictions into pyhmmer Sequence objects
		# for hmm in hmm_file:
		# 	print(hmm.cutoffs)
		hits = hmmer.hmmsearch(hmm_file, digseqs, E=evalue)
		with open(hmmout, 'wb') as outfile:
			outfile.write(tbl_head)
			for hitlist in hits:
				hitlist.write(outfile, header=False) 
				for hit in hitlist:
					if hit.included:
						Contig =  hit.name.decode().rsplit("_", maxsplit=1)[0]
						eval = "%.3g" % hit.evalue
						results.append(Result(
							Contig ,
							hit.name.decode() ,
							hitlist.query.name.decode(),
							round(hit.score, 2),
							eval
						))
	return results