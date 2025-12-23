'''
This is a modified version of TIGTOG originally written by Dr. Anh D. Ha and is available on <https://github.com/anhd-ha/TIGTOG>
'''

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
import joblib
import csv, re
from pyfaidx import Fasta
from pathlib import Path
from collections import Counter


input_gen = Path("/home/abdeali/viralR_test_output/Chlamy_punui/contig_536_vregion_1.fna")
vreg_tab = pd.read_csv("/home/abdeali/viralR_test_output/Chlamy_punui/Chlamy_punui_contig_viralregions.annot.tsv", sep= "\t", header=0)

# print(vreg_tab)

vreg = Fasta(input_gen)

def get_gc(vreg) :
    seq_len = len(vreg[0])
    seq_np = np.asarray(vreg[0])
    cnt = np.count_nonzero(seq_np == b"G") + np.count_nonzero(seq_np == b"C")
    GC_perc = round(cnt * 100 / seq_len, 2)
    return GC_perc
    

def density(table: pd.DataFrame, vreg) :
    prot_lens = pd.Series
    prot_lens = table["pend"] - table["pstart"]
    sum_len = sum(prot_lens)
    dens = round(100 * sum_len / len(vreg[0]), 2)
    return dens

def large_seq(vreg) :
    if len(vreg[0]) > 5_000_000 :
        return True
    else:
        return False

def imp_names() :
    list_g = list()
    g_names = open(Path("~/hmm/db/names_imp_GVOGs_both_levels.csv").expanduser(),"r")
    for line in g_names.readlines(): 
        g = line.rstrip()
        g2 = re.sub(".hmm$", "", g)
        g2 = re.sub("_", "" , g2)
        list_g.append(g2)
    return list_g

imp_name_list = imp_names()

def parse_hits(table) :
    gvg_cnt = table["HMM_hit"].value_counts()
    print(gvg_cnt)

hit_dict = parse_hits(vreg_tab)
print(hit_dict)
    



dens = density(vreg_tab, vreg)


