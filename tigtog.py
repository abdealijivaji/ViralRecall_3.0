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
from collections import Counter, defaultdict


input_gen = Path("/home/abdeali/viralR_test_output/Chlamy_punui/contig_536_vregion_1.fna")
vreg_tab = pd.read_csv("/home/abdeali/viralR_test_output/Chlamy_punui/Chlamy_punui_contig_viralregions.annot.tsv", sep= "\t", header=0)

# print(vreg_tab)

vreg = Fasta(input_gen)

def get_gc(vreg : Fasta) -> float :
    seq_len = len(vreg[0])
    seq_np = np.asarray(vreg[0])
    cnt = np.count_nonzero(seq_np == b"G") + np.count_nonzero(seq_np == b"C")
    GC_perc = round(cnt * 100 / seq_len, 2)
    return GC_perc
    

def density(table: pd.DataFrame, vreg: Fasta) -> float :
    prot_lens = pd.Series
    prot_lens = table["pend"] - table["pstart"]
    sum_len = sum(prot_lens)
    dens = round(100 * sum_len / len(vreg[0]), 2)
    return dens

dens = density(vreg_tab, vreg)

def large_seq(vreg: Fasta) -> bool :
    if len(vreg[0]) > 5_000_000 :
        return True
    else:
        return False

def imp_names(list_file : Path) -> list :
    list_g = list()
    g_names = open(list_file,"r")
    for line in g_names.readlines(): 
        g = line.rstrip()
        g2 = re.sub("hmm$", "trim", g)
        # g2 = re.sub("_", "" , g2)
        list_g.append(g2)
    return list_g

# imp_name_list = imp_names()

def parse_hits(table, imp_names) -> dict :
    gvg_hits = {i : 0 for i in imp_names}
    for hmm, count in table["HMM_hit"].value_counts().items() :
        hmm = re.sub("GVOGm", "GVOGm_", hmm)
        hmm = hmm + ".trim"
        if hmm in imp_names:
            gvg_hits[hmm] = count
    return gvg_hits


# hitdict = parse_hits(vreg_tab, imp_name_list)
# print(imp_name_list)

def create_df(gc_perc, hitdict) :
    
    df = pd.DataFrame(hitdict, index=["seq",])
    df.insert(0, "GC_content", gc_perc)
    print(df)

# create_df(get_gc(vreg), hitdict)


def tax_predict(clf_file, input_df) :
    
    clf = joblib.load(clf_file)
    pred = clf.predict(input_df)
    conf_pred = clf.predict_proba(input_df)
    class_labels = clf.classes_
    for i , prediction in enumerate(conf_pred) :
        max_index = np.argmax(prediction)
        confidence = prediction[max_index]
    
    return pred, confidence


def run_tigtog(vreg: Fasta, table: pd.DataFrame) :

    GC_perc = get_gc(vreg)

    both_levels =  imp_names(Path("/home/abdeali/hmm/db/names_imp_GVOGs_both_levels.csv"))
    order_levels = imp_names(Path("/home/abdeali/hmm/db/names_imp_gvog_order.csv"))
    fam_levels = imp_names(Path("/home/abdeali/hmm/db/names_imp_gvog_fam.csv"))
    
    hitdict = parse_hits(table, both_levels)
    input_df = create_df(GC_perc, hitdict)

    ord_file = Path("/home/abdeali/hmm/clf/clf_Order_final.joblib")
    fam_file = Path("/home/abdeali/hmm/clf/clf_Fam_final.joblib")

    Ord_pred, Ord_prob = tax_predict(ord_file, input_df)
    print(f"Predicted order is {Ord_pred} with {Ord_prob} % probability")
    


# run_tigtog(vreg, vreg_tab)

joblib.load("/home/abdeali/hmm/clf/clf_Order_final.joblib")