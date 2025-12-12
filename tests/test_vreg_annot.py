import pandas as pd
import numpy as np
from pathlib import Path
import pytest

from viralrecall import vreg_annot


def test_above_threshold_ind():
    s = pd.Series([1, 5, 10, 3])
    res = vreg_annot.above_threshold_ind(s, 4)
    assert res.tolist() == [False, True, True, False]


def test_contiguous_true_ranges_numpy():
    s = pd.Series([False, True, True, False, True, True, True, False])
    ranges = vreg_annot.contiguous_true_ranges_numpy(s)
    assert ranges == [(1, 2), (4, 6)]


def test_sliding_window_mean_small_window():
    df = pd.DataFrame({
        'contig': ['c1'] * 5,
        # pstart in seconds
        'pstart': [0, 1, 2, 3, 4],
        'bitscore': [0, 10, 20, 30, 40],
    })
    # window = 3 seconds, ensure rolling collects at least 3 data points
    window = pd.offsets.Second(3)
    means = vreg_annot.sliding_window_mean(df, window)
    assert isinstance(means, pd.Series)
    # There should be non-zero values where rolling applies
    assert means.max() > 0


def make_test_df():
    # create a df for extract_reg
    df = pd.DataFrame({
        'contig': ['c1']*6,
        'pstart': [10,20,30,40,50,60],
        'pend': [19,29,39,49,59,69],
        'HMM_hit': ['Mirus_A','Mirus_B','no_hit','Mirus_C','Mirus_D','Mirus_E'],
        'rollscore': [0,10,12,5,15,20]
    })
    return df


def test_extract_reg_simple():
    df = make_test_df()
    # We set phagesize so that it passes single region test, window small, minscore=10
    filt_contig_list = ['c1']
    res = vreg_annot.extract_reg(window=10, phagesize=1, minscore=10, filt_contig_list=filt_contig_list, df=df, minhit=1)
    assert 'c1' in res


def test_merge_annot(tmp_path):
    # create gvog_annotation.tsv file
    d = tmp_path / 'hmm'
    d.mkdir()
    tsv = d / 'gvog_annotation.tsv'
    tsv.write_text('HMM_hit\tannotation\nMirus_A\tSomeAnnotation\n')
    df = pd.DataFrame({'HMM_hit':['Mirus_A','no_hit']})
    merged = vreg_annot.merge_annot(df, d)
    assert 'annotation' in merged.columns


def test_str_hits():
    s = pd.Series(['Mirus_A','NCLDV_marker','no_hit','Mirus_B'])
    assert vreg_annot.str_hits(s, 'Mirus_').startswith('Mirus_')
    assert vreg_annot.str_hits(s, 'XYZ_') == 'NA'
