"""

CodonW-slim regression tests

"""

import os
import io
import warnings

import numpy as np
import pandas as pd

import pytest
import Bio.SeqIO

import codonw


def read_fasta(fn):
    """Reads from a file (string or file handle) assuming format specified
    """
    seqs = [(r.id, str(r.seq)) for r in Bio.SeqIO.parse(fn, "fasta")]
    seqs = list(zip(*seqs))
    seqs = pd.Series(seqs[1], index=seqs[0], name="seq")
    return seqs

# location of *this* script
path = os.path.dirname(os.path.realpath(__file__))
seq_fn = "{}/input.fna".format(path)
test_seqs = read_fasta(seq_fn)


def compare_df(df_ref, df_test, name):
    # compare
    print("Checking {} ...".format(name), end=" ")
    missing = []
    for c in df_ref.columns:
        if c in df_test.columns:
            print(c, end=", ")
            # Try printing error values before throwing error
            close_vals = np.isclose(df_ref[c], df_test[c], atol=.001, rtol=.01)
            if (~close_vals).sum() > 0:
                print(df_ref.loc[~close_vals, c], df_test.loc[~close_vals, c])
            np.testing.assert_allclose(df_ref[c], df_test[c], atol=.001, rtol=.01)
        else:
            missing.append(str(c))
    print()

    if missing:
        warnings.warn("Missing: {} from {}".format(", ".join(missing), name), UserWarning)
    
    return

def test_index_regression():
    df_test = test_seqs.to_frame()
    df_test['seqw'] = df_test['seq'].apply(lambda x: codonw.CodonSeq(x))
    df_test['CAI'] = df_test['seqw'].apply(lambda x: x.cai())
    df_test['CBI'] = df_test['seqw'].apply(lambda x: x.cbi())
    df_test['Fop'] = df_test['seqw'].apply(lambda x: x.fop())
    df_test['Nc'] = df_test['seqw'].apply(lambda x: x.enc())
    df_test['Gravy'] = df_test['seqw'].apply(lambda x: x.hydropathy())
    df_test['Aromo'] = df_test['seqw'].apply(lambda x: x.aromaticity())
    df_test[['T3s', 'C3s', 'A3s', 'G3s']] = df_test['seqw'].apply(
        lambda x: x.silent_base_usage()).apply(pd.Series)
            
    # read ref and test output
    df_ref = pd.read_csv("{}/ref/input.out".format(path), index_col="title",
                         sep='[\\s,]+', engine='python')
    df_ref.index = df_ref.index.str.split(' ').str[0]
    df_ref.dropna(axis=1, inplace=True)

    compare_df(df_ref, df_test, "indicies")
    return


@pytest.mark.parametrize("blk", [
    "aau",
    "raau",
    "base",
])
def test_bulk_regression_a(blk):
    ref_fn = "{p}/ref/input.{b}.blk".format(p=path, b=blk)
    df_ref = pd.read_csv(ref_fn, index_col=0, sep='[\\s,]+', engine='python')
    df_ref.index = df_ref.index.str.split(' ').str[0]
    df_ref.dropna(axis=1, inplace=True)

    df_test = test_seqs.to_frame()
    df_test['seqw'] = df_test['seq'].apply(lambda x: codonw.CodonSeq(x))
    if blk == "aau":
        df_ref.rename(columns=codonw.aa3_aa1.to_dict(), inplace=True)
        def func(x): return x.aa_usage()
    elif blk == "raau":
        df_ref.rename(columns=codonw.aa3_aa1.to_dict(), inplace=True)
        def func(x): return x.raau()
    elif blk == "base":
        def func(x): return x.gc()
    
    df_test = df_test['seqw'].apply(lambda x: func(x))
    compare_df(df_ref, df_test, blk)
    return

@pytest.mark.parametrize("blk", [
    "cu",
    "rscu"
])
def test_bulk_regression_b(blk):
    ref_fn = "{p}/ref/input.{b}.blk".format(p=path, b=blk)
    def read_table64(f):
        return pd.read_csv(f, sep='[\\s,]+', header=None,
                        usecols=range(8), index_col=False, engine='python')
    df_ref = read_table64(ref_fn)

    df_test = test_seqs.to_frame()
    df_test['seqw'] = df_test['seq'].apply(lambda x: codonw.CodonSeq(x))

    if blk == "cu":
        def func(x): return x.codon_usage()
    elif blk == "rscu":
        def func(x): return x.rscu()
    
    df_out = df_test['seqw'].apply(lambda x: func(x).values.reshape([4, 16]))    
    df_out = pd.DataFrame(np.vstack(df_out.values))

    compare_df(df_ref, df_out, blk)
    return


def test_bulk_regression_dinuc():
    ref_fn = "{p}/ref/input.{b}.blk".format(p=path, b='dinuc')
    def read_3frame_col(f):
        df = pd.read_csv(f, sep='[\\s,]+', header=0, index_col='title', engine='python')
        df.drop(columns=[c for c in df.columns if 'frame' in c], inplace=True)
        return df
    df_ref = read_3frame_col(ref_fn)

    df_test = test_seqs.to_frame()
    df_test['seqw'] = df_test['seq'].apply(lambda x: codonw.CodonSeq(x))

    def dinuc_ser(x):
        return pd.Series(x.dinuc().values.reshape([16 * 4]),
                         index=df_ref.columns)

    df_out = df_test['seqw'].apply(dinuc_ser)

    compare_df(df_ref, df_out, "dinuc")
    return
