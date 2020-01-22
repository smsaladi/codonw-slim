"""
Codonw(-barebones) regression tests
"""

import os
import io
import subprocess

import numpy as np
import pandas as pd

import pytest

# location of *this* script
path = os.path.dirname(os.path.realpath(__file__))

with open(os.devnull, 'w') as devnull:
    # compile 
    subprocess.check_call("make -C {p}/../ clean && make -j4 -C {p}/../".format(p=path),
                          shell=True, stderr=devnull)


def test_index_regression():
    out_fn = "{}/test.out".format(path)

    with open(os.devnull, 'w') as devnull:
        cmd = "{p}/../bin/codonw {p}/input.fna /dev/null /dev/null " \
              "-all_indices -aro -hyd -nomenu -silent -nowarn -noblk -machine > {o}".format(p=path, o=out_fn)
        subprocess.check_call(cmd, shell=True, stderr=devnull)
            
    # read ref and test output
    df_ref = pd.read_csv("{}/ref/input.out".format(path), index_col="title",
                         sep='[\\s,]+', engine='python')
    df_ref.index = df_ref.index.str.split(' ').str[0]
    df_ref.dropna(axis=1, inplace=True)

    df_test = pd.read_csv(out_fn, index_col="title")

    # compare
    print("Checking...", end=" ")
    for c in df_ref.columns:
        np.testing.assert_allclose(df_ref[c], df_test[c])
        print(c, end=", ")
    print()

    return


@pytest.mark.parametrize("blktype", [
    "aau",
    "raau",
    "base",
    "dinuc",
    "cu",
    "rscu"
])
def test_bulk_regression(blktype):
    blk_fn = "{p}/{b}.blk".format(p=path, b=blktype)
    with open(os.devnull, 'w') as devnull:
        subprocess.check_call(
            "{p}/../bin/codonw {p}/input.fna /dev/null {bfn} "
            "-nomenu -silent -nowarn -machine -{b}".format(p=path, bfn=blk_fn, b=blktype),
            shell=True, stdout=devnull)

    ref_fn = "{p}/ref/input.{b}.blk".format(p=path, b=blktype)
    # Certain bulk files need to be parsed differently
    if blktype in ["aau", "raau", "base"]:
        df_ref = pd.read_csv(ref_fn, index_col=0, sep='[\\s,]+', engine='python')
        df_ref.index = df_ref.index.str.split(' ').str[0]
        df_ref.dropna(axis=1, inplace=True)
        df_test = pd.read_csv(blk_fn, index_col=0, sep='[\\s,]+', engine='python')

    elif blktype in ["dinuc"]:
        def read_3frame_col(f):
            df = pd.read_csv(f, sep='[\\s,]+', header=0, index_col='title', engine='python')
            df.drop(columns=[c for c in df.columns if 'frame' in c], inplace=True)
            return df
        df_ref = read_3frame_col(ref_fn)

        with open(blk_fn, 'r') as fh:
            blk_txt = fh.read()
            blk_txt = blk_txt.replace(',\n', '\n')
        df_test = read_3frame_col(io.StringIO(blk_txt))

    elif blktype in ["cu", "rscu"]:
        def read_table64(f):
            return pd.read_csv(f, sep='[\\s,]+', header=None,
                            usecols=range(8), index_col=False, engine='python')
        df_ref = read_table64(ref_fn)
        df_test = read_table64(blk_fn)

    else:
        raise ValueError("unknown blktype: `" + blktype + "`")

    # compare
    print("Checking...", end=" ")
    for c in df_ref.columns:
        np.testing.assert_allclose(df_ref[c], df_test[c])
        print(c, end=", ")
    print()

    # if sucessful, remove test output file
    os.unlink(blk_fn)

    return
