"""
Quick regression test of codonw
"""

import os
import subprocess

import numpy as np
import pandas as pd

# location of *this* script
path = os.path.dirname(os.path.realpath(__file__))

with open(os.devnull, 'w') as devnull:
    # compile 
    subprocess.check_call("make -C {}/../ clean && make -j4 -C {}/../".format(path, path),
                          shell=True, stderr=devnull)

    # run test file
    subprocess.check_call(
        "{}/../bin/codonw {}/input.dat test blk "
        "-enc -gc -gc3s -L_aa -cai -dinuc -aau -coa_expert > {}/test.out".format(path, path, path),
        shell=True, stderr=devnull)

def test_regression():
    # read ref and test output
    df_ref = pd.read_csv("{}/ref/out.txt".format(path), index_col="title")
    df_ref.index = df_ref.index.str.split(' ').str[0]
    df_ref.dropna(axis=1, inplace=True)

    df_test = pd.read_csv("{}/test.out".format(path), index_col="title")

    # compare
    for c in df_ref.columns:
        assert np.array_equal(df_ref[c], df_test[c])

