"""
Codonw(-barebones) regression tests
"""

import os
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
    out_fn = "{}/input.out".format(path)

    with open(os.devnull, 'w') as devnull:
        subprocess.check_call(
            "{p}/../bin/codonw {p}/input.fna {o} /dev/null "
            "-all_indices -aro -hyd -nomenu -silent -nowarn -noblk -machine".format(p=path, o=out_fn),
            shell=True, stdout=devnull, stderr=devnull)
            
    # read ref and test output
    df_ref = pd.read_csv(out_fn, index_col="title",
                         sep='[\\s,]+', engine='python')
    df_ref.index = df_ref.index.str.split(' ').str[0]
    df_ref.dropna(axis=1, inplace=True)

    df_test = pd.read_csv("{}/test.out".format(path), index_col="title")

    # compare
    print("Checking...", end=" ")
    for c in df_ref.columns:
        assert np.array_equal(df_ref[c], df_test[c])
        print(c, end=", ")
    print()

    return


@pytest.mark.parametrize("blktype", [
    "dinuc",
    "aau",
    "raau",
    "cu",
    "cutab",
    "cutot",
    "rscu",
    "base"
])
def test_bulk_regression(blktype):
    with open(os.devnull, 'w') as devnull:
        subprocess.check_call(
            "{p}/../bin/codonw {p}/input.fna /dev/null {p}/{b}.blk "
            "-nomenu -silent -nowarn -machine -{b}".format(p=path, b=blktype),
            shell=True, stdout=devnull)

    return
