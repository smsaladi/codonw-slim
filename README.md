[![Build Status](https://travis-ci.org/smsaladi/codonw-slim.svg?branch=master)](https://travis-ci.org/smsaladi/codonw-slim)
[![PyPI version](https://badge.fury.io/py/codonw-slim.svg)](https://badge.fury.io/py/codonw-slim)
![PyPI - Downloads](https://img.shields.io/pypi/dm/codonw-slim)

codonw-slim
===========

CodonW is a package for codon usage analysis written by John Peden in
Paul Sharp's group (University of Nottingham).
It was originally designed to simplify multivariate analysis of codon usage
with other metrics related to codon usage calculated.

codonw-slim refactors the original codebase to add Python bindings to the
underlying methods written in C and to focus on the *other* codon usage
metrics. See below for a list of metrics and their usage.

A detailed description of each metric, with references, can be found in
in the docstrings of Python methods (`codonw/codonwlib/codonw.pyx`).
All of the multivariate analysis code has been removed since this sort of
analysis is more easily done in a higher level language
(e.g. [FactoMineR](https://cran.r-project.org/web/packages/FactoMineR/index.html)).
The interative interface has also been removed.

No error checking of nucleotide sequences is done, e.g. for start,
stop codons, internal stops codons, non-translatable, and partial codons.
Users should do this themself to inputs provided. For more information about
how amino acids and codons have are represented internally (`Recoding.md`).

The source code and releases for codonw-slim can be obtained from
https://www.github.com/smsaladi/codonw-slim. Please report bugs and improvements
via pull requests at this repository. All modifications must pass regression
testing.


## Build and Installation

```bash
pip install codonw-slim
```

## Usage

The following metrics are available:

* codon adaptation index (CAI)
    - `CodonSeq.cai`
* frequency of optimal codons (Fop)
    - `CodonSeq.fop`
* codon bias index (CBI)
    - `CodonSeq.cbi`
* the effective number of codons (Nc)
    - `CodonSeq.enc`
* hydropathicity of protein
    - `CodonSeq.hydropathy`
* aromaticity score
    - `CodonSeq.aromaticity`
* Silent base composition (including GC3s)
    - `CodonSeq.silent_base_usage`
* Codon & Amino acid usage (count and relative)
    - `CodonSeq.codon_usage`
    - `CodonSeq.aa_usage`
    - `CodonSeq.rscu`
    - `CodonSeq.raau`
* Base composition by codon position
    - `CodonSeq.bases`
* Base composition in all frames,
    Length of gene,
    Number of synonymous codons,
    G+C content (overall and by codon position),
    G+C content of synonymous codons at the 3rd position,
    G+C content of non-synonymous codons at the 3rd position,
    Number of synonymous codons,
    Number of amino acids
    - `CodonSeq.bases2`
* Dinucleotide count by frame
    - `CodonSeq.dinuc`

As written above, each is a method of the `codonw.CononSeq` object, e.g.

```python
import codonw
cseq = codonw.CodonSeq("ATGAATATGCTCATTGTCGGTAGAGTTGTTGCTAGTGTTGGGGGAAGCGGACTTCAAACG")
cseq.cai()
```

The return type can be a simple value, `pd.Series`, or `pd.DataFrame`.

The genetic codes can be specified by setting the `CodonSeq.genetic_code`
property with a `pd.Series` whose index is a codon and value is the single
letter amino acid. Instantiate an object and see `CodonSeq.genetic_code`
for more details.

Some indicies have an option of reference values to choose from (e.g. `CodonSeq.fop`).
Several references values can be chosen by specifying the corresponding integer.
If you'd like to have user-provided reference values, please implement this
functionality and make a pull-request.


## Why the name codonW?

Excerpted directly from John Peden's CodonW README...

Well first you must realise that "clustal" (a very popular multiple
alignment program by Des Higgins) was originally written in Paul's lab in
Trinity College Dublin. Clustal has since been rewritten from FORTRAN into
C and undergone several name changes clustal-> clustalv-> clustalw ->
clustalx. There was also a program called "codons" written in FORTRAN by
Andrew Lloyd (a post-doc in Paul's lab), this was the original inspiration
for codonW. An early version of codonW, written in C, was called codonv.
When the code was enhanced to include multivariate analysis, what better
name than codonW.
