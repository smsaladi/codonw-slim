"""
Wrappers around CodonW C functions
"""

from ctypes import c_int, c_long, c_float

import numpy as np
cimport numpy as np
np.import_array()

cimport codonwlib

cpdef np.ndarray[dtype=float, ndim=1, mode="c"] rscu(char* seq, int genetic_code=0):
    """Calculates Relative Synonymous Codon Usage
    """
    cdef codonwlib.GENETIC_CODE_STRUCT ref_code = codonwlib.cu_ref[genetic_code]

    cdef long codon_tot = 0
    cdef int valid_stops = 0
    cdef np.ndarray[dtype=long, ndim=1, mode="c"] ncod = np.zeros([65], dtype=c_long)
    cdef np.ndarray[dtype=long, ndim=1, mode="c"] naa = np.zeros([23], dtype=c_long)
    codonwlib.codon_usage_tot(seq, &codon_tot, &valid_stops, &ncod[0], &naa[0], &ref_code)

    cdef np.ndarray[dtype=int, ndim=1, mode="c"] dds = np.zeros([65], dtype=c_int)
    codonwlib.how_synon(&dds[0], &ref_code)

    cdef np.ndarray[dtype=float, ndim=1, mode="c"] rscu_vals = np.zeros([65], dtype=c_float)
    codonwlib.rscu_usage(&ncod[0], &naa[0], &rscu_vals[0], &dds[0], &ref_code)

    return rscu_vals

cpdef double cai(char* seq, int genetic_code=0, int cai_ref=0):
    """Calculates Codon Adaptation Index
    """
    cdef codonwlib.GENETIC_CODE_STRUCT ref_code = codonwlib.cu_ref[genetic_code]
    cdef codonwlib.CAI_STRUCT ref_cai = codonwlib.cai_ref[cai_ref]

    cdef long codon_tot = 0
    cdef int valid_stops = 0
    cdef np.ndarray[dtype=long, ndim=1, mode="c"] ncod = np.zeros([65], dtype=c_long)
    cdef np.ndarray[dtype=long, ndim=1, mode="c"] naa = np.zeros([23], dtype=c_long)
    codonwlib.codon_usage_tot(seq, &codon_tot, &valid_stops, &ncod[0], &naa[0], &ref_code)

    cdef np.ndarray[dtype=int, ndim=1, mode="c"] dds = np.zeros([65], dtype=c_int)
    codonwlib.how_synon(&dds[0], &ref_code)

    cdef double cai = 0;
    codonwlib.cai(&ncod[0], &cai, &dds[0], &ref_cai, &ref_code)

    return cai


"""
def raau(seq):
    int ret = codonwlib.raau_usage(long nnaa[], double raau[])

def silent_base_usage(seq, code):
    int ret = codonwlib.base_sil_us(long *nncod, long *nnaa, double base_sil[], int *ds, int *da, GENETIC_CODE_STRUCT *pcu)

def cbi(seq, code):
    int ret = codonwlib.cbi(long *nncod, long *nnaa, float *fcbi, int *ds, int *da, GENETIC_CODE_STRUCT *pcu, FOP_STRUCT *pcbi)

def fop(seq, code):
    int ret = codonwlib.fop(long *nncod, float *ffop, int *ds, bool factor_in_rare, GENETIC_CODE_STRUCT *pcu, FOP_STRUCT *pfop)

def enc(seq, code):
    int ret = codonwlib.enc(long *nncod, long *nnaa, float *enc_tot, int *da, GENETIC_CODE_STRUCT *pcu)

def gc(seq):
    int ret = codonwlib.gc(int *ds, long bases[5], long base_tot[5], long base_1[5], long base_2[5], long base_3[5], long *tot_s, long *totalaa, GENETIC_CODE_STRUCT *pcu)

def dinuc(seq):
    int ret = codonwlib.dinuc(long dinuc_tot[4])
    
def hydro(seq):
    int ret = codonwlib.hydro(long *nnaa, float *hydro, AMINO_PROP_STRUCT *pap)

def aromo(seq):
    int ret = codonwlib.aromo(long *nnaa, float *aromo, AMINO_PROP_STRUCT *pap)

"""