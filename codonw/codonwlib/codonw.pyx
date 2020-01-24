"""
Wrappers around CodonW C functions
"""

from libcpp cimport bool
from cython.operator cimport dereference
from ctypes import c_int, c_long, c_float, c_double

import numpy as np
cimport numpy as np
np.import_array()
import pandas as pd

cimport codonwlib

def convert_char(arr):
    return [c.decode('UTF-8') for c in arr]

ref_codons = convert_char(codonwlib.amino_acids.cod)
ref_aa1 = convert_char(codonwlib.amino_acids.aa1)
ref_aa3 = convert_char(codonwlib.amino_acids.aa3)

aa1_aa3 = pd.Series(index=ref_aa1, data=ref_aa3)
aa3_aa1 = pd.Series(index=ref_aa3, data=ref_aa1)

cdef class CodonSeq:
    # Use memory view to arrays
    # https://suzyahyah.github.io/cython/programming/2018/12/01/Gotchas-in-Cython.html
    cdef public genetic_code
    cdef public codonwlib.GENETIC_CODE_STRUCT ref_code
    cdef public int[::1] dds
    cdef public int[::1] dda

    cdef public char* seq
    cdef public long codon_tot
    cdef public int valid_stops
    cdef public long[::1] ncod
    cdef public long[::1] naa

    def __init__(self, char* seq, int genetic_code=0):
        self.dds = np.zeros([65], dtype=c_int)
        self.dda = np.zeros([23], dtype=c_int)
        self.genetic_code = genetic_code
        self.ref_code = codonwlib.cu_ref[self.genetic_code]

        codonwlib.how_synon(&self.dds[0], &self.ref_code)
        codonwlib.how_synon_aa(&self.dda[0], &self.ref_code)

        self.codon_tot = 0
        self.valid_stops = 0
        self.ncod = np.zeros([65], dtype=c_long)
        self.naa = np.zeros([22], dtype=c_long)

        self.seq = seq
        codonwlib.codon_usage_tot(self.seq,
            &self.codon_tot, &self.valid_stops, &self.ncod[0], &self.naa[0], &self.ref_code)
        
        return

    cpdef double cai(self, int cai_ref=0):
        """Calculates Codon Adaptation Index
        """
        cdef codonwlib.CAI_STRUCT ref_cai = codonwlib.cai_ref[cai_ref]
        cdef double cai_val = 0
        cdef int ret = codonwlib.cai(&self.ncod[0], &cai_val, &self.dds[0], &ref_cai, &self.ref_code)
        return cai_val

    cpdef float cbi(self, int cai_ref=0):
        """Calculate codon bias index
        """
        cdef float cbi_val
        cdef int ret = codonwlib.cbi(&self.ncod[0], &self.naa[0], &cbi_val, \
            &self.dds[0], &self.dda[0], &self.ref_code, &codonwlib.fop_ref[cai_ref])
        return cbi_val

    cpdef float fop(self, bool factor_in_rare=False, int fop_ref=0):
        """Calculate fraction of optimal codons
        """
        cdef float fop_val
        cdef int ret = codonwlib.fop(&self.ncod[0], &fop_val, \
            &self.dds[0], factor_in_rare, &self.ref_code, &codonwlib.fop_ref[fop_ref])
        return fop_val

    cpdef float enc(self):
        """Calculate effective number of codons (E_nc)
        """
        cdef float enc_val
        cdef int ret = codonwlib.enc(&self.ncod[0], &self.naa[0], &enc_val, \
            &self.dda[0], &self.ref_code)
        return enc_val

    cpdef float hydropathy(self):
        """Calculate mean hydropathy
        """
        cdef float hydro_val
        cdef int ret = codonwlib.hydro(&self.naa[0], &hydro_val, \
            <float (*)>codonwlib.amino_prop.hydro)
        return hydro_val

    cpdef float aromaticity(self):
        """Calculate mean aromaticity
        """
        cdef float aromo_val
        cdef int ret = codonwlib.aromo(&self.naa[0], &aromo_val, \
            <int (*)>codonwlib.amino_prop.aromo)
        return aromo_val

    cpdef np.ndarray[dtype=double, ndim=1, mode="c"] silent_base_usage(self):
        """Calculate silent base usage
        """
        cdef np.ndarray[dtype=double, ndim=1, mode="c"] base_sil_vals = np.zeros([4], dtype=c_double)
        cdef int ret = codonwlib.base_sil_us(&self.ncod[0], &self.naa[0], &base_sil_vals[0],
                                        &self.dds[0], &self.dda[0], &self.ref_code)
        return base_sil_vals

    def codon_usage(self):
        """Codon tabulation
        """
        return pd.Series(self.ncod[1:65], index=ref_codons[1:65])
        
    def aa_usage(self):
        """Amino acid tabulation
        """
        return pd.Series(self.naa, index=ref_aa1)

    cpdef np.ndarray[dtype=float, ndim=1, mode="c"] _rscu(self):
        """Calculate Relative Synonymous Codon Usage
        """
        cdef np.ndarray[dtype=float, ndim=1, mode="c"] rscu_vals = np.zeros([65], dtype=c_float)
        cdef int ret = codonwlib.rscu_usage(&self.ncod[0], &self.naa[0], &rscu_vals[0], &self.dds[0], &self.ref_code)
        return rscu_vals

    def rscu(self):
        """Calculate Relative Synonymous Codon Usage
        """
        return pd.Series(self._rscu()[1:65], index=ref_codons[1:65])

    cpdef np.ndarray[dtype=double, ndim=1, mode="c"] _raau(self):
        """Calculate Relative Amino Acid Usage
        """
        cdef np.ndarray[dtype=double, ndim=1, mode="c"] raau_vals = np.zeros([22], dtype=c_double)
        cdef int ret = codonwlib.raau_usage(&self.naa[0], &raau_vals[0])
        return raau_vals

    def raau(self):
        """Calculate Relative Amino Acid Usage
        """
        return pd.Series(self._raau(), index=ref_aa1)

    cpdef np.ndarray[dtype=long, ndim=2, mode="c"] _bases(self):
        """Calculates base composition
        """
        cdef long tot_s
        cdef long totalaa
        cdef np.ndarray[dtype=long, ndim=2, mode="c"] bases = np.zeros([5, 5], dtype=c_long)
        cdef np.ndarray[dtype=double, ndim=1, mode="c"] metrics = np.zeros([18], dtype=c_double)

        cdef int ret = codonwlib.gc(&self.dds[0], &self.ncod[0],
            &bases[4, 0], &bases[3, 0], &bases[0, 0], &bases[1, 0], &bases[2, 0],
            &tot_s, &totalaa, &metrics[0], &self.ref_code)

        return bases[0:6, 1:5]

    def bases(self):
        """Calculates base composition
        """
        v = pd.DataFrame(self._bases(),
            columns=['T', 'C', 'A', 'G'],
            index=['1', '2', '3', 'all', 'syn'])
        return v

    cpdef np.ndarray[dtype=double, ndim=1, mode="c"] _gc(self):
        """Calculates various %GC-related metrics
        """
        cdef long tot_s
        cdef long totalaa
        cdef np.ndarray[dtype=long, ndim=2, mode="c"] bases = np.zeros([5, 5], dtype=c_long)
        cdef np.ndarray[dtype=double, ndim=1, mode="c"] metrics = np.zeros([20], dtype=c_double)

        cdef int ret = codonwlib.gc(&self.dds[0], &self.ncod[0],
            &bases[4, 0], &bases[3, 0], &bases[0, 0], &bases[1, 0], &bases[2, 0],
            &tot_s, &totalaa, &metrics[2], &self.ref_code)

        metrics[0] = <double>totalaa;
        metrics[1] = <double>tot_s;
        return metrics

    def gc(self):
        """Calculates various %GC-related metrics
        """
        v = pd.Series(self._gc(),
            index=['Len_aa', 'Len_sym',
                   'GC', 'GC3s', 'GCn3s',
                   'GC1', 'GC2', 'GC3',
                   'T1', 'T2', 'T3',
                   'C1', 'C2', 'C3',
                   'A1', 'A2', 'A3',
                   'G1', 'G2', 'G3'])
        return v

    cpdef np.ndarray[dtype=double, ndim=2, mode="c"] _dinuc(self, bool pct):
        """Calculate Dinucleotide Usage
        """
        cdef np.ndarray[dtype=long, ndim=2, mode="c"] dinuc_frames = np.zeros([4, 16], dtype=c_long)
        cdef np.ndarray[dtype=long, ndim=1, mode="c"] dinuc_tot = np.zeros([4], dtype=c_long)
        cdef int fram = 0

        cdef int ret = codonwlib.dinuc_count(self.seq,
            <long (*)[16]>&dinuc_frames[0, 0], &dinuc_tot[0], &fram)

        dinuc_frames[3, :] = np.sum(dinuc_frames, axis=0)

        cdef np.ndarray[dtype=double, ndim=2, mode="c"] dinuc_frames_out = \
            np.zeros([4, 16], dtype=c_double)

        if pct:
            dinuc_frames_out = dinuc_frames / np.reshape(dinuc_tot, [4, 1])
        else:
            dinuc_frames_out = dinuc_frames.astype(np.double)

        return dinuc_frames_out

    def dinuc(self, pct=True):
        """Calculate Dinucleotide Usage
        If `pct`, report percentages
        """
        
        cdef np.ndarray[dtype=double, ndim=2, mode="c"] frames = self._dinuc(pct)

        if pct:
            def conv(x): return x
        else:
            def conv(x): return x.astype(long)

        v = pd.DataFrame(conv(frames),
            columns=['TT', 'TC', 'TA', 'TG',
                   'CT', 'CC', 'CA', 'CG',
                   'AT', 'AC', 'AA', 'AG',
                   'GT', 'GC', 'GA', 'GG'],
            index=['1:2', '2:3', '3:1', 'all'])
        
        return v

"""
Want to expose codonwlib.
    GENETIC_CODE_STRUCT *cu_ref
    FOP_STRUCT *fop_ref
    CAI_STRUCT *cai_ref
    AMINO_STRUCT *amino_acids
    AMINO_PROP_STRUCT *amino_prop
"""