"""
Wrappers around CodonW C functions
"""

from libcpp cimport bool
from cython.operator cimport dereference
from ctypes import c_int, c_long, c_float, c_double

import numpy as np
cimport numpy as np
np.import_array()

cimport codonwlib

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
        self.naa = np.zeros([23], dtype=c_long)

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
        cdef int ret = codonwlib.cbi(&self.ncod[0], &self.naa[0], &cbi_val, &self.dds[0], &self.dda[0], &self.ref_code, &codonwlib.fop_ref[cai_ref])
        return cbi_val

    cpdef float fop(self, bool factor_in_rare=False, int fop_ref=0):
        """Calculate fraction of optimal codons
        """
        cdef float fop_val
        cdef int ret = codonwlib.fop(&self.ncod[0], &fop_val, &self.dds[0], factor_in_rare, &self.ref_code, &codonwlib.fop_ref[fop_ref])
        return fop_val

    cpdef float enc(self):
        """Calculate effective number of codons (E_nc)
        """
        cdef float enc_val
        cdef int ret = codonwlib.enc(&self.ncod[0], &self.naa[0], &enc_val, &self.dda[0], &self.ref_code)
        return enc_val

    cpdef float hydro(self):
        """Calculate mean hydropathy
        """
        cdef float hydro_val
        cdef int ret = codonwlib.hydro(&self.naa[0], &hydro_val, <float (*)>codonwlib.amino_prop.hydro)
        return hydro_val

    cpdef float aromo(self):
        """Calculate mean aromaticity
        """
        cdef float aromo_val
        cdef int ret = codonwlib.aromo(&self.naa[0], &aromo_val, <int (*)>codonwlib.amino_prop.aromo)
        return aromo_val

    cpdef np.ndarray[dtype=double, ndim=1, mode="c"] silent_base_usage(self):
        """Calculate silent base usage
        """
        cdef np.ndarray[dtype=double, ndim=1, mode="c"] base_sil_vals = np.zeros([4], dtype=c_double)
        cdef int ret = codonwlib.base_sil_us(&self.ncod[0], &self.naa[0], &base_sil_vals[0],
                                        &self.dds[0], &self.dda[0], &self.ref_code)
        return base_sil_vals

    cpdef np.ndarray[dtype=float, ndim=1, mode="c"] rscu(self):
        """Calculate Relative Synonymous Codon Usage
        """
        cdef np.ndarray[dtype=float, ndim=1, mode="c"] rscu_vals = np.zeros([65], dtype=c_float)
        cdef int ret = codonwlib.rscu_usage(&self.ncod[0], &self.naa[0], &rscu_vals[0], &self.dds[0], &self.ref_code)
        return rscu_vals

    cpdef np.ndarray[dtype=double, ndim=1, mode="c"] raau(self):
        """Calculate Relative Amino Acid Usage
        """
        cdef np.ndarray[dtype=double, ndim=1, mode="c"] raau_vals = np.zeros([65], dtype=c_double)
        cdef int ret = codonwlib.raau_usage(&self.naa[0], &raau_vals[0])
        return raau_vals

    cpdef gc(self):
        """Calculate various %GC-elated metrics
        """
        cdef long tot_s
        cdef long totalaa
        cdef np.ndarray[dtype=long, ndim=1, mode="c"] bases = np.zeros([5], dtype=c_long)
        cdef np.ndarray[dtype=long, ndim=1, mode="c"] base_tot = np.zeros([5], dtype=c_long)
        cdef np.ndarray[dtype=long, ndim=1, mode="c"] base_1 = np.zeros([5], dtype=c_long)
        cdef np.ndarray[dtype=long, ndim=1, mode="c"] base_2 = np.zeros([5], dtype=c_long)
        cdef np.ndarray[dtype=long, ndim=1, mode="c"] base_3 = np.zeros([5], dtype=c_long)

        cdef int ret = codonwlib.gc(&self.dds[0], &self.ncod[0],
            &bases[0], &base_tot[0], &base_1[0], &base_2[0], &base_3[0],
            &tot_s, &totalaa, &self.ref_code)
        
        return (bases, base_tot, base_1, base_2, base_3)

    cpdef np.ndarray[dtype=long, ndim=1, mode="c"] dinuc(self):
        """Calculate Dinucleotide Usage
        """
        cdef np.ndarray[dtype=long, ndim=2, mode="c"] din = np.zeros([3, 16], dtype=c_long)
        cdef np.ndarray[dtype=long, ndim=1, mode="c"] dinuc_tot = np.zeros([4], dtype=c_long)
        cdef int fram = 0;

        cdef int ret = codonwlib.dinuc_count(self.seq, <long (*)[16]>&din[0, 0], &dinuc_tot[0], &fram)

        return dinuc_tot

"""
Want to expose codonwlib.
    GENETIC_CODE_STRUCT *cu_ref
    FOP_STRUCT *fop_ref
    CAI_STRUCT *cai_ref
    AMINO_STRUCT *amino_acids
    AMINO_PROP_STRUCT *amino_prop
"""