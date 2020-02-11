"""

CodonW C functions are wrapped into methods of the `CodonSeq` class here.
The method documentation is pulled from `README_indicies` as well as
from the C code.

"""

# cython: c_string_type=str, c_string_encoding=ascii

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

"""
Would be great to expose these internals as pd.Series...
    FOP_STRUCT *fop_ref
    CAI_STRUCT *cai_ref
    AMINO_PROP_STRUCT *amino_prop
"""

def get_reference_code(idx):
    cseq = CodonSeq("ATG", idx)
    return cseq.genetic_code

cdef class CodonSeq:
    # Use memory view to arrays
    # https://suzyahyah.github.io/cython/programming/2018/12/01/Gotchas-in-Cython.html
    cdef public codonwlib.GENETIC_CODE_STRUCT ref_code
    cdef public int[::1] dds
    cdef public int[::1] dda

    cdef public object seq
    cdef public long codon_tot
    cdef public int valid_stops
    cdef public long[::1] ncod
    cdef public long[::1] naa

    def __init__(self, object seq, genetic_code=0):
        """Initializes an object of class CodonSeq

        `seq`: the nucleotide sequence to be analyzed/for which metrics are desired

        `genetic_code`: the genetic code to be used
            0. Universal Genetic code [default]
            1. Vertebrate Mitochondrial code
            2. Yeast Mitochondrial code
            3. Filamentous fungi Mitochondrial code
            4. Insects and Plathyhelminthes Mitochondrial code
            5. Nuclear code of Cilitia
            6. Nuclear code of Euplotes
            7. Mitochondrial code of Echinoderms

        """
        self.dds = np.zeros([65], dtype=c_int)
        self.dda = np.zeros([23], dtype=c_int)

        if isinstance(genetic_code, int):
            self.ref_code = codonwlib.cu_ref[genetic_code]
        else:
            self.genetic_code = genetic_code

        codonwlib.how_synon(&self.dds[0], &self.ref_code)
        codonwlib.how_synon_aa(&self.dda[0], &self.ref_code)

        self.codon_tot = 0
        self.valid_stops = 0
        self.ncod = np.zeros([65], dtype=c_long)
        self.naa = np.zeros([22], dtype=c_long)

        self.seq = seq.encode()
        codonwlib.codon_usage_tot(<char *>self.seq,
            &self.codon_tot, &self.valid_stops, &self.ncod[0], &self.naa[0], &self.ref_code)
        
        return

    # Read/Set genetic code through pd.Series
    @property
    def genetic_code(self):
        cdef np.ndarray[dtype=object, ndim=1, mode="c"] aa = np.array(ref_aa1, dtype=object)
        return pd.Series(aa[self.ref_code.ca], index=ref_codons,
                         name=self.ref_code.des.decode('UTF-8'))

    @genetic_code.setter
    def genetic_code(self, ser):
        if 'UNK' not in ser:
            ser['UNK'] = 'X'
        
        # place in required order and map amino acids letter to code
        ser = ser[ref_codons]
        aa_to_idx = pd.Series(np.arange(len(ref_aa1)), index=ref_aa1)

        self.ref_code = codonwlib.GENETIC_CODE_STRUCT(b"", b"")
        self.ref_code.ca = aa_to_idx[ser].values
        return


    cpdef double cai(self, int cai_ref=0):
        """Calculates Codon Adaptation Index

        `cai_ref`: The relative adaptiveness of codon
            0. Escherichia coli - No reference [default]
            1. Bacillus subtilis - No reference
            2. Saccharomyces cerevisiae - Sharp and Cowe (1991) Yeast 7:657-678

        If you'd like to have user-provided reference values, please implement this
        functionality and make a pull-request.


        CAI is a measurement of the relative adaptiveness of the codon usage of a
        gene towards the codon usage of highly expressed genes. The relative
        adaptiveness (w) of each codon is the ratio of the usage of each codon, to
        that of the most abundant codon for the same amino acid.
        
        The CAI index is defined as the geometric mean of these relative
        adaptiveness values. Non-synonymous codons and termination codons (dependent
        on genetic code) are excluded.

        To prevent a codon absent from the reference set but present in other genes
        from having a relative adaptiveness value of zero, which would cause CAI to
        evaluate to zero for any genes which used that codon; it was suggested that
        absent codons should be assigned a frequency of 0.5 (Sharp and Li 1987).
        
        An alternative suggestion was that such codons should be adjusted to
        0.01, where otherwise it would be less than this value
        [(Bulmer 1988)](https://doi.org/10.1046/j.1420-9101.1988.1010015.x).
        
        The CAI is calculated as using a natural log summation. To prevent a codon having
        a relative adaptiveness value of zero, which could result in a CAI of zero, codons
        with a value of < 0.0001 are adjusted to 0.01.

        [Sharp and Li 1987](https://doi.org/10.1093/nar/15.3.1281)
        """
        cdef codonwlib.CAI_STRUCT ref_cai = codonwlib.cai_ref[cai_ref]
        cdef double cai_val = 0
        cdef int ret = codonwlib.cai(&self.ncod[0], &cai_val, &self.dds[0], &ref_cai, &self.ref_code)
        return cai_val

    cpdef float cbi(self, int cai_ref=0):
        """Calculate codon bias index

        `cai_ref`: The relative adaptiveness of codon
            0. Escherichia coli - No reference [default]
            1. Bacillus subtilis - No reference
            2. Saccharomyces cerevisiae - Sharp and Cowe (1991) Yeast 7:657-678

        If you'd like to have user-provided reference values, please implement this
        functionality and make a pull-request.


        Codon bias index is another measure of directional codon bias, it measures
        the extent to which a gene uses a subset of optimal codons. CBI is similar
        to Fop, with expected usage used as a scaling factor.

            CBI = (Nopt-Nran)/(Nopt-Nran)
            
        where Nopt = number of optimal codons
              Ntot = number of synonymous codons
              Nran = expected number of optimal codons for randomly assigned codons 
        
        In a gene with extreme codon bias, CBI will equal 1.0. In a gene with random
        codon usage CBI will equal 0.0. It is possible for the number of optimal
        codons to be less than expected by random change (i.e. Nopt < Nran). This
        results in a negative value for CBI.

        [Bennetzen and Hall 1982](https://europepmc.org/article/MED/7037777)
        """
        cdef float cbi_val
        cdef int ret = codonwlib.cbi(&self.ncod[0], &self.naa[0], &cbi_val, \
            &self.dds[0], &self.dda[0], &self.ref_code, &codonwlib.fop_ref[cai_ref])
        return cbi_val

    cpdef float fop(self, bool factor_in_rare=False, int fop_ref=0):
        """Calculate fraction of optimal codons

        `fop_ref`:
            0. Escherichia coli
                - [Ikemura 1985](https://doi.org/10.1093/oxfordjournals.molbev.a040335)
                  (updated by Irish National Centre for BioInformatics 1991)
            1. Bacillus subtilis
                - [Sharp, et al. 1990](https://doi.org/10.1016/B978-0-12-274162-3.50013-X)
            2. Dictyostelium discoideum
                - [Sharp & Devine 1989](https://doi.org/10.1093/nar/17.13.5029)
            3. Aspergillus nidulans
                - [Lloyd & Sharp 1991](https://doi.org/10.1007/bf00290679)
            4. Saccharomyces cerevisiae
                - [Sharp & Cowe 1991](https://doi.org/10.1002/yea.320070702)
            5. Drosophila melanogaster
                - [Shields, et al. 1988](https://doi.org/10.1093/oxfordjournals.molbev.a040525)
            6. Caenorhabditis elegans
                - [Stenico, Lloyd, & Sharp 1994](https://doi.org/10.1093/nar/22.13.2437)
            7. Neurospora crassa
                - Lloyd & Sharp 1993 (Citation cannot be found)
    
            If you'd like to have user-provided reference values, please implement this
            functionality and make a pull-request.

        `factor_in_rare`: 
            If non-optimal codons are identified in the set of optimal codons selected, use
            the following formulation: Fop = (opt-rare)/total


        This index, is the ratio of optimal codons to synonymous codons:

            Fop = opt/total

        Fop values are always between 0 (where no optimal codons are used)
        and 1 (where only optimal codons are used). When `factor_in_rare` is True,
        negative values are adjusted to zero.

        [Ikemura 1981](https://doi.org/10.1016/0022-2836(81)90003-6)
        """
        cdef float fop_val
        cdef int ret = codonwlib.fop(&self.ncod[0], &fop_val, \
            &self.dds[0], factor_in_rare, &self.ref_code, &codonwlib.fop_ref[fop_ref])
        return fop_val

    cpdef float enc(self):
        """Calculate effective number of codons

        This index is a simple measure of overall codon bias and is analogous to the
        effective number of alleles measure used in population genetics. Knowledge
        of the optimal codons or a reference set of highly expressed genes is
        unnecessary. Initially the homozygosity for each amino acid is estimated
        from the squared codon frequencies.

        The reported value of Nc is always between 20 (when only one codon is
        effectively used for each amino acid) and 61 (when codons are used randomly).
        If the calculated Nc is greater than 61 (because codon usage is more evenly
        distributed than expected), it is adjusted to 61.

        If amino acids are rare or missing, adjustments must be made. When
        there are no amino acids in a synonymous family, Nc is not calculated
        as the gene is either too short or has extremely skewed amino acid
        usage. An exception to this is made for genetic codes where isoleucine is the
        only 3-fold synonymous amino acid, and is not used in the protein gene.

       [Wright 1990](https://doi.org/10.1016/0378-1119(90)90491-9)
        """
        cdef float enc_val
        cdef int ret = codonwlib.enc(&self.ncod[0], &self.naa[0], &enc_val, \
            &self.dda[0], &self.ref_code)
        return enc_val

    cpdef float hydropathy(self):
        """Calculate mean hydropathy

        The general average hydropathicity or (GRAVY) score, for the hypothetical
        translated gene product. It is the arithmetic mean hydropathy values
        assigned to each amino acid in a protein.

        [Kyte & Doolittle 1982](https://doi.org/10.1016/0022-2836(82)90515-0)
        """
        cdef float hydro_val
        cdef int ret = codonwlib.hydro(&self.naa[0], &hydro_val, \
            <float (*)>codonwlib.amino_prop.hydro)
        return hydro_val

    cpdef float aromaticity(self):
        """Calculate mean aromaticity

        The frequency of aromatic amino acids (Phe, Tyr, Trp) in the hypothetical
        translated gene product.
        """
        cdef float aromo_val
        cdef int ret = codonwlib.aromo(&self.naa[0], &aromo_val, \
            <int (*)>codonwlib.amino_prop.aromo)
        return aromo_val


    cpdef np.ndarray[dtype=double, ndim=1, mode="c"] silent_base_usage_(self):
        cdef np.ndarray[dtype=double, ndim=1, mode="c"] base_sil_vals = np.zeros([4], dtype=c_double)
        cdef int ret = codonwlib.base_sil_us(&self.ncod[0], &self.naa[0], &base_sil_vals[0],
                                        &self.dds[0], &self.dda[0], &self.ref_code)
        return base_sil_vals

    def silent_base_usage(self):
        """Calculate silent base usage

        Calculates the base composition at silent sites normalised by the possible
        usage at that silent site with changing the amino acid composition.
        For example, the index A3s is the frequency of codons have an A at their
        synonymous third position, relative to the amino acids that could have a
        synonym with A in the synonymous third codon position.

        It is inspired by GC3s but is more complicated to calculate as not every
        AA can use any base at the third position. When calculating GC3s each
        synonymous amino acid has at least one synonym with G or C in the third
        position. Two or three fold synonymous amino acids do not have an equal
        choice between bases in the synonymous third position.
        It's correlated with GC3s but not directly comparable. 
        """
        return pd.Series(self.silent_base_usage_(), index=['G3s', 'C3s', 'A3s', 'T3s'])
        

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
        cdef np.ndarray[dtype=double, ndim=1, mode="c"] raau_vals = np.zeros([22], dtype=c_double)
        cdef int ret = codonwlib.raau_usage(&self.naa[0], &raau_vals[0])
        return raau_vals

    def raau(self):
        """Calculate Relative Amino Acid Usage
        """
        return pd.Series(self._raau(), index=ref_aa1)


    cpdef np.ndarray[dtype=long, ndim=2, mode="c"] _bases(self):
        cdef long tot_s
        cdef long totalaa
        cdef np.ndarray[dtype=long, ndim=2, mode="c"] bases = np.zeros([5, 5], dtype=c_long)
        cdef np.ndarray[dtype=double, ndim=1, mode="c"] metrics = np.zeros([18], dtype=c_double)

        cdef int ret = codonwlib.gc(&self.dds[0], &self.ncod[0],
            &bases[4, 0], &bases[3, 0], &bases[0, 0], &bases[1, 0], &bases[2, 0],
            &tot_s, &totalaa, &metrics[0], &self.ref_code)

        return bases[0:6, 1:5]

    def bases(self):
        """Calculates base composition (overall and by position)
        """
        v = pd.DataFrame(self._bases(),
            columns=['T', 'C', 'A', 'G'],
            index=['1', '2', '3', 'all', 'syn'])
        return v


    cpdef np.ndarray[dtype=double, ndim=1, mode="c"] _gc(self):
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

    def bases2(self):
        """Calculates additional metrics related to nucleotide base composition

        These metrics include the following and are returned as a pd.Series
            * base composition in all frames
            * length of gene
            * Number of synonymous codons
            * G+C content (overall and by codon position)
            * G+C content of synonymous codons at the 3rd position
            * G+C content of non-synonymous codons at the 3rd position
            * Number of synonymous codons
            * Number of amino acids
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
        cdef np.ndarray[dtype=long, ndim=2, mode="c"] dinuc_frames = np.zeros([4, 16], dtype=c_long)
        cdef np.ndarray[dtype=long, ndim=1, mode="c"] dinuc_tot = np.zeros([4], dtype=c_long)
        cdef int fram = 0

        cdef int ret = codonwlib.dinuc_count(<char *>self.seq,
            <long (*)[16]>&dinuc_frames[0, 0], &dinuc_tot[0], &fram)

        dinuc_frames[3, :] = np.sum(dinuc_frames, axis=0)
            
        if pct:
            return dinuc_frames / np.reshape(dinuc_tot, [4, 1])
        
        return dinuc_frames.astype(np.double)

    def dinuc(self, pct=True):
        """Calculate Dinucleotide Usage
        
        `pct`:
            If True, report percentages

        The frequency of all 16 dinucleotides, in total, and across
        all three possible reading frames, i.e. `1:2`, `2:3`, `3:1`.
        """
        
        cdef np.ndarray[dtype=double, ndim=2, mode="c"] frames = self._dinuc(pct)

        if pct:
            def convert(x): return x
        else:
            def convert(x): return x.astype(long)

        v = pd.DataFrame(convert(frames),
            columns=['TT', 'TC', 'TA', 'TG',
                   'CT', 'CC', 'CA', 'CG',
                   'AT', 'AC', 'AA', 'AG',
                   'GT', 'GC', 'GA', 'GG'],
            index=['1:2', '2:3', '3:1', 'all'])
        
        return v
