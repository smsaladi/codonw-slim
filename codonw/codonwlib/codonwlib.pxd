"""
Exports to be used in the Python interface
"""

from libcpp cimport bool

cdef extern from "include/codonW.h":
    ctypedef struct GENETIC_CODE_STRUCT:
        char *des
        char *typ
        int ca[65]

    ctypedef struct FOP_STRUCT:
        char *des
        char *ref
        char fop_cod[65]

    ctypedef struct CAI_STRUCT:
        char *des
        char *ref
        float cai_val[65]

    ctypedef struct AMINO_STRUCT:
        char *aa1[22]
        char *aa3[22]
        char *cod[65]

    ctypedef struct AMINO_PROP_STRUCT:
        float hydro[22]
        int aromo[22]

    ctypedef struct REF_STRUCT:
        GENETIC_CODE_STRUCT *cu
        FOP_STRUCT *fop
        CAI_STRUCT *cai

        AMINO_STRUCT *amino_acids
        AMINO_PROP_STRUCT *amino_prop

    REF_STRUCT Z_ref

    int ident_codon(char *codon)
    int how_synon(int dds[], GENETIC_CODE_STRUCT *pcu)
    int how_synon_aa(int dda[], GENETIC_CODE_STRUCT *pcu)

    int codon_usage_tot(char *seq, long *codon_tot, int *valid_stops, long ncod[], long naa[], GENETIC_CODE_STRUCT *pcu)
    int rscu_usage(long *nncod, long *nnaa, float rscu[], int *ds, GENETIC_CODE_STRUCT *pcu)
    int raau_usage(long nnaa[], double raau[])
    int base_sil_us(long *nncod, long *nnaa, double base_sil[], int *ds, int *da, GENETIC_CODE_STRUCT *pcu)
    int cai(long *nncod, double *sigma, int *ds, CAI_STRUCT *pcai, GENETIC_CODE_STRUCT *pcu)
    int cbi(long *nncod, long *nnaa, float *fcbi, int *ds, int *da, GENETIC_CODE_STRUCT *pcu, FOP_STRUCT *pcbi)
    int fop(long *nncod, float *ffop, int *ds, bool factor_in_rare, GENETIC_CODE_STRUCT *pcu, FOP_STRUCT *pfop)
    int enc(long *nncod, long *nnaa, float *enc_tot, int *da, GENETIC_CODE_STRUCT *pcu)
    int gc(int *ds, long *ncod, long bases[5], long base_tot[5], long base_1[5], long base_2[5], long base_3[5], long *tot_s, long *totalaa, GENETIC_CODE_STRUCT *pcu)
    int dinuc(long dinuc_tot[4])
    int hydro(long *nnaa, float *hydro, AMINO_PROP_STRUCT *pap)
    int aromo(long *nnaa, float *aromo, AMINO_PROP_STRUCT *pap)
