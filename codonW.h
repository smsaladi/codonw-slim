/**************************************************************************/
/* CodonW codon usage analysis package                                    */
/* Copyright (C) 2005            John F. Peden                            */
/* This program is free software; you can redistribute                    */
/* it and/or modify it under the terms of the GNU General Public License  */
/* as published by the Free Software Foundation; version 2 of the         */
/* License,                                                               */
/*                                                                        */
/* This program is distributed in the hope that it will be useful, but    */
/* WITHOUT ANY WARRANTY; without even the implied warranty of             */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the           */
/* GNU General Public License for more details.                           */
/* You should have received a copy of the GNU General Public License along*/
/* with this program; if not, write to the Free Software Foundation, Inc.,*/
/* 675 Mass Ave, Cambridge, MA 02139, USA.                                */
/*                                                                        */
/*                                                                        */
/* The author can be contacted by email (jfp#hanson-codonw@yahoo.com Anti-*/
/* Spam please change the # in my email to an _)                          */
/*                                                                        */
/* For the latest version and information see                             */
/* http://codonw.sourceforge.net 					  */
/**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <stdbool.h>

#define GARG_EXACT 0x800             /* used in function gargs  */
#define GARG_NEXT 0x1000             /* used in function gargs  */
#define GARG_THERE 0x2000            /* used in function gargs  */
#define GARG_SUBSQ 0x4000            /* used in function gargs  */
#define MAX_ARGS 100                 /* used in function gargs  */
#define MAX_FILENAME_LEN 90
#define MAX_MESSAGE_LEN 300

#define NUM_GENETIC_CODES 8
#define NUM_FOP_SPECIES 8
#define NUM_CAI_SPECIES 3

/* define the structures used within codonW                               */
typedef struct
{
  char *des;
  char *typ;
  int ca[65];
} GENETIC_CODE_STRUCT; /* genetic code information */

typedef struct
{
  char *aa1[22]; /* 1 letter AA code         */
  char *aa3[22]; /* 3 letter AA code         */
  char *cod[65]; /* 3 letter name of codons  */
} AMINO_STRUCT;

typedef struct
{
  float hydro[22]; /* hydropathicity values    */
  int aromo[22];   /* aromaticity values       */
} AMINO_PROP_STRUCT;

typedef struct
{
  char *des;        /* store a description      */
  char *ref;        /*       a reference        */
  char fop_cod[65]; /* the optimal codons       */
} FOP_STRUCT;

typedef struct
{
  char *des;         /* store a description      */
  char *ref;         /*       a reference        */
  float cai_val[65]; /* the CAI w values         */
} CAI_STRUCT;

typedef struct
{
  char bulk;    /* used to ident blk output */
  char totals;  /* concatenate genes ?      */
  char warn;    /* show sequence warning    */

  char fop;      /* calc index fop           */
  char cai;      /* calc index CAI           */
  char cbi;      /* calc index CBI           */
  char bases;    /* calc base composition    */
  char gc3s;     /* calc gc at sil.3rd base  */
  char gc;       /* calc gc                  */
  char enc;      /* calc enc                 */
  char sil_base; /* calc silent base compo   */
  char L_sym;    /* No of synonymous codons  */
  char L_aa;     /* No of amino acids        */
  char hyd;      /* calc hydropathicity      */
  char aro;      /* calc aromaticity         */

  char separator; /* column separator         */

  char code;   /* which genetic code       */
  char f_type; /* which predefined fop val */
  char c_type; /* which predefined CAI val */

  FILE *inputfile;   /* input file               */
  FILE *outputfile;  /* .out file                */
  FILE *tidyoutfile; /* .blk file                */
  FILE *cuout;       /* codon usage output       */
  FILE *fopfile;     /* fop input values         */
  FILE *caifile;     /* cai input values         */
  FILE *cbifile;     /* cbi input values         */
  FILE *logfile;     /* log file name            */
  FILE *my_err;      /* pointer for err stream   */

  AMINO_STRUCT *paa;        /* Settings for indicies */
  GENETIC_CODE_STRUCT *pcu;  
  FOP_STRUCT *pfop;
  FOP_STRUCT *pcbi;
  CAI_STRUCT *pcai;
  AMINO_PROP_STRUCT *pap;
  int *da;
  int *ds;
} MENU_STRUCT;

typedef struct {
  GENETIC_CODE_STRUCT *cu;
  FOP_STRUCT *fop;
  CAI_STRUCT *cai;

  AMINO_STRUCT *amino_acids;
  AMINO_PROP_STRUCT *amino_prop;
} REF_STRUCT;

extern REF_STRUCT Z_ref;
extern MENU_STRUCT Z_menu;

/****************** Function type declarations *****************************/

// defined in commline.c
int proc_comm_line(int *argc, char ***arg_list, MENU_STRUCT *pm);

// defined in codons.c
int tidy(FILE *finput, FILE *foutput, FILE *fblkout);
int my_exit(int exit_value, char *message);
FILE *open_file(char *filename, char *mode);

// defined in codon_us.c
int clean_up(long *ncod, long *naa, long din[3][16], int *fram, int *valid_stops);
int initialize_point(char code, char fop_type, char cai_type, MENU_STRUCT *pm, REF_STRUCT *ref);

int count_codons(long* ncod, long *loc_cod_tot);
int dinuc_count(char *seq, long din[3][16], int *fram);

int codon_usage_tot(char *seq, long *codon_tot, int *valid_stops, long ncod[], long naa[], GENETIC_CODE_STRUCT *pcu);
int codon_usage_out(FILE *fblkout, long *ncod, char *info, MENU_STRUCT *pm);
int rscu_usage_out(FILE *fblkout, long *ncod, long *naa, char* title, MENU_STRUCT *pm);
int raau_usage_out(FILE *fblkout, long *naa, char* title, MENU_STRUCT *pm);
int aa_usage_out(FILE *fblkout, long *naa, char* title, MENU_STRUCT *pm);
int cai_out(FILE *foutput, long *ncod, MENU_STRUCT *pm);
int cbi_out(FILE *foutput, long *ncod, long *naa, MENU_STRUCT *pm);
int fop_out(FILE *foutput, long *ncod, MENU_STRUCT *pm);
int hydro_out(FILE *foutput, long *naa, char* title, MENU_STRUCT *pm);
int aromo_out(FILE *foutput, long *naa, char* title, MENU_STRUCT *pm);
int cutab_out(FILE *fblkout, long *nncod, long *nnaa, char* title, MENU_STRUCT *pm);
int dinuc_out(long din[3][16], FILE *fblkout, char *title, char sep);
int enc_out(FILE *foutput, long *ncod, long *naa, MENU_STRUCT *pm);
int gc_out(FILE *foutput, FILE *fblkout, long *ncod, int which, char* title, MENU_STRUCT *pm);
int base_sil_us_out(FILE *foutput, long *ncod, long *naa, MENU_STRUCT *pm);


int rscu_usage(long *nncod, long *nnaa, float rscu[], int *ds, GENETIC_CODE_STRUCT *pcu);
int raau_usage(long nnaa[], double raau[]);
int base_sil_us(long *nncod, long *nnaa, double base_sil[], int *ds, int *da, GENETIC_CODE_STRUCT *pcu);
int cai(long *nncod, double *sigma, int *ds, CAI_STRUCT *pcai, GENETIC_CODE_STRUCT *pcu);
int cbi(long *nncod, long *nnaa, float *fcbi, int *ds, int *da, GENETIC_CODE_STRUCT *pcu, FOP_STRUCT *pcbi);
int fop(long *nncod, float *ffop, int *ds, bool factor_in_rare, GENETIC_CODE_STRUCT *pcu, FOP_STRUCT *pfop);
int enc(long *nncod, long *nnaa, float *enc_tot, int *da, GENETIC_CODE_STRUCT *pcu);
int gc(int *ds, long *ncod, long bases[5], long base_tot[5], long base_1[5], long base_2[5], long base_3[5], long *tot_s, long *totalaa, GENETIC_CODE_STRUCT *pcu);
int dinuc(long din[3][16], long dinuc_tot[4]);
int hydro(long *nnaa, float *hydro, AMINO_PROP_STRUCT *pap);
int aromo(long *nnaa, float *aromo, AMINO_PROP_STRUCT *pap);
