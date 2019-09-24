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
  char level;                     /* either expert or standard*/
  int axis;                       /* how many axis to generate*/
  int rows;                       /* how many genes in dataset*/
  int colm;                       /* how many columns in data */
  int fop_gene;                   /* No of genes to use to ident opt codon*/
  char add_row[MAX_FILENAME_LEN]; /* file with supp sequences */
  float inertia;                  /* total data inertia       */
  char codons[65];                /* codon to be analysed     */
  char amino[22];                 /* amino acids to be COA'ed */
} COA_STRUCT;

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
  char coa;       /* calculate a COA or not ? */

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
  FILE *fcoa_in;
  FILE *fcoa_out;

  AMINO_STRUCT *paa;        /* Settings for indicies */
  GENETIC_CODE_STRUCT *pcu;  
  FOP_STRUCT *pfop;
  FOP_STRUCT *pcbi;
  CAI_STRUCT *pcai;
  COA_STRUCT *pcoa;
  AMINO_PROP_STRUCT *pap;
  int *da;
  int *ds;
} MENU_STRUCT;

typedef struct {
  GENETIC_CODE_STRUCT *cu;
  FOP_STRUCT *fop;
  CAI_STRUCT *cai;

  COA_STRUCT *coa;
  AMINO_STRUCT *amino_acids;
  AMINO_PROP_STRUCT *amino_prop;
} REF_STRUCT;

extern REF_STRUCT Z_ref;
extern MENU_STRUCT Z_menu;

extern char *title;
extern long ncod[65];
extern long naa[23];
extern long din[3][16];
extern long codon_tot;
extern long num_sequence;
extern long num_seq_int_stop;
extern long tot;
extern int last_aa;
extern int valid_stops;
extern int fram;

/****************** Function type declarations *****************************/

// defined in commline.c
int proc_comm_line(int *argc, char ***arg_list, MENU_STRUCT *pm);

// defined in codons.c
int tidy(FILE *finput, FILE *foutput, FILE *fblkout,
         FILE *fcoaout);
int my_exit(int exit_value, char *message);
int fileclose(FILE **file_pointer);
FILE *open_file(char *filename, char *mode);

// defined in codon_us.c
int clean_up(long *ncod, long *naa);
int initialize_point(char code, char fop_type, char cai_type, MENU_STRUCT *pm, REF_STRUCT *ref);
int initialize_coa(COA_STRUCT *pcoa, GENETIC_CODE_STRUCT *pcu, int ds[]);

long codon_error(int last_aa, int valid_stops, char *title,
                     char error_level, MENU_STRUCT *pm);
int dinuc_count(char *seq, long tot);

int codon_usage_tot(char *seq, long *codon_tot, long ncod[], long naa[], MENU_STRUCT *pm);
int codon_usage_out(FILE *fblkout, long *ncod, int last_aa,
                    int valid_stops, char *info, MENU_STRUCT *pm);
int rscu_usage_out(FILE *fblkout, long *ncod, long *naa, MENU_STRUCT *pm);
int raau_usage_out(FILE *fblkout, long *naa, MENU_STRUCT *pm);
char coa_raw_out(FILE *fcoaout, long *ncod, long *naa, char *title, MENU_STRUCT *pm);
int aa_usage_out(FILE *fblkout, long *naa, MENU_STRUCT *pm);
int cai_out(FILE *foutput, long *ncod, MENU_STRUCT *pm);
int cbi_out(FILE *foutput, long *ncod, long *naa, MENU_STRUCT *pm);
int fop_out(FILE *foutput, long *ncod, MENU_STRUCT *pm);
int hydro_out(FILE *foutput, long *naa, MENU_STRUCT *pm);
int aromo_out(FILE *foutput, long *naa, MENU_STRUCT *pm);
int cutab_out(FILE *fblkout, long *nncod, long *nnaa, MENU_STRUCT *pm);
int dinuc_out(FILE *fblkout, char *title, char sep);
float enc_out(FILE *foutput, long *ncod, long *naa, MENU_STRUCT *pm);
void gc_out(FILE *foutput, FILE *fblkout, int which, MENU_STRUCT *pm);
void gen_cusort_fop(int *sortax1, int lig, FILE *fnam, FILE *summ, MENU_STRUCT *pm);
void base_sil_us_out(FILE *foutput, long *ncod, long *naa, MENU_STRUCT *pm);

// defined in coresp.c
void DiagoRC(FILE *summary, MENU_STRUCT *pm, COA_STRUCT *pcoa);
void textbin(char *filein, char *fileout, MENU_STRUCT *pm, COA_STRUCT *pcoa, GENETIC_CODE_STRUCT *pcu, int *ds);
void colmout(char *nfice, char *nfics, AMINO_STRUCT *paa, FILE *summary, MENU_STRUCT *pm);
void rowout(char *nfice, char *nfics, char *ncout, char sp, FILE *summary, MENU_STRUCT *pm, COA_STRUCT *pcoa);
void PrepAFC(char *nfic, COA_STRUCT *pcoa);
void inertialig(char *inertia_out, char *filen, FILE *summary, MENU_STRUCT *pm, COA_STRUCT *pcoa);
void inertiacol(char *inertia_out, FILE *summary, MENU_STRUCT *pm, COA_STRUCT *pcoa, AMINO_STRUCT *paa);
void suprow(int num_seq, char *nficvp, char *nfictasup,
            char *nficlisup, char *option, char sp, FILE *summary, COA_STRUCT *pcoa);
