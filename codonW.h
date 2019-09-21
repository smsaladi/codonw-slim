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

#define ARB_UNIT 100                 /* used to define the array*/
#define MAX_GENE (ARB_UNIT * 3)      /* seq, which holds readin */
#define LINE_LENGTH (ARB_UNIT + 100) /* sequence data           */
#define GARG_EXACT 0x800             /* used in function gargs  */
#define GARG_NEXT 0x1000             /* used in function gargs  */
#define GARG_THERE 0x2000            /* used in function gargs  */
#define GARG_SUBSQ 0x4000            /* used in function gargs  */
#define MAX_ARGS 100                 /* used in function gargs  */
#define MAX_FILENAME_LEN 90 /* max filename             */

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

  char codonW;   /* am I codonW              */
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

  char seperator; /* column separator         */
  char coa;       /* calculate a COA or not ? */

  char code;   /* which genetic code       */
  char f_type; /* which predefined fop val */
  char c_type; /* which predefined CAI val */

  char seq_type;                           /* DNA or Protein or CU     */
  char curr_infilename[MAX_FILENAME_LEN];  /* input filename           */
  char curr_outfilename[MAX_FILENAME_LEN]; /* .out filename            */
  char curr_tidyoutname[MAX_FILENAME_LEN]; /* .blk filename            */
  char fop_filen[MAX_FILENAME_LEN];        /* user fop filename        */
  char cai_filen[MAX_FILENAME_LEN];        /* user CAI filename        */
  char cbi_filen[MAX_FILENAME_LEN];        /* user CBI filename        */
  char curr_logfilename[MAX_FILENAME_LEN]; /* used for logging errors  */

  char junk[BUFSIZ + 1]; /* used to store char info  */
  char messages[300];    /* used to constuct messgs  */

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
} MENU_STRUCT;

extern AMINO_STRUCT *paa;
extern GENETIC_CODE_STRUCT *pcu;
extern FOP_STRUCT *pfop;
extern FOP_STRUCT *pcbi;
extern CAI_STRUCT *pcai;
extern MENU_STRUCT *pm;
extern COA_STRUCT *pcoa;
extern AMINO_PROP_STRUCT *pap;

extern CAI_STRUCT cai[];
extern GENETIC_CODE_STRUCT cu[];
extern FOP_STRUCT fop[];

extern COA_STRUCT coa;
extern AMINO_STRUCT amino_acids;
extern AMINO_PROP_STRUCT amino_prop;
extern MENU_STRUCT Z_menu;

extern char title[100];
extern char long_seq;
extern char last_base;

extern long int ncod[65];
extern long int naa[23];
extern long int din[3][16];
extern long int codon_tot;
extern long int master_ic;
extern long int fl_pos_start;
extern long int fl_pos_curr;
extern long int GC_TOT;
extern long int AT_TOT;
extern long int AA_TOT;
extern long int IUBC_TOT;
extern long int GAP_TOT;
extern long int num_sequence;
extern long int num_seq_int_stop;
extern long int non_std_char;
extern long int tot;
extern int last_aa;
extern int reg;
extern int valid_stops;
extern int valid_start;
extern int fram;
extern int *da;
extern int *ds;
extern int NumGeneticCodes;
extern int NumFopSpecies;
extern int NumCaiSpecies;

/****************** Function type declarations *****************************/

FILE *open_file(char *filename, char *mode);

int *how_synon(void);
int *how_synon_aa(void);
int *how_synon(void);
int *how_synon_aa(void);

int codon_usage_tot(char *seq, long int how_many);
int ident_codon(char *codon);
int codon_usage_out(FILE *fblkout, long int *ncod, int last_aa,
                    int valid_stops, char *info);
int rscu_usage_out(FILE *fblkout, long int *ncod, long int *naa);
int raau_usage_out(FILE *fblkout, long int *naa);
int aa_usage_out(FILE *fblkout, long int *naa);
int cai_out(FILE *foutput, long int *ncod);
int cbi_out(FILE *foutput, long int *ncod, long int *naa);
int fop_out(FILE *foutput, long int *ncod);
int hydro_out(FILE *foutput, long int *naa);
int aromo_out(FILE *foutput, long int *naa);
int toutput(FILE *fblkout, char *seq);
int output_long(FILE *fblkout, char *seq);
int cutab_out(FILE *fblkout, long *ncod, long *naa);
int dinuc_out(FILE *fblkout, char *title);
int fileclose(FILE **file_pointer);
int clean_up(long int *ncod, long int *naa);
int initilize_point(char code, char fop_type, char cai_type);
int initilize_coa(char code);
int proc_comm_line(int *argc, char ***arg_list);
int my_exit(int exit_value, char *message);

int dinuc_count(char *seq, long int tot);
int tidy(FILE *finput, FILE *foutput, FILE *fblkout,
         FILE *fcoaout);

long int codon_error(int last_aa, int valid_stops, char *title,
                     char error_level);

float enc_out(FILE *foutput, long int *ncod, long int *naa);
double inertot(void);

char *get_aa(int one_or_3_letter, char *the_dna_word);
char *garg(int argc, char *argv[], const char *targ, int mode);
char coa_raw_out(FILE *fcoaout, long *ncod, long *naa, char *title);

void sorted_by_axis1(double *ax1, int *sortax1, int lig);
void highlow(long int *low, long int *high, FILE *summ);

void asummary(void);
void vecalloc(double **vec, int n);
void vecalloc(double **vec, int n);
void writevec(double *v1, FILE *fic);
void lecmat(double **tab, char *nfic);
void freetab(double **tab);
void freevec(double *vec);
void taballoc(double ***tab, int l1, int c1);
void lecvec(double *v1, char *nfic);
void ecrmat(double **tab, char *nfic);
void ecrvec(double *v1, char *nfic);
void scalmat(double **tab, double r);
void scalvec(double *v1, double r);
void sqrvec(double *v1);
void prodmatAAtB(double **a, double **b);
void prodmatABC(double **a, double **b, double **c);
void prodmatAtAB(double **a, double **b);
void ecrmatred(double **tab, int c1, char *nfic);
void readvec(double *v1, FILE *fic);
void lecvalpro(double *v1, char *nfic);
void editvalpro(FILE *ficlist, double *vp, int n, double s);
void DiagoRC(FILE *summary);
void gc_out(FILE *foutput, FILE *fblkout, int which);
void base_sil_us_out(FILE *foutput, long int *ncod, long int *naa);
void bintext(char *nfice, char *nfics);
void select_coa(char choice);
void textbin(char *filein, char *fileout);
void colmout(char *nfice, char *nfics, AMINO_STRUCT *paa, FILE *summary);
void output(char *seq, FILE *foutput, FILE *fblkout, FILE *fcoaout);
void rowout(char *nfice, char *nfics, char *ncout, FILE *summary);
void PrepAFC(char *nfic);
void inertialig(char *inertia_out, char *filen, FILE *summary);
void inertiacol(char *inertia_out, FILE *summary);
void selectcol(char *nfic, double *col, int numcol);
void gen_cusort_fop(int *sortax1, int lig, FILE *fnam, FILE *summ);
void dot(int y, long int period);
void DiagoComp(int n0, double **w, double *d, int *rang);
void suprow(int num_seq, char *nficvp, char *nfictasup,
            char *nficlisup, char *option, FILE *summary);
