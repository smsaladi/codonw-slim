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
/*                                                                        */
/* -----------------------        codon_us.C     ------------------------ */
/* This file contains most of the codon usage analysis subroutines        */
/* except for the COA analysis                                            */
/* Internal subroutines and functions                                     */
/* codon_usage_out    Write out Codon Usage to file                       */
/* rscu_usage_out     Write out RSCU                                      */
/* raau_usage_out     Write out normalised amino acid usage               */
/* aa_usage_out       Write out amino acid usage                          */
/* gc_out             Writes various analyses of base usage               */
/* cutab_out          Write a nice tabulation of the RSCU+CU+AA           */
/* dinuc_count        Count the dinucleotide usage                        */
/* dinuc_out          Write out dinucleotide usage                        */
/*                                                                        */
/*                                                                        */
/* External subroutines to codon_us.c                                     */
/* my_exit            Controls exit from CodonW closes any open files     */
/* tidy               reads the input data                                */
/* output             called from tidy to decide what to do with the data */
/* open_file          Open files, checks for existing files               */
/*                                                                        */
/**************************************************************************/

/*
Codon error checking 
Check for start, stop codons, internal stop, non-translatable and partial codons
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <stdbool.h>

#include "../include/codonW.h"

/****************** Codon Usage Out           *****************************/
/* Writes codon usage output to file. Note this subroutine is only called */
/* when machine readable output is selected, otherwise cutab_out is used  */
/**************************************************************************/
int codon_usage_out(FILE *fblkout, long *nncod, char *ttitle, MENU_STRUCT *pm)
{
   GENETIC_CODE_STRUCT *pcu = pm->pcu; 

   long ccodon_tot = 0;
   int x;
   char sp = pm->separator;

   count_codons(nncod, &ccodon_tot);

   /*example of output                                                     */
   /*0,0,0,0,3,2,2,0,0,0,0,0,0,3,0,0,                                      */
   /*0,0,0,4,3,4,1,7,0,0,0,0,3,1,3,1,Codons=100                              */
   /*0,0,0,0,10,6,3,0,0,0,0,0,1,1,12,0,Universal Genetic code              */
   /*0,0,0,3,7,5,7,9,0,1,1,1,8,4,5,0,MLSPCOPER.PE1                         */

   for (x = 1; x < 65; x++)
   {

      fprintf(fblkout, "%ld%c", nncod[x], sp);

      switch (x)
      {
      case 16:
         fprintf(fblkout, "\n");
         break;
      case 32:
         fprintf(fblkout, "Codons=%ld\n", ccodon_tot);
         break;
      case 48:
         fprintf(fblkout, "%.30s\n", pcu->des);
         break;
      case 64:
         fprintf(fblkout, "%.20s\n", ttitle);
         break;
      default:
         break;
      }
   }
   return 0;
}
/******************  RSCU  Usage out          *****************************/
/* Writes Relative synonymous codon usage output to file. Note this subrou*/
/* tine is only called if machine readable output is selected             */
/* If human readable format was selected then what the user really wanted */
/* was cutab so this is automatically selected in codons.c                */
/* RSCU values are genetic codon dependent                                */
/**************************************************************************/
int rscu_usage(long *nncod, long *nnaa, float rscu[], int *ds, GENETIC_CODE_STRUCT *pcu)
{
   int x;

   /* ds points to an array[64] of synonym values i.e. how synon its AA is  */
   for (x = 1; x < 65; x++) {
      if (nnaa[pcu->ca[x]] != 0)
         rscu[x] = ((float)nncod[x] / (float)nnaa[pcu->ca[x]]) * ((float)*(ds + x));
      else
         rscu[x] = 0.0;
   }

   
   return 0;
}

int rscu_usage_out(FILE *fblkout, long *nncod, long *nnaa, char* title, MENU_STRUCT *pm)
{
   float rscu[65];
   rscu_usage(nncod, nnaa, rscu, pm->ds, pm->pcu);

   int x;
   char sp = pm->separator;
   /* ds points to an array[64] of synonym values i.e. how synon its AA is  */

   for (x = 1; x < 65; x++)
   {
      fprintf(fblkout, "%5.3f%c", rscu[x], sp);

      if (x == 64)
         fprintf(fblkout, "%-20.20s", title);

      if (!(x % 16))
         fprintf(fblkout, "\n");
   }

   return 0;
}

/******************   RAAU output             *****************************/
/* Writes Relative amino acid usage output to file. Amino Acid usage is   */
/* normalised for gene length                                             */
/**************************************************************************/
int raau_usage(long nnaa[], double raau[])
{
   int i;

   // total No. of AAs
   long aa_tot = 0;

   for (i = 1; i < 22; i++)
      if (i != 11)
         aa_tot += nnaa[i];
   
   for (i = 0; i < 22; i++)
      if (i == 11)
         raau[i] = 0;
      else if (aa_tot)
         raau[i] = (double)nnaa[i] / (double)aa_tot;
      else /* No AminoAcids! */
         raau[i] = 0;

   return 0;
}

int raau_usage_out(FILE *fblkout, long *nnaa, char* title, MENU_STRUCT *pm)
{
   AMINO_STRUCT *paa = pm->paa;

   static char first_line = true;
   int i, x;
   char sp;

   sp = '\t';

   if (first_line)
   { /* if true write a header*/
      fprintf(fblkout, "%s", "Gene_name");

      for (i = 0; i < 22; i++)
         fprintf(fblkout, "%c%s", sp, paa->aa3[i]); /* three letter AA names*/
      fprintf(fblkout, "\n");
      first_line = false;
   }

   fprintf(fblkout, "%.30s", title);

   double raau[22];
   raau_usage(nnaa, raau);

   for (x = 0; x < 22; x++)
      fprintf(fblkout, "%c%.4f", sp, raau[x]);

   fprintf(fblkout, "\n");

   return 0;
}

/******************   AA usage output         *****************************/
/* Writes amino acid usage output to file.                                */
/**************************************************************************/
int aa_usage_out(FILE *fblkout, long *nnaa, char* title, MENU_STRUCT *pm)
{
   AMINO_STRUCT *paa = pm->paa;

   static char first_line = true;
   int i;
   char sp = pm->separator;

   if (first_line)
   {
      fprintf(fblkout, "%s", "Gene_name");

      for (i = 0; i < 22; i++)
         fprintf(fblkout, "%c%s", sp, paa->aa3[i]); /* 3 letter AA code     */

      fprintf(fblkout, "\n");
      first_line = false;
   }
   fprintf(fblkout, "%.20s", title);

   for (i = 0; i < 22; i++)
      fprintf(fblkout, "%c%li", sp, nnaa[i]);

   fprintf(fblkout, "\n");
   return 0;
}

/*******************   G+C output          *******************************/
/* This function is a real work horse, initially it counts base composit */
/* ion in all frames, length of gene, num synonymous codons, number of   */
/* non synonymous codons. Then dependent on the value for which used in  */
/* switch statement. We return various analyses of this data             */
/* if which ==1 then the output is very detailed, base by base etc.      */
/* if which ==2 then the output is for GC content only                   */
/* if which ==3 then the output is for GC3s (GC at synonymous 3rd posit) */
/* if which ==4 then the output is for L_sym                             */
/* if which ==5 then the output is for L_aa                              */
/* The output from this subroutine is in a tabular format if human read- */
/* able output is selected, and in columns if machine readable. Also the */
/* number of values reported changes as it is assumed the user has access*/
/* to a spreadsheet type programme if they are requesting tabular output */
/*************************************************************************/
int gc(int *ds, long *ncod, long bases[5], long base_tot[5], long base_1[5], long base_2[5], long base_3[5], long *tot_s, long *totalaa, double gc_metrics[], GENETIC_CODE_STRUCT *pcu)
{
   long id;
   // long bases[5]; /* base that are synonymous GCAT     */
   *tot_s = 0;
   *totalaa = 0;
   int x, y, z;

   for (x = 0; x < 5; x++)
   {
      bases[x] = 0; /* initialise array values to zero    */
      base_tot[x] = 0;
      base_1[x] = 0;
      base_2[x] = 0;
      base_3[x] = 0;
   }

   for (x = 1; x < 5; x++)
      for (y = 1; y < 5; y++)
         for (z = 1; z < 5; z++)
         { /* look at all 64 codons              */
            id = (x - 1) * 16 + y + (z - 1) * 4;

            if (pcu->ca[id] == 11)
               continue;             /* skip if a stop codon               */
            base_tot[x] += ncod[id]; /* we have a codon xyz therefore the  */
            base_1[x] += ncod[id];   /* frequency of each position for base*/
            base_tot[y] += ncod[id]; /* x,y,z are equal to the number of   */
            base_2[y] += ncod[id];   /* xyz codons .... easy               */
            base_tot[z] += ncod[id]; /* will be fooled a little if there   */
            base_3[z] += ncod[id];   /* non translatable codons, but these */
                                     /* are ignored when the avg is calc   */
            *totalaa += ncod[id];

            if (*(ds + id) == 1)
               continue; /* if not synon  skip codon           */

            bases[z] += ncod[id]; /* count no of codons ending in Z     */

            *tot_s += ncod[id]; /* count tot no of silent codons      */
         }


   /* Calculate metrics */
   typedef double lf;
   double metrics_local[] = {
      (lf)(base_tot[2] + base_tot[4]) / (lf)(*totalaa * 3),
      (lf)(bases[2] + bases[4]) / (lf)*tot_s,
      (lf)(base_tot[2] + base_tot[4] - bases[2] - bases[4]) / (lf)(*totalaa * 3 - *tot_s),
      (lf)(base_1[2] + base_1[4]) / (lf)(*totalaa),
      (lf)(base_2[2] + base_2[4]) / (lf)(*totalaa),
      (lf)(base_3[2] + base_3[4]) / (lf)(*totalaa),
      (lf)base_1[1] / (lf)*totalaa,
      (lf)base_2[1] / (lf)*totalaa,
      (lf)base_3[1] / (lf)*totalaa,
      (lf)base_1[2] / (lf)*totalaa,
      (lf)base_2[2] / (lf)*totalaa,
      (lf)base_3[2] / (lf)*totalaa,
      (lf)base_1[3] / (lf)*totalaa,
      (lf)base_2[3] / (lf)*totalaa,
      (lf)base_3[3] / (lf)*totalaa,
      (lf)base_1[4] / (lf)*totalaa,
      (lf)base_2[4] / (lf)*totalaa,
      (lf)base_3[4] / (lf)*totalaa
   };

   // Copy into output array
   for (x = 0; x < 18; x++)
      gc_metrics[x] = metrics_local[x];

   return 0;
}

int gc_out(FILE *foutput, FILE *fblkout, long *nncod, int which, char* title, MENU_STRUCT *pm)
{
   long bases[5]; /* base that are synonymous GCAT     */
   long base_tot[5];
   long base_1[5];
   long base_2[5];
   long base_3[5];
   long tot_s = 0;
   long totalaa = 0;
   double metrics[18];
   int i;

   gc(pm->ds, nncod, bases, base_tot, base_1, base_2, base_3, &tot_s, &totalaa, &metrics, pm->pcu);

   static char header = false;
   char sp = pm->separator;
   
   typedef double lf;

   if (!tot_s || !totalaa)
   {
      fprintf(pm->my_err, "Warning %.20s appear to be too short\n", title);
      fprintf(pm->my_err, "No output was written to file   \n");
      return 1;
   }

   switch ((int)which)
   {
   case 1: /* exhaustive output for analysis     */
      if (!header)
      { /* print a first line                 */
         fprintf(fblkout,
                  "Gene_description%cLen_aa%cLen_sym%cGC%cGC3s%cGCn3s%cGC1%cGC2"
                  "%cGC3%cT1%cT2%cT3%cC1%cC2%cC3%cA1%cA2%cA3%cG1%cG2%cG3\n",
                  sp, sp, sp, sp, sp, sp, sp, sp, sp, sp, sp, sp, sp, sp, sp, sp, sp, sp, sp, sp);
         header = true;
      }
      /* now print the information          */
      fprintf(fblkout, "%-.20s%c", title, sp);
      fprintf(fblkout, "%ld%c%ld", totalaa, sp, tot_s);
      for (i = 0; i < 18; i++)
         fprintf(fblkout, "%c%5.3f", sp, metrics[i]);
      fprintf(fblkout, "\n");
      break;
   case 2: /* a bit more simple ... GC content   */
      fprintf(foutput, "%5.3f%c", (lf)((base_tot[2] + base_tot[4]) / (lf)(totalaa * 3)), sp);
      break;
   case 3: /* GC3s                               */
      fprintf(foutput, "%5.3f%c", (lf)(bases[2] + bases[4]) / (lf)tot_s, sp);
      break;
   case 4: /* Number of synonymous codons        */
      fprintf(foutput, "%3li%c", tot_s, sp);
      break;
   case 5: /* Total length in translatable AA    */
      fprintf(foutput, "%3li%c", totalaa, sp);
      break;
   }

   return 0;
}

/**********************   cutab_out     ***********************************/
/* Generates a formatted table of codon, RSCU and amino acid usage        */
/* ds points to an array[64] of synonymous values                         */
/* it reveals how many synonyms there are for each aa                     */
/**************************************************************************/
int cutab_out(FILE *fblkout, long *nncod, long *nnaa, char* title, MENU_STRUCT *pm)
{
   AMINO_STRUCT *paa = pm->paa;
   GENETIC_CODE_STRUCT *pcu = pm->pcu;
   int *ds = pm->ds;

   int last_row[4];
   int x;
   char sp = pm->separator;

   for (x = 0; x < 4; x++)
      last_row[x] = 0;

   long codon_tot;
   count_codons(nncod, &codon_tot);

   for (x = 1; x < 65; x++)
   {
      if (last_row[x % 4] != pcu->ca[x])
         fprintf(fblkout, "%s%c%s%c", paa->aa3[pcu->ca[x]], sp, paa->cod[x], sp);
      else
         fprintf(fblkout, "%c%s%c", sp, paa->cod[x], sp);
      /* Sample of output *******************************************************/
      /*Phe UUU    0 0.00 Ser UCU    1 0.24 Tyr UAU    1 0.11 Cys UGU    1 0.67 */
      /*    UUC   22 2.00     UCC   10 2.40     UAC   17 1.89     UGC    2 1.33 */
      /*Leu UUA    0 0.00     UCA    1 0.24 TER UAA    0 0.00 TER UGA    1 3.00 */
      /*    UUG    1 0.12     UCG    6 1.44     UAG    0 0.00 Trp UGG    4 1.00 */
      /**************************************************************************/
      fprintf(fblkout, "%i%c%.2f%c",
               (int)nncod[x],
               sp, (nncod[x]) ? ((float)nncod[x] / (float)nnaa[pcu->ca[x]]) * (float)(*(ds + x)) : 0, sp);

      last_row[x % 4] = pcu->ca[x];

      if (!(x % 4))
         fprintf(fblkout, "\n");
      if (!(x % 16))
         fprintf(fblkout, "\n");
   }
   fprintf(fblkout, "%li codons in %16.16s (used %22.22s)\n\n",
           (long)codon_tot, title, pcu->des);
   return 0;
}

/********************  Dinuc_count  Dinuc_out  ****************************/
/* Count the frequency of all 16 dinucleotides in all three possible      */
/* reading frames.                                                        */
/*                                                                        */
/* Outputs the frequency of dinucleotides, either in fout rows per seq    */
/* if the output is meant to be in a human readable form, each row repre- */
/* senting a reading frame. The fourth row is the total of the all the    */
/* reading frames. Machine readable format writes all the data into a     */
/* single row                                                             */
/**************************************************************************/
int dinuc_count(char *seq, long din[3][16], long dinuc_tot[4], int *fram)
{
   int last, cur = 0;
   int i, x;

   long ttot = (long)strlen(seq);
   for (i = 0; i < ttot; i++)
   {
      last = cur;
      switch (seq[i])
      {
      case 't':
      case 'T':
      case 'u':
      case 'U':
         cur = 1;
         break;
      case 'c':
      case 'C':
         cur = 2;
         break;
      case 'a':
      case 'A':
         cur = 3;
         break;
      case 'g':
      case 'G':
         cur = 4;
         break;
      default:
         cur = 0;
         break;
      }
      if (cur == 0 || i == 0)
         continue; /* true if either of the base is not  */
                   /* a standard UTCG, or the current bas*/
                   /* is the start of the sequence       */
      din[*fram][((last - 1) * 4 + cur) - 1]++;
      if (++(*fram) == 3)
         *fram = 0; /* resets the frame to zero           */
   }
   
   for (x = 0; x < 4; x++)
      dinuc_tot[x] = 0;

   for (x = 0; x < 3; x++)
      for (i = 0; i < 16; i++)
      {
         dinuc_tot[x] += din[x][i]; /* count dinuc usage in each frame   */
         dinuc_tot[3] += din[x][i]; /* and total dinuc usage,            */
      }
   
   return 0;
}

int dinuc_out(char *seq, FILE *fblkout, char *ttitle, char sp) {
   static char called = false;
   char bases[5] = {'T', 'C', 'A', 'G'};
   int i, x, y;

   long din[3][16];
   long dinuc_tot[4];
   int fram = 0;

   for (x = 0; x < 3; x++)
      for (i = 0; i < 16; i++)
         din[x][i] = 0;

   for (x = 0; x < 4; x++)
      dinuc_tot[x] = 0;

   dinuc_count(seq, din, dinuc_tot, &fram);

   if (!called)
   { /* write out the first row as a header*/
      called = true;

      fprintf(fblkout, "%s", "title");
      for (y = 0; y < 4; y++)
      {
         fprintf(fblkout, "%c%s", sp, "frame");
         for (x = 0; x < 4; x++)
            for (i = 0; i < 4; i++)
               fprintf(fblkout, "%c%c%c", sp, bases[x], bases[i]);
      }

      fprintf(fblkout, "\n");
   } /* matches if (!called)               */

   /*Sample output   truncated  **********************************************/
   /*title         frame TT    TC    TA    TG    CT    CC    CA    CG    AT  */
   /*MLSPCOPER.PE1__ 1:2 0.024 0.041 0.016 0.008 0.049 0.041 0.033 0.098 ... */
   /*MLSPCOPER.PE1__ 2:3 0.000 0.195 0.000 0.098 0.000 0.138 0.008 0.073 ... */
   /*MLSPCOPER.PE1__ 3:1 0.008 0.016 0.000 0.033 0.033 0.107 0.172 0.262 ... */
   /*MLSPCOPER.PE1__ all 0.011 0.084 0.005 0.046 0.027 0.095 0.071 0.144 ... */
   /*MLSPCOPER.PE2__ 1:2 0.026 0.026 0.009 0.009 0.053 0.035 0.053 0.061 ... */
   /**************************************************************************/
   for (x = 0; x < 4; x++)
   {
      if (x == 0)
         fprintf(fblkout, "%-.15s%c", ttitle, sp);

      switch (x)
      {
      case 0:
         fprintf(fblkout, "%c1:2", sp);
         break;
      case 1:
         fprintf(fblkout, "%c2:3", sp);
         break;
      case 2:
         fprintf(fblkout, "%c3:1", sp);
         break;
      case 3:
         fprintf(fblkout, "%call", sp);
         break;
      }

      if (x == 3)
      {
         for (i = 0; i < 16; i++)
            if (dinuc_tot[x])
               fprintf(fblkout, "%c%5.3f", sp,
                       (float)(din[0][i] + din[1][i] + din[2][i]) /
                           (float)dinuc_tot[x]);
            else
               fprintf(fblkout, "%c%5.3f", sp, 0.00);
      }
      else
      {
         for (i = 0; i < 16; i++)
            if (dinuc_tot[x])
               fprintf(fblkout, "%c%5.3f", sp,
                       (float)din[x][i] / (float)dinuc_tot[x]);
            else
               fprintf(fblkout, "%c%5.3f", 0.00);
      }

      if (x == 3)
         fprintf(fblkout, "\n");
   }
   return 0;
}

