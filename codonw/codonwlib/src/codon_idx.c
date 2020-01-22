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
/* -----------------------       codon_idx.C     ------------------------ */
/* This file contains most of the codon usage analysis subroutines        */
/* except for the COA analysis                                            */
/* Internal subroutines and functions                                     */
/* base_sil_us_out    Write out base composition at silent sites          */
/* cai_out            Write out CAI usage                                 */
/* cbi_out            Write out codon bias index                          */
/* fop_out            Write out Frequency of Optimal codons               */
/* enc_out            Write out Effective Number of codons                */
/* gc_out             Writes various analyses of base usage               */
/* hydro_out          Write out Protein hydropathicity                    */
/* aromo_out          Write out Protein aromaticity                       */
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

#include "codonW.h"

/******************  Base Silent output     *******************************/
/* Calculates and write the base composition at silent sites              */
/* normalised as a function of the possible usage at that silent site with*/
/* changing the amino acid composition of the protein. It is inspired by  */
/* GC3s but is much more complicated to calculate as not every AA has the */
/* option to use any base at the third position                           */
/* All synonymous AA can select between a G or C though                   */
/**************************************************************************/
int base_sil_us(long *nncod, long *nnaa, double base_sil[], int *ds, int *da, GENETIC_CODE_STRUCT *pcu)
{
   int id, i, x, y, z;
   long bases_s[4]; /* synonymous GCAT bases               */
   long cb[4]; /* codons that could have been GCAT    */

   int done[4];

   for (x = 0; x < 4; x++)
   {
      cb[x] = 0;
      bases_s[x] = 0;
   } /* blank the arrays                    */

   for (x = 1; x < 5; x++)
      for (y = 1; y < 5; y++)
         for (z = 1; z < 5; z++)
         { /* look at all 64 codons               */
            id = (x - 1) * 16 + y + (z - 1) * 4;

            if (*(ds + id) == 1 || pcu->ca[id] == 11)
               continue;                 /* if no synon skip to next       codon */
            bases_s[z - 1] += nncod[id]; /* count No. codon ending in base X     */
         }

   for (i = 1; i < 22; i++)
   {
      for (x = 0; x < 4; x++) /* don't want to count bases in 6 fold  */
         done[x] = false;     /* sites twice do we so we remember     */

      if (i == 11 || *(da + i) == 1)
         continue; /* if stop codon skip, or AA not synony */

      for (x = 1; x < 5; x++) /* else add aa to could have ended count */
         for (y = 1; y < 5; y++)
            for (z = 1; z < 5; z++)
            {
               id = (x - 1) * 16 + y + (z - 1) * 4;
               /* assign codon values in range 1-64                           */
               if (pcu->ca[id] == i && done[z - 1] == false)
               {
                  /* encode AA i which we know to be synon so add could_be_x ending*/
                  /* by the Number of that amino acid                              */
                  cb[z - 1] += nnaa[i];
                  done[z - 1] = true; /* don't look for any more or we might   */
                                      /* process leu+arg+ser twice             */
               }
            }
   }

   /* Now the easy bit ... just output the results                */
   for (i = 0; i < 4; i++)
   {
      if (cb[i] > 0)
         base_sil[i] = (double)bases_s[i] / (double)cb[i];
      else
         base_sil[i] = 0;
   }

   return 0;
}

int base_sil_us_out(FILE *foutput, long *nncod, long *nnaa, MENU_STRUCT *pm)
{
   double base_sil[4];

   base_sil_us(nncod, nnaa, base_sil, pm->ds, pm->da, pm->pcu);

   char sp = pm->separator;

   for (int i = 0; i < 4; i++)
      fprintf(foutput, "%6.4f%c", base_sil[i], sp);

   return 0;
}

/*****************Codon Adaptation Index output   *************************/
/* Codon Adaptation Index (CAI) (Sharp and Li 1987). CAI is a measurement */
/* of the relative adaptiveness of the codon usage of a gene towards the  */
/* codon usage of highly expressed genes. The relative adaptiveness (w) of*/
/* each codon is the ratio of the usage of each codon, to that of the most*/
/* abundant codon for the same amino acid. The relative adaptiveness of   */
/* codons for albeit a limited choice of species, can be selected from the*/
/* Menu. The user can also input a personal choice of values. The CAI     */
/* index is defined as the geometric mean of these relative adaptiveness  */
/* values. Non-synonymous codons and termination codons (genetic code     */
/* dependent) are excluded. To aid computation, the CAI is calculated as  */
/* using a natural log summation, To prevent a codon having a relative    */
/* adaptiveness value of zero, which could result in a CAI of zero;       */
/* these codons have fitness of zero (<.0001) are adjusted to 0.01        */
/**************************************************************************/
int cai(long *nncod, double *sigma, int *ds, CAI_STRUCT *pcai, GENETIC_CODE_STRUCT *pcu)
{
   long totaa = 0;
   int x;
   
   for (x = 1, *sigma = 0; x < 65; x++)
   {
      if (pcu->ca[x] == 11 || *(ds + x) == 1)
         continue;
      if (pcai->cai_val[x] < 0.0001)     /* if value is effectively zero       */
         pcai->cai_val[x] = 0.01F;       /* make it .01 */
      *sigma += (double)*(nncod + x) * log((double)pcai->cai_val[x]);
      totaa += *(nncod + x);
   }

   if (totaa)
   { /* catch floating point overflow error*/
      *sigma = *sigma / (double)totaa;
      *sigma = exp(*sigma);
   }
   else
      *sigma = 0;

   return 0;
}

int cai_out(FILE *foutput, long *nncod, MENU_STRUCT *pm)
{
   double sigma;

   fprintf(stderr, "Using %s (%s) w values to calculate CAI\n",
           pm->pcai->des, pm->pcai->ref);

   cai(nncod, &sigma, pm->ds, pm->pcai, pm->pcu);

   char sp = pm->separator;
   fprintf(foutput, "%5.3f%c", sigma, sp);

   return 0;
}

/*****************Codon Bias Index output        **************************/
/* Codon bias index is a measure of directional codon bias, it measures   */
/* the extent to which a gene uses a subset of optimal codons.            */
/* CBI = ( Nopt-Nran)/(Nopt-Nran) Where Nopt = number of optimal codons;  */
/* Ntot = number of synonymous codons; Nran = expected number of optimal  */
/* codons if codons were assigned randomly. CBI is similar to Fop as used */
/* by Ikemura, with Nran used as a scaling factor. In a gene with extreme */
/* codon bias, CBI will equal 1.0, in a gene with random codon usage CBI  */
/* will equal 0.0. Note that it is possible for Nopt to be less than Nran.*/
/* This results in a negative value for CBI.                              */
/* ( Bennetzen and Hall 1982 )                                            */
/**************************************************************************/
int cbi(long *nncod, long *nnaa, float *fcbi, int *ds, int *da, GENETIC_CODE_STRUCT *pcu, FOP_STRUCT *pcbi)
{
   long tot_cod = 0;
   long opt = 0;
   float exp_cod = 0.0F;
   int x;
   static char first_call_cbi = true;
   static char has_opt_info[22];

   if (first_call_cbi)
   { /* have we been called already   */
      /* initilise has_opt_info             */
      for (x = 1; x < 22; x++)
         has_opt_info[x] = 0;

      for (x = 1; x < 65; x++)
      {
         if (pcu->ca[x] == 11 || *(ds + x) == 1)
            continue;
         if (pcbi->fop_cod[x] == 3)
            has_opt_info[pcu->ca[x]]++;
      }

      first_call_cbi = false; /*      this won't be called again      */
   }

   for (x = 1; x < 65; x++)
   {
      if (!has_opt_info[pcu->ca[x]])
         continue;
      switch ((int)pcbi->fop_cod[x])
      {
      case 3:
         opt += nncod[x];
         tot_cod += nncod[x];
         exp_cod += (float)nnaa[pcu->ca[x]] / (float)da[pcu->ca[x]];
         break;
      case 2:
      case 1:
         tot_cod += *(nncod + x);
         break;
      default:
         fprintf(stderr, " Serious error in CBI information found"
                         " an illegal CBI value of %c for codon %i"
                         " permissible values are \n 1 for non-optimal"
                         " codons\n 2 for common codons\n"
                         " 3 for optimal codons\n",
                 pcbi->fop_cod[x], x);
         return 1;
      } /*                   end of switch     */
   }    /*                   for (    )        */

   if (tot_cod - exp_cod)
      *fcbi = (opt - exp_cod) / (tot_cod - exp_cod);
   else
      *fcbi = 0.0F;

   return 0;
}

int cbi_out(FILE *foutput, long *nncod, long *nnaa, MENU_STRUCT *pm)
{
   float fcbi;

   fprintf(stderr, "Using %s (%s) \noptimal codons to calculate CBI\n",
           pm->pcbi->des, pm->pcbi->ref);

   cbi(nncod, nnaa, &fcbi, pm->ds, pm->da, pm->pcu, pm->pcbi);

   char sp = pm->separator;
   fprintf(foutput, "%5.3f%c", fcbi, sp); /* CBI     QED     */

   return 0;
}

/****************** Frequency of OPtimal codons output  ********************/
/* Frequency of Optimal codons (Fop) (Ikemura 1981). This index, is ratio  */
/* of optimal codons to synonymous codons (genetic code dependent). Optimal*/
/* codons for several species are in-built and can be selected using Menu 3*/
/* By default, the optimal codons of E. coli are assumed. The user may also*/
/* enter a personal choice of optimal codons. If rare synonymous codons    */
/* have been identified, there is a choice of calculating the original Fop */
/* index or a modified index. Fop values for the original index are always */
/* between 0 (where no optimal codons are used) and 1 (where only optimal  */
/* codons are used). When calculating the modified Fop index, any negative */
/* values are adjusted to zero.                                            */
/***************************************************************************/
// factor_in_rare
// If non-optimal codons are identified in the set of optimal codons selected, use
// in the calculation of a modified Fop, (Fop=(opt-rare)/total). Otherwise, the original
// formulae (Fop=opt/total) is used.
int fop(long *nncod, float *ffop, int *ds, bool factor_in_rare, GENETIC_CODE_STRUCT *pcu, FOP_STRUCT *pfop)
{
   long nonopt = 0;
   long std = 0;
   long opt = 0;
   int x;
   
   /* initilise has_opt_info             */
   bool has_opt_info[22];
   for (x = 1; x < 22; x++)
      has_opt_info[x] = false;
   for (x = 1; x < 65; x++)
   {
      if (pcu->ca[x] == 11 || *(ds + x) == 1)
         continue;

      if (pfop->fop_cod[x] == 3)
         has_opt_info[pcu->ca[x]] = true;

      if (pfop->fop_cod[x] == 1 && factor_in_rare == true)
         has_opt_info[pcu->ca[x]] = true;
   }

   for (x = 1; x < 65; x++)
   {
      if (!has_opt_info[pcu->ca[x]])
         continue;

      switch ((int)pfop->fop_cod[x])
      {
      case 3:
         opt += *(nncod + x);
         break;
      case 2:
         std += *(nncod + x);
         break;
      case 1:
         nonopt += *(nncod + x);
         break;
      default:
         fprintf(stderr, " Serious error in fop information found"
                         " an illegal fop value of %c for codon %i"
                         " permissible values are \n 1 for non-optimal"
                         " codons\n 2 for common codons\n"
                         " 3 for optimal codons\n",
                 pfop->fop_cod[x], x);
         fprintf(stderr, "opt %ld, std %ld, nonopt %ld\n", opt, std, nonopt);
         return 1;
      }
   }
   /* only ask this once  ...            */

   if (factor_in_rare && (opt + nonopt + std))
      *ffop = (float)(opt - nonopt) / (float)(opt + nonopt + std);
   else if ((opt + nonopt + std))
      *ffop = (float)opt / (float)(opt + nonopt + std);
   else
      *ffop = 0.0;

   return 0;
}

int fop_out(FILE *foutput, long *nncod, MENU_STRUCT *pm) {
   float ffop;

   fprintf(stderr, "Using %s (%s)\noptimal codons to calculate Fop\n",
            pm->pfop->des, pm->pfop->ref);

   bool factor_in_rare = false;
   int retval = fop(nncod, &ffop, pm->ds, factor_in_rare, pm->pcu, pm->pfop);

   char sp = pm->separator;
   fprintf(foutput, "%5.3f%c", ffop, sp);

   return 0;
}

/***************  Effective Number of Codons output   *********************/
/* The effective number of codons (NC) (Wright 1990). This index is a     */
/* simple measure of overall codon bias and is analogous to the effective */
/* number of alleles measure used in population genetics. Knowledge of the*/
/* optimal codons or a reference set of highly expressed genes is not     */
/* needed when calculating this index. Initially the homozygosity for each*/
/* amino acid is estimated from the squared codon frequencies.            */
/**************************************************************************/
int enc(long *nncod, long *nnaa, float *enc_tot, int *da, GENETIC_CODE_STRUCT *pcu)
{
   int numaa[9];
   int fold[9];
   int error_t = false;
   int i, z, x;
   double totb[9];
   double averb = 0, bb = 0, k2 = 0, s2 = 0;
   *enc_tot = 0.0F;

   /* don't assume that 6 is the largest possible amino acid family assume 9*/
   for (i = 0; i < 9; i++)
   {
      fold[i] = 0; /* initialise arrays to zero             */
      totb[i] = 0.0;
      numaa[i] = 0;
   }

   for (i = 1; i < 22; i++)
   { /* for each amino acid                  */
      if (i == 11)
         continue; /* but not for stop codons              */

      if (*(nnaa + i) <= 1) /* if this aa occurs once then skip     */
         bb = 0;
      else
      {
         for (x = 1, s2 = 0; x < 65; x++)
         {
            /* Try all codons but we are only looking for those that encode*/
            /* amino amid i, saves having to hard wire in any assumptions  */
            if (pcu->ca[x] != i)
               continue; /* skip is not i       */

            if (*(nncod + x) == 0) /* if codons not used then              */
               k2 = 0.0;           /* k2 = 0                               */
            else
               k2 = pow(((double)*(nncod + x) / (double)*(nnaa + i)),
                        (double)2);

            s2 += k2; /* sum of all k2's for aa i             */
         }
         bb = (((double)*(nnaa + i) * s2) - 1.0) / /* homozygosity        */
              (double)(*(nnaa + i) - 1.0);
      }

      if (bb > 0.0000001)
      {
         totb[*(da + i)] += bb; /* sum of all bb's for amino acids  */
                                /* which have z alternative codons  */
         numaa[*(da + i)]++;    /* where z = *(da+i)                */
      }
      /* numaa is no of aa that were z    */
      fold[*(da + i)]++; /* fold z=4 can have 9 in univ code */
   }                     /* but some aa may be absent from   */
                         /* gene therefore numaa[z] may be 0 */
   *enc_tot = (float)fold[1];

   for (z = 2, averb = 0, error_t = false; z <= 8; z++)
   {
      /* look at all values of z if there  */
      if (fold[z])
      { /* are amino acids that are z fold   */
         if (numaa[z] && totb[z] > 0)
            averb = totb[z] / numaa[z];
         else if (z == 3 && numaa[2] && numaa[4] && fold[z] == 1)
            /* special case                      */
            averb = (totb[2] / numaa[2] + totb[4] / numaa[4]) * 0.5;
         else
         {
            fprintf(stderr, "%i amino acids with %i synonymous codons\n", numaa[z], z);
            fprintf(stderr, "\t -- Nc was not calculated\n");
            return 1;
         }
         *enc_tot += (float)fold[z] / (float)averb;
         /* the calculation                   */
      }
   }

   return 0;
}

int enc_out(FILE *foutput, long *nncod, long *nnaa, MENU_STRUCT *pm)
{
   char sp = pm->separator;
   float enc_tot;

   int retval = enc(nncod, nnaa, &enc_tot, pm->da, pm->pcu);

   if (retval == 1)
      fprintf(foutput, "*****%c", sp);
   else if (enc_tot <= 61)
      fprintf(foutput, "%5.2f%c", enc_tot, sp);
   else
      fprintf(foutput, "61.00%c", sp);

   return 0;
}


/*********************  hydro_out        **********************************/
/* The general average hydropathicity or (GRAVY) score, for the hypothet- */
/* ical translated gene product. It is calculated as the arithmetic mean  */
/* of the sum of the hydropathic indices of each amino acid. This index   */
/* was used to quantify the major COA trends in the amino acid usage of   */
/* E. coli genes (Lobry, 1994).                                           */
/* Calculates and outputs total protein hydropathicity based on the Kyte  */
/* and Dolittle Index of hydropathicity (1982)                            */
/* nnaa               Array with frequency of amino acids                 */
/* paa                points to a struct containing Amino Acid values     */
/* pap->hydro         Pointer to hydropathicity values for each AA        */
/**************************************************************************/
int hydro(long *nnaa, float *hydro, AMINO_PROP_STRUCT *pap)
{
   long a2_tot = 0;
   *hydro = 0.0F;
   int i;

   for (i = 1; i < 22; i++)
      if (i != 11)
         a2_tot += nnaa[i];

   /* whow   .. no amino acids what happened     */
   if (!a2_tot)
      return 1;

   for (i = 1; i < 22; i++)
      if (i != 11)
         *hydro += ((float)nnaa[i] / (float)a2_tot) * (float)pap->hydro[i];

   return 0;
}

int hydro_out(FILE *foutput, long *nnaa, char* title, MENU_STRUCT *pm)
{
   float out;
   char sp = pm->separator;

   int retval = hydro(nnaa, &out, pm->pap);
   if(retval == 1) {
      fprintf(pm->my_err, "Warning %.20s appear to be too short\n", title);
      fprintf(pm->my_err, "No output was written to file\n");
      return 1;
   }
      
   fprintf(foutput, "%8.6f%c", out, sp);
   return 0;

}

/**************** Aromo_out ***********************************************/
/* Aromaticity score of protein. This is the frequency of aromatic amino  */
/* acids (Phe, Tyr, Trp) in the hypothetical translated gene product      */
/* nnaa               Array with frequency of amino acids                 */
/* paa                points to a struct containing Amino Acid values     */
/* pap->aromo         Pointer to aromaticity values for each AA           */
/**************************************************************************/
int aromo(long *nnaa, float *aromo, AMINO_PROP_STRUCT *pap)
{
   long a1_tot = 0;
   *aromo = 0.0F;
   int i;

   for (i = 1; i < 22; i++)
      if (i != 11)
         a1_tot += nnaa[i];

   if (!a1_tot)
      return 1;

   for (i = 1; i < 22; i++)
      if (i != 11)
         *aromo += ((float)nnaa[i] / (float)a1_tot) * (float)pap->aromo[i];
   
   return 0;
}

int aromo_out(FILE *foutput, long *nnaa, char* title, MENU_STRUCT *pm)
{
   float out;
   char sp = pm->separator;

   int retval = aromo(nnaa, &out, pm->pap);
   if(retval == 1) {
      fprintf(pm->my_err, "Warning %.20s appear to be too short\n", title);
      fprintf(pm->my_err, "No output was written to file\n");
      return 1;
   }
      
   fprintf(foutput, "%8.6f%c", out, sp);
   return 0;
}