/*************************************************************************

CodonW codon usage analysis package

    Copyright (C) 2005            John F. Peden
    Copyright (C) 2020            Shyam Saladi

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; version 2 of the License.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
675 Mass Ave, Cambridge, MA 02139, USA.

*************************************************************************

This file contains functions used to calculate single-value indicies
related to a given gene sequence. Functions *_out are not used in the
Python bindings and may be removed in the future.

************************************************************************/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <stdbool.h>

#include "codonW.h"

/****************** Silent Base Usage     *******************************/
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

/***************** Codon Adaptation Index   *************************/
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

/*****************     Codon Bias Index     **************************/
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

/****************** Frequency of OPtimal codons  ********************/
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

/***************  Effective Number of Codons   *********************/
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
         /* the calculation                   */
         *enc_tot += (float)fold[z] / (float)averb;
         if (*enc_tot > 61)
            *enc_tot = 61;
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
   else
      fprintf(foutput, "%5.2f%c", enc_tot, sp);
      
   return 0;
}


/*********************  Hydropathy        **********************************/
int hydro(long *nnaa, float *hydro, float hydro_ref[22])
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
         *hydro += ((float)nnaa[i] / (float)a2_tot) * hydro_ref[i];

   return 0;
}

int hydro_out(FILE *foutput, long *nnaa, char* title, MENU_STRUCT *pm)
{
   float out;
   char sp = pm->separator;

   int retval = hydro(nnaa, &out, pm->pap->hydro);
   if(retval == 1) {
      fprintf(pm->my_err, "Warning %.20s appear to be too short\n", title);
      fprintf(pm->my_err, "No output was written to file\n");
      return 1;
   }
      
   fprintf(foutput, "%8.6f%c", out, sp);
   return 0;

}

/**************** Aromaticity ***********************************************/
int aromo(long *nnaa, float *aromo, int aromo_ref[22])
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
         *aromo += ((float)nnaa[i] / (float)a1_tot) * (float)aromo_ref[i];
   
   return 0;
}

int aromo_out(FILE *foutput, long *nnaa, char* title, MENU_STRUCT *pm)
{
   float out;
   char sp = pm->separator;

   int retval = aromo(nnaa, &out, pm->pap->aromo);
   if(retval == 1) {
      fprintf(pm->my_err, "Warning %.20s appear to be too short\n", title);
      fprintf(pm->my_err, "No output was written to file\n");
      return 1;
   }
      
   fprintf(foutput, "%8.6f%c", out, sp);
   return 0;
}
