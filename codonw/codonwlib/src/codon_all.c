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
/* initialize_point    assigns genetic code dependent parameters to structs*/
/* codon_usage_tot    Counts codon and amino acid usage                   */
/* ident_codon        Converts codon into a numerical value in range 1-64 */
/* how_synon          Calculates how synonymous each codon is             */
/* how_synon_aa       Calculates how synonymous each AA is                */
/* clean_up           Re-zeros various internal counters and arrays       */
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

/********************* Initilize Pointers**********************************/
/* Various pointers to structures are assigned here dependent on the      */
/* genetic code chosen.                                                   */
/* paa                points to a struct containing Amino Acid names      */
/* pap                points to amino acid properties                     */
/* pcai               points to Adaptation values used to calc CAI        */
/* pfop               points to a struct describing optimal codons        */
/* pcbi               points to the same structure as pfop                */
/* pcu                points to data which has the translation of codons  */
/* ds                 is a struct describing how synonymous a codon is    */
/* da                 is a struct describing the size of each AA family   */
/*                    included/excluded from any COA analysis             */
/**************************************************************************/
int initialize_point(char code, char fop_species, char cai_species, MENU_STRUCT *pm, REF_STRUCT *ref)
{
   pm->paa = ref->amino_acids;
   pm->pap = ref->amino_prop;
   pm->pcai = &(ref->cai[cai_species]);
   pm->pfop = &(ref->fop[fop_species]);
   pm->pcbi = &(ref->fop[fop_species]);
   pm->pcu = &(ref->cu[code]);

   static int dds[65];
   how_synon(dds, pm->pcu);
   pm->ds = dds;

   static int dda[22];
   how_synon_aa(dda, pm->pcu);
   pm->da = dda;

   fprintf(pm->my_err, "Genetic code set to %s %s\n", pm->pcu->des, pm->pcu->typ);

   return 0;
}

/*******************How Synonymous is this codon  *************************/
/* Calculate how synonymous a codon is by comparing with all other codons */
/* to see if they encode the same AA                                      */
/**************************************************************************/
int how_synon(int dds[], GENETIC_CODE_STRUCT *pcu)
{
   int x, i;

   for (x = 0; x < 65; x++)
      dds[x] = 0;

   for (x = 1; x < 65; x++)
      for (i = 1; i < 65; i++)
         if (pcu->ca[x] == pcu->ca[i])
            dds[x]++;

   return 0;
}
/*******************How Synonymous is this AA     *************************/
/* Calculate how synonymous an amino acid is by checking all codons if    */
/* they encode this same AA                                               */
/**************************************************************************/
int how_synon_aa(int dda[], GENETIC_CODE_STRUCT *pcu)
{
   int x;

   for (x = 0; x < 22; x++)
      dda[x] = 0;

   for (x = 1; x < 65; x++)
      dda[pcu->ca[x]]++;

   return 0;
}

/****************** Codon Usage Counting      *****************************/
/* Counts the frequency of usage of each codon and amino acid this data   */
/* is used throughout CodonW                                              */
/* pcu->ca contains codon to amino acid translations for the current code */
/* and is assigned in initialise point                                    */
/**************************************************************************/
int codon_usage_tot(char *seq, long *codon_tot, int *valid_stops, long ncod[], long naa[], GENETIC_CODE_STRUCT *pcu)
{
   char codon[4];
   int icode;
   unsigned int i;
   unsigned int seqlen = (int)strlen(seq);

   for (i = 0; i < seqlen - 2; i += 3)
   {
      strncpy(codon, (seq + i), 3);
      icode = ident_codon(codon);
      ncod[icode]++;             /*increment the codon count */
      naa[pcu->ca[icode]]++; /*increment the AA count    */
      (*codon_tot)++;
   }

   if (seqlen % 3)
   {             /*if last codon was partial */
      icode = 0; /*set icode to zero and     */
      ncod[0]++; /*increment untranslated    */
   }             /*codons                    */

   if (pcu->ca[icode] == 11)
      (*valid_stops)++;

   return icode;
}

/****************** Ident codon               *****************************/
/* Converts each codon into a numerical array (codon) and converts this   */
/* array into a numerical value in the range 0-64, zero is reserved for   */
/* codons that contain at least one unrecognised base                     */
/**************************************************************************/
int ident_codon(char *codon)
{
   int icode = 0;
   int x;

   for (x = 0; x < 3; x++)
   {
      switch (codon[x])
      {
      case 'T':
      case 't':
      case 'U':
      case 'u':
         codon[x] = (char)1;
         continue;
      case 'C':
      case 'c':
         codon[x] = (char)2;
         continue;
      case 'A':
      case 'a':
         codon[x] = (char)3;
         continue;
      case 'G':
      case 'g':
         codon[x] = (char)4;
         continue;
      case '\0':
         return 0;
      default:
         codon[x] = (char)0;
         break;
      }
   }
   if (codon[0] * codon[1] * codon[2] != 0)
      icode = (codon[0] - 1) * 16 + codon[1] + (codon[2] - 1) * 4;
   else
      icode = 0;

   return icode;
}

int count_codons(long* ncod, long *loc_cod_tot) {
   int i;

   *loc_cod_tot = 0;
   for (i = 1; i < 65; i++)
      (*loc_cod_tot) += ncod[i];

   return 0;
}


/******************  Clean up               *******************************/
/* Called after each sequence has been completely read from disk          */
/* It re-zeros all the main counters, but is not called when concatenating*/
/* sequences together                                                     */
/**************************************************************************/
int clean_up(long *nncod, long *nnaa, int *valid_stops)
{
   int x;
   int i;

   for (x = 0; x < 65; x++)
      nncod[x] = 0;
   for (x = 0; x < 23; x++)
      nnaa[x] = 0;
   *valid_stops = 0;

   return 0;
}
