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
/* -----------------------        Codons.C       ------------------------ */
/* This file contains main() function and drives CodonW.                  */
/*                                                                        */
/* External subroutines and functions                                     */
/* proc_comm_line     process command line arguments                      */
/* initilize_point    assigns genetic code dependent parameters to structs*/
/* initilize_coa      selects the default codons to exclude from the      */
/*                    Correspondence Analysis                             */
/* clean_up           Re-zeros various internal counters and arrays       */
/* open_file          Open files, checks for existing files               */
/* fileclose          Closes files and returns a NULL pointer or exits    */
/* textbin            Converts codon usage to binary data file            */
/* dot(,X)            prints a period every X times it is called          */
/* PrepAFC            Prepare for the COA                                 */
/* DiagoRC            This routine generates the COA                      */
/* colmout            write the output from COA to file                   */
/* rowout             save as above except records the gene information   */
/* inertialig         analyse row inertia and records the results to file */
/* inertiacol         analyse column inertia and record the results       */
/* suprow             add supplementary genes into COA                    */
/* get_aa             converts a three base codon into a 1 or 3 letter AA */
/* codon_error        Called after all codons read, checks data was OK    */
/* rscu_usage_out     Write out RSCU                                      */
/* codon_usage_out    Write out Codon Usage                               */
/* raau_usage_out     Write out normalised amino acid usage               */
/* dinuc_count        Count the dinucleotide usage                        */
/* dinuc_out          Write out dinucleotide usage                        */
/* aa_usage_out       Write out amino acid usage                          */
/* gc_out             Writes various analyses of base usage               */
/* cutab_out          Write a nice tabulation of the RSCU+CU+AA           */
/* base_sil_us_out    Write out base composition at silent sites          */
/* cai_out            Write out CAI usage                                 */
/* cbi_out            Write out codon bias index                          */
/* fop_out            Write out Frequency of Optimal codons               */
/* enc_out            Write out Effective Number of codons                */
/* hydro_out          Write out Protein hydropathicity                    */
/* aromo_out          Write out Protein aromaticity                       */
/* coa_raw_out        Write out raw codon usage for use by COA analysis   */
/*                                                                        */
/*                                                                        */
/* Internal subroutines to Codon.c                                        */
/* my_exit            Controls exit from CodonW closes any open files     */
/* tidy               reads the input data                                */
/* output             called from tidy to decide what to do with the data */
/* file_close         Closes open files                                   */
/*                                                                        */
/**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <stdbool.h>

#include "kseq.h"
KSEQ_INIT(int, read)  

#include "codonW.h"

#if defined(__MWERKS__)
#include <console.h>
#endif

/**************************   MAIN   **************************************/
/* The main function processes commandline arguments to decide analysis to*/
/* run. process_sequence_input() to read in the data files, and count codon usage*/
/* depending on the requested output options output calls various subrou */
/* tines. If COA has been requested it also calls these subroutines and   */
/* recording useful information to summary.coa.                           */
/**************************************************************************/

int main(int argc, char *argv[])
{
  FILE *finput = NULL, *foutput = NULL, *fblkout = NULL;
  FILE *fcoaout = NULL;
  FILE *fsummary = NULL;

  clean_up(ncod, naa); /* zero count of codons and amino acids   */

#if defined(__MWERKS__) /* Macintosh code-warrior */
  argc = ccommand(&argv);
#endif

  pm = &Z_menu;
  pm->my_err = stderr;

  initilize_point(pm->code, pm->f_type, pm->c_type);
  initilize_coa(pm->code);

  fprintf(stderr, "Welcome to CodonW\n");
  proc_comm_line(&argc, &argv);

  finput = pm->inputfile;
  foutput = pm->outputfile;
  fblkout = pm->tidyoutfile;

  fileclose(&pm->fcoa_out);
  if (pm->coa)
    if ((pm->fcoa_out = open_file("coa_raw", "w")) == NULL)
      my_exit(1, "coa_raw");
  fcoaout = pm->fcoa_out;

  fprintf(stderr, "%s", finput);
  num_sequence = process_sequence_input(finput, foutput, fblkout, fcoaout);

  /* num_seq_int_stop value is calculated in codon_usage_out               */
  if (num_seq_int_stop > 0 && pm->warn)
  {
    if (pm->totals && (num_seq_int_stop >= valid_stops))
      fprintf(pm->my_err, "\tWARNING\t At least one sequence in your"
                          " input file has\ninternal stop codons (found %i"
                          " internal stops) \tWARNING\n",
              num_seq_int_stop);
    else
      fprintf(pm->my_err, "\tWARNING\t %i sequences had internal "
                          "stop codons \tWARNING\n",
              num_seq_int_stop);
  }

  if (pm->coa) {
    printf("Generating correspondence analysis\n");
    dot(0, 10);
    run_coa_summary();
  }

  my_exit(0, "");
  return 0;
}

int run_coa_summary()
{
  FILE *finput = NULL, *foutput = NULL, *fblkout = NULL;
  FILE *fcoaout = NULL;
  FILE *fsummary = NULL;

  if ((fsummary = open_file("summary.coa", "w")) == NULL)
    my_exit(1, "summary.coa");
  
  /* set the number of genes in the analysis to the number read in by tidy   */
  pcoa->rows = num_sequence;
  
  /* if COA has been selected then during the reading in phase raw codon usag*/
  /* will have been written to the file coa_raw                              */
  /* text bin converts this to binary data for the COA analysis program      */
  textbin("coa_raw", "cbrawin");

  fprintf(fsummary, "\t\tSummary of Correspondence Analysis \n\n"
                    "The input file was %s it contained %i genes\n"
                    "The number of axes generated was %i\n"
                    "A COA was requested of %s%s usage\n\n\n"
                    "Most of the output presented in this file "
                    "has also been written to separate files\n"
                    "genes.coa\tThe position of the genes on the "
                    "first %i axis\n"
                    "%s.coa\tThe position of the %i %s on the %i "
                    "principle axes\n\n\n",
          pm->curr_infilename,
          pcoa->rows,
          ((pcoa->rows < pcoa->colm) ? pcoa->rows : pcoa->colm) - 1,
          (pm->coa == 'r') ? "relative synonymous " : "",
          (pm->coa == 'a') ? "amino acid" : "codon",
          pcoa->axis,
          (pm->coa == 'a') ? "amino" : "codon",
          pcoa->colm,
          (pm->coa == 'a') ? "amino acids" : "codons",
          pcoa->axis);

  /* allocate memory for the rows and columns, scale both, and write out the*/
  /* resulting matrix to the file cbrawin                                   */
  PrepAFC("cbrawin");

  /* Now do the analysis, calculate the data inertia and all the vectors    */
  DiagoRC(fsummary);

  /* colmout records the position of the columns on each of the factors/axes*/
  if (pm->coa == 'a')
    colmout("cbfcco", "amino.coa", paa, fsummary);
  else
    colmout("cbfcco", "codon.coa", paa, fsummary);

  /* rowout records the position of the genes on each of the axis           */
  rowout("cbfcli", "genes.coa", "coa_raw", fsummary);

  /* pcoa->level == e for exhaustive analysis of inertia                    */
  if (pcoa->level == 'e')
  {
    fprintf(fsummary, "\n\n\nYou requested detailed output from the COA"
                      "\n\nThe absolute and relative inertia "
                      "of each gene and %s (see also inertia.coa)\n",
            (pm->coa == 'a') ? "amino acids" : "codons");
    /* inertialig must preceed inertiacol, records inertia of genes to file   */
    /* it opens the raw codon usage file and loads the raw data to memory     */
    inertialig("inertia.coa", "coa_raw", fsummary);
    /* uses the preloaded raw codon usage, to calculate inertia and other data*/
    /* such as contribution of each column to each factor and to the extent   */
    /* each column is explained by each factor and what the residual variation*/
    /* is                                                                     */
    inertiacol("inertia.coa", fsummary);
  }

  /* if pcoa->add_row is real string, then it will be the name of the file  */
  /* containing additional sequence data, that will be excluded from the COA*/
  /* but factored in, using the original COA vectors and then all other     */
  /* calculation can proceed as with the original data                      */
  if (strlen(pcoa->add_row))
  {
    if ((finput = open_file(pcoa->add_row, "r")) == NULL)
      my_exit(6, "add_row");

    if ((fcoaout = open_file("coa1_raw", "w")) == NULL)
      my_exit(1, "coa1_raw");

    int num_seq = process_sequence_input(finput, NULL, NULL, fcoaout); /* write only fcoaout */

    /* close the files now we are finished                                    */
    fileclose(&fcoaout);
    fileclose(&finput);

    /* covert to binary, use additional raw data file, note not coa_raw this  */
    textbin("coa1_raw", "cb1raw");
    
    /* now call the routine suprow and add these additional genes, we will    */
    /* process this data for inertia and append the gene and col. coordinates */
    /* to the original gene.coa and codon.coa (or amino.coa)                  */
    suprow(num_seq, "cbfcvp", "cb1raw", "genes.coa", "coa1_raw", fsummary);
    fileclose(&fsummary);

  }

  return 0;
}

/* Titles are cleaned up by removing separator if present    */
static int clean_title(char* title, char sep) {
    for (int i = 0; i < (int)strlen(title); i++)
      if (title[i] == sep)
        title[i] = '_';
    return 0;
}

static int strtoupper(char *str) 
{
  for(int i = 0; i < strlen(str); i++)
    str[i] = toupper(str[i]);
  return 0;
}

/**********************  Subroutines     **********************************/
/*  reads input data from a fasta/q formatted sequence file               */
/**************************************************************************/

int process_sequence_input(FILE *finput, FILE *foutput, FILE *fblkout, FILE *fcoaout)
{
	int l;

	kseq_t *seq = kseq_init(fileno(finput));
	while ((l = kseq_read(seq)) >= 0) {
    // strtoupper(seq->seq.s);

    if (pm->totals) /* accumulate sequence          */        
      codon_usage_tot(seq->seq.s);
    else
      print_output(seq->seq.s, seq->name.s, foutput, fblkout, fcoaout);

    num_sequence++;
  }
	kseq_destroy(seq);

  if (pm->totals)
    print_output("", "", foutput, fblkout, fcoaout);

  return (int)num_sequence;
}

// int calculate_features(char *seq, FEATS *features)
// {
//   return 0;
// }

/*************************  output       **********************************/
/* Called from after subroutine tidy has read the sequence into memory    */
/* or  more accurately counted the codon and amino acid usage. This sub-  */
/* routine, via a switch checks which parameters and indices have been    */
/* requested and write these to file, it handles all output except for COA*/
/**************************************************************************/

int print_output(char *seq, char *title, FILE *foutput, FILE *fblkout, FILE *fcoaout)
{
  char sp = pm->separator;

  clean_title(title, pm->separator);

  valid_stops = 0;
  last_aa = codon_usage_tot(seq);
  /* codon_error, if 4th parameter is 1, then checks for valid start and  */
  /* internal stop codon, if 4th parmater is 2, checks that the last codon*/
  /* is a stop or was partial, and for non-translatable codons            */
  codon_error(last_aa, valid_stops, title, (char)1);
  codon_error(last_aa, valid_stops, title, (char)2);

  /* if we are concatenating sequences then change the title to avger_of  */
  if (pm->totals)
    strcpy(title, "Average_of_genes");

  if (strchr("OCASDLBX", (int)pm->bulk) != NULL)
  {
    switch ((int)pm->bulk)
    {
    case 'S':
      rscu_usage_out(fblkout, ncod, naa);
      break;
    case 'C':
      codon_usage_out(fblkout, ncod, last_aa, valid_stops, title);
      break;
    case 'L':
      raau_usage_out(fblkout, naa);
      break;
    case 'D':
      dinuc_count(seq, tot);
      dinuc_out(fblkout, title);
      break;
    case 'A':
      aa_usage_out(fblkout, naa);
      break;
    case 'B':
      gc_out(foutput, fblkout, 1);
      break;
    case 'O':
      cutab_out(fblkout, ncod, naa);
      break;
    }
  }

  /* if an index has been requested then this is true                     */
  if (pm->sil_base || pm->cai || pm->fop || pm->enc || pm->gc3s ||
      pm->gc || pm->cbi || pm->L_sym || pm->L_aa || pm->coa ||
      pm->hyd || pm->aro)
  {
    /* if this is the first sequence then write a header line           */

    if (num_sequence == 0 || pm->totals)
    {

      fprintf(foutput, "%-.25s%c", "title", sp);
      if (pm->sil_base)
        fprintf(foutput, "%s%c%s%c%s%c%s%c", "T3s", sp, "C3s", sp, "A3s", sp,
                "G3s", sp);
      if (pm->cai)
        fprintf(foutput, "%s%c", "CAI", sp);
      if (pm->cbi)
        fprintf(foutput, "%s%c", "CBI", sp);
      if (pm->fop)
        fprintf(foutput, "%s%c", "Fop", sp);
      if (pm->enc)
        fprintf(foutput, "%s%c", "Nc", sp);
      if (pm->gc3s)
        fprintf(foutput, "%s%c", "GC3s", sp);
      if (pm->gc)
        fprintf(foutput, "%s%c", "GC", sp);
      if (pm->L_sym)
        fprintf(foutput, "%s%c", "L_sym", sp);
      if (pm->L_aa)
        fprintf(foutput, "%s%c", "L_aa", sp);
      if (pm->hyd)
        fprintf(foutput, "%s%c", "Gravy", sp);
      if (pm->aro)
        fprintf(foutput, "%s%c", "Aromo", sp);

      fprintf(foutput, "\n");
    }

    /* if output format is human readable print the fixed width sequence  */
    /* name, else print only the name of the sequence                     */
    fprintf(foutput, "%-.25s%c", title, sp);

    /*Need to use if statements as we allow more than one index to be calc*/
    /* per sequence read in                                               */
    if (pm->sil_base)
      base_sil_us_out(foutput, ncod, naa);
    if (pm->cai)
      cai_out(foutput, ncod);
    if (pm->cbi)
      cbi_out(foutput, ncod, naa);
    if (pm->fop)
      fop_out(foutput, ncod);
    if (pm->enc)
      enc_out(foutput, ncod, naa);
    if (pm->gc3s)
      gc_out(foutput, fblkout, 3);
    if (pm->gc)
      gc_out(foutput, fblkout, 2);
    if (pm->L_sym)
      gc_out(foutput, fblkout, 4);
    if (pm->L_aa)
      gc_out(foutput, fblkout, 5);
    if (pm->hyd)
      hydro_out(foutput, naa);
    if (pm->aro)
      aromo_out(foutput, naa);
    if (pm->coa)
      coa_raw_out(fcoaout, ncod, naa, title);

    fprintf(foutput, "\n");
  }

  clean_up(ncod, naa);

  return 0;
}

/************************* my_exit       **********************************/
/* Called to clean up open files and generate an intelligent exit message */
/* Also warns if no analysis has been run, the user did not select R from */
/* the main menu. If COA was selected then it reminds the user to look    */
/* at the file summary.coa, and deletes any stray binary files            */
/**************************************************************************/

int my_exit(int error_num, char *message)
{

  fileclose(&pm->inputfile);
  fileclose(&pm->outputfile);
  fileclose(&pm->tidyoutfile);

  fileclose(&pm->cuout);
  fileclose(&pm->fopfile);
  fileclose(&pm->cbifile);
  fileclose(&pm->caifile);
  fileclose(&pm->logfile);
  fileclose(&pm->fcoa_in);
  fileclose(&pm->fcoa_out);

  if (pm->inputfile = fopen("cbrawin", "r"))
  {
    fclose(pm->inputfile);
    remove("cbrawin");
  }
  if (pm->inputfile = fopen("cbfcco", "r"))
  {
    fclose(pm->inputfile);
    remove("cbfcco");
  }
  if (pm->inputfile = fopen("cbfcli", "r"))
  {
    fclose(pm->inputfile);
    remove("cbfcli");
  }
  if (pm->inputfile = fopen("cbfcpc", "r"))
  {
    fclose(pm->inputfile);
    remove("cbfcpc");
  }
  if (pm->inputfile = fopen("cbfcpl", "r"))
  {
    fclose(pm->inputfile);
    remove("cbfcpl");
  }
  if (pm->inputfile = fopen("cbfcta", "r"))
  {
    fclose(pm->inputfile);
    remove("cbfcta");
  }
  if (pm->inputfile = fopen("cbfcvp", "r"))
  {
    fclose(pm->inputfile);
    remove("cbfcvp");
  }
  if (pm->inputfile = fopen("cb1rawin", "r"))
  {
    fclose(pm->inputfile);
    remove("cb1rawin");
  }

  switch ((int)error_num)
  {

  case 0:
    /* silent exit */
    exit(0);
    break;
  case 1:
    fprintf(stderr, "failed to open file for input/output <%s>\n", message);
    exit(1);
    break;
  case 2:
    printf("user requested exit <%s>\n", message);
    exit(0);
    break;
  case 3:
    printf("failed to allocate memory <%s>\n", message);
    exit(1);
    break;
  case 4:
    printf("Write to disk failed ! <%s>\n", message);
    exit(1);
    break;
  case 5:
    printf("Read from disk failed! <%s>\n", message);
    exit(1);
    break;
  case 6:
    printf("failed to open file for reading <%s>\n", message);
    exit(1);
    break;
  case 7:
    printf("failed to close file <%s>\n", message);
    exit(1);
  case 99:
    printf(" Controlled exit <%s>\n", message);
    exit(0);
    break;
  default:
    printf("for unknown reason\n");
    exit(1);
    break;
  }
  return 0;
}

/************************** file_close   **********************************/
/* Fileclose function checks whether the filepointer is open, if so it    */
/* attempts to close the open file handle and assigns a null pointer      */
/* to that handle                                                         */
/**************************************************************************/

int fileclose(FILE **file_pointer)
{
  if (*file_pointer != NULL)
  {
    if (fclose(*file_pointer) == EOF)
    {
      fprintf(stderr, "Failed to close file %i \n", errno);
      perror("Unexpected condition in fileclose");
      exit(7);
    }
    *file_pointer = NULL; /* make sure file_pointer is null*/
  }
  return 1;
}


/************** open_file **************************************************/
/* This subroutine is a front end to fopen open.                           */
/***************************************************************************/

FILE *open_file(char *file_needed, char *write_perm)
{
    FILE *input = NULL;

    if (!(input = fopen(file_needed, write_perm)))
        my_exit(2, "File not found");

    return input;
}
