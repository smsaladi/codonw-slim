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
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "codonW.h"

/************** process_command_line  *************************************/
/* The command line is passed to this function for processing. The name of*/
/* the programme is read, and based on this, CodonW will emulate several  */
/* useful codon usage analysis programmes routinely used in our laboratory*/
/* all other command line arguments are read. Unrecognised arguments are  */
/* reported to the user, arguments not preceded by a dash are assumed to  */
/* be filenames. The input, output and bulk output files to be precise    */
/**************************************************************************/
int proc_comm_line(int *pargc, char ***pargv)
{
    char *p;
    char c;
    int n;
    char root[MAX_FILENAME_LEN];

    /* first call to garg initialises the function  with the command line*/
    /* parameters and the number of arguments, subsequent calls strip    */
    /* these off one by one                                              */

    /* has the user asked for help               ????????????            */
    if ((p = garg(*pargc, *pargv, "-h", GARG_EXACT)) ||
        (p = garg(0, NULL, "-help", GARG_EXACT)))
    {
        printf(
            "codonW [inputfile] [outputfile] [bulkoutfile] [options]\n"
            "General options and defaults:\n"
            " -h(elp)\tThis help message\n"
            " -nowarn\tPrevent warnings about sequences being displayed\n"
            " -silent\tOverwrite files silently\n"
            " -totals\tConcatenate all genes in inputfile\n"
            " -code N\tGenetic code as defined under menu 3 option 5\n"
            " -f_type N\tFop/CBI codons as defined by menu 3 option 6\n"
            " -c_type N\tCai fitness values as defined by menu 3 option 7\n"
            " -t (char)\tColumn separator to be used in output files "
            "(comma,tab,space)\n"
            "\nCodon usage indices and Amino acid indices \n"
            " -cai\t\tcalculate Codon Adaptation Index (CAI)\n"
            " -fop\t\tcalculate Frequency of OPtimal codons index (FOP)\n"
            " -cbi\t\tcalculate Codon Bias Index (CBI)\n"
            " -enc\t\tEffective Number of Codons (ENc)\n"
            " -gc\t\tG+C content of gene (all 3 codon positions)\n"
            " -gcs3\t\tGC of synonymous codons 3rd positions\n"
            " -sil_base\tBase composition at synonymous third codon "
            "positions\n"
            " -L_sym\t\tNumber of synonymous codons\n"
            " -L_aa\t\tTotal number of synonymous and non-synonymous codons\n"
            " -all_indices\t\tAll the above indices\n"
            " -aro\t\tCalculate aromaticity of protein\n"
            " -hyd\t\tCalculate hydropathicity of protein\n"
            " -cai_file  {file}\tUser input file of CAI values\n"
            " -cbi_file  {file}\tUser input file of CBI values\n"
            " -fop_file  {file}\tUser input file of Fop values\n"
            "\nCorrespondence analysis (COA) options \n"
            " -coa_cu \tCOA of codon usage frequencies\n"
            " -coa_rscu\tCOA of Relative Synonymous Codon Usage\n"
            " -coa_aa\tCOA of amino acid usage frequencies\n"
            " -coa_expert\tGenerate detailed(expert) statistics on COA\n"
            " -coa_axes N\tSelect number of axis to record\n"
            " -coa_num N\tSelect number of genes to use to identify "
            "optimal codons\n"
            "\t\tvalues can be whole numbers or a percentage (5 or 10%%)\n"
            "\nBulk output options | only one can be selected per analysis\n"
            " -aau\t\tAmino Acid Usage (AAU)\n"
            " -raau\t\tRelative Amino Acid Usage (RAAU)\n"
            " -cu\t\tCodon Usage (CU) (default)\n"
            " -cutab\t\tTabulation of codon usage\n"
            " -cutot\t\tTabulation of dataset's codon usage\n"
            " -rscu\t\tRelative Synonymous Codon Usage (RSCU)\n"
            " -fasta\t\tfasta format\n"
            " -tidy\t\tfasta format\n"
            " -reader\tReader format (codons are separated by spaces)\n"
            " -transl\tConceptual translation of DNA to amino acid\n"
            " -base\t\tDetailed report of codon G+C composition\n"
            " -dinuc\t\tDinucleotide usage of the three codon pos.\n"
            " -noblk\t\tNo bulk output to be written to file\n"
            "\nWhere {file} represents an input filename, and N an integer"
            " value\n");
        my_exit(0, ""); /* after writing out help quit         */
    }

    /* -silent stops warnings about file about to be overwritten               */
    if (garg(0, NULL, "-silent", GARG_THERE))
        pm->verbose = FALSE;

    /* -total causes sequences to be concatenated and treated as one sequence */
    if (garg(0, NULL, "-total", GARG_THERE))
        pm->totals = TRUE;

    /* -code determines the genetic code                                       */
    if (p = garg(0, NULL, "-code", GARG_NEXT | GARG_EXACT))
    {
        strcpy(pm->junk, p);
        n = 0;
        while (isdigit((int)pm->junk[n]) && pm->junk[n] != '\0')
            n++;
        if (n != (int)strlen(pm->junk) || atoi(pm->junk) < 0 || atoi(pm->junk) > NumGeneticCodes)
        {
            printf("FATAL: The value for genetic code %s is invalid\n",
                   pm->junk);
            my_exit(99, "Fatal error in genetic code value");
        }
        else
        {
            pm->code = (char)atoi(p); /* define genetic code */
            initilize_point(pm->code, pm->f_type, pm->c_type);
        }
    }

    /* -f_type selects which of the predefined fop values to use               */
    /* NB. The fop is selected with the integer value corresponding to the menu*/
    /* choice under the defaults menu. It must be in the range 1-NumFopSpecies */

    if (p = garg(0, NULL, "-f_type", GARG_NEXT | GARG_EXACT))
    {
        strcpy(pm->junk, p);
        n = 0;
        while (isdigit((int)pm->junk[n]) && pm->junk[n] != '\0')
            n++;
        if (n != (int)strlen(pm->junk) || atoi(pm->junk) < 0 ||
            atoi(pm->junk) >= NumFopSpecies)
        {
            printf("FATAL: The value for fop_type %s is not valid\n",
                   pm->junk);
            my_exit(99, "Fatal error in Fop value");
        }
        else
        {
            pm->f_type = (char)atoi(p); /* define organism type for Fop  */
            initilize_point(pm->code, pm->f_type, pm->c_type);
        }
    }

    /* -d_type selects which of the predefined CAI values to use               */
    /* NB. The CAI is selected with the integer value corresponding to the menu*/
    /* choice under the defaults menu. It must be in the range 1-NumCAISpecies */
    if (p = garg(0, NULL, "-c_type", GARG_NEXT | GARG_EXACT))
    {
        strcpy(pm->junk, p);
        n = 0;
        while (isdigit((int)pm->junk[n]) && pm->junk[n] != '\0')
            n++;
        if (n != (int)strlen(pm->junk) || atoi(pm->junk) < 0 ||
            atoi(pm->junk) >= NumCaiSpecies)
        {
            printf("FATAL: The value for cai_type %s is not valid\n",
                   pm->junk);
            my_exit(99, "Fatal error in CAI type value");
        }
        else
        {
            pm->c_type = (char)atoi(p); /* define organism type for CAI  */
            initilize_point(pm->code, pm->f_type, pm->c_type);
        }
    }

    /* Command line arguments for the indices menu (4)                        */
    /* The presence of any of these flags, cause the relevant indices to be   */
    /* calculated                                                             */
    /* Indices are CAI, FOP, CBI, Nc, GC, GC3s, Lsyn, Laa, silent_base        */
    /* composition, hydropathicity, aromaticity                               */
    if (p = garg(0, NULL, "-cai", GARG_EXACT))
        pm->cai = TRUE;
    if (p = garg(0, NULL, "-fop", GARG_EXACT))
        pm->fop = TRUE;
    if (p = garg(0, NULL, "-cbi", GARG_EXACT))
        pm->cbi = TRUE;
    if (p = garg(0, NULL, "-enc", GARG_EXACT))
        pm->enc = TRUE;
    if (p = garg(0, NULL, "-gc", GARG_EXACT))
        pm->gc = TRUE;
    if (p = garg(0, NULL, "-gc3s", GARG_EXACT))
        pm->gc3s = TRUE;
    if (p = garg(0, NULL, "-sil_base", GARG_EXACT))
        pm->sil_base = TRUE;
    if (p = garg(0, NULL, "-L_sym", GARG_EXACT))
        pm->L_sym = TRUE;
    if (p = garg(0, NULL, "-L_aa", GARG_EXACT))
        pm->L_aa = TRUE;
    if (p = garg(0, NULL, "-hyd", GARG_EXACT))
        pm->hyd = TRUE;
    if (p = garg(0, NULL, "-aro", GARG_EXACT))
        pm->aro = TRUE;
    /* Turns on all the above indices                                         */
    if (p = garg(0, NULL, "-all_indices", GARG_EXACT))
    {
        pm->cai = TRUE;
        pm->fop = TRUE;
        pm->cbi = TRUE;
        pm->enc = TRUE;
        pm->gc = TRUE;
        pm->gc3s = TRUE;
        pm->sil_base = TRUE;
        pm->L_sym = TRUE;
        pm->L_aa = TRUE;
        pm->hyd = TRUE;
        pm->aro = TRUE;
    }

    /* This section in used to input the filenames for personal choices of Fop */
    /* CBI or CAI values. The name is tested to make sure the file is readable */
    /* the pointer to the file is then assign to the relevant pointer in the   */
    /* struct Z_menu and then processed properly in codon_us.c                 */

    /* Fop                                                                     */
    if (p = garg(0, NULL, "-fop_file", GARG_NEXT | GARG_EXACT))
    {
        if ((pm->fopfile = open_file("", p, "r", FALSE)) == NULL)
        {
            printf("Could not open Fop file - %s\n", p);
            my_exit(1, "commline open fop file");
        }
        else
            strncpy(pm->fop_filen, pm->junk, MAX_FILENAME_LEN - 1);
        /* idiot catch, if you load personal fop values you want to calculate fop */
        pm->fop = TRUE;
    }

    /* CAI                                                                     */
    if (p = garg(0, NULL, "-cai_file", GARG_NEXT | GARG_EXACT))
    {
        if ((pm->caifile = open_file("", p, "r", FALSE)) == NULL)
        {
            printf("Could not open CAI file - %s\n", p);
            my_exit(1, "commline failed error");
        }
        else
            strncpy(pm->cai_filen, pm->junk, MAX_FILENAME_LEN - 1);
        pm->cai = TRUE; /* idiot catch          */
    }
    /* CBI                                                                     */
    if (p = garg(0, NULL, "-cbi_file", GARG_NEXT | GARG_EXACT))
    {
        if ((pm->cbifile = open_file("", p, "r", FALSE)) == NULL)
        {
            printf("Could not open CBI file - %s\n", p);
            my_exit(1, "Commline failed to open file");
        }
        else
            strncpy(pm->cbi_filen, pm->junk, MAX_FILENAME_LEN - 1);
        pm->cbi = TRUE; /* idiot catch            */
    }

    /* This section changes the default correspondence menu choices normally   */
    /* set in menu menu 5.                                                     */
    /* Note only one of -coa_cu -coa_rscu -coa_aa can be chosen                */
    if (p = garg(0, NULL, "-coa_cu", GARG_EXACT))
        pm->coa = 'c';
    if (p = garg(0, NULL, "-coa_rscu", GARG_EXACT))
        pm->coa = 'r';
    if (p = garg(0, NULL, "-coa_aa", GARG_EXACT))
        pm->coa = 'a';
    if (p = garg(0, NULL, "-coa_expert", GARG_EXACT)) /* detailed inertia */
        (coa.level = 'e');                            /* analysis         */

    /* These are options selectable under the advanced COA menu                */
    /* This first option -coa_axes changes the number of axis recorded to file */
    if (p = garg(0, NULL, "-coa_axes", GARG_NEXT | GARG_EXACT))
    {
        if (isdigit((int)*p))
        {
            n = (char)atoi(p);
            /* just check that correspondence analysis has been selected         */
            if (pm->coa == 'a' && (n > 20 || n < 0) || (n < 0 || n > 59))
            {
                fprintf(pm->my_err, "Value %d is out of range for Number COA Axis "
                                    "adjusting to max value\n",
                        n);
                if (pm->coa == 'a')
                    pcoa->axis = 20;
                else
                    pcoa->axis = 59;
            }
            else
            {
                pcoa->axis = (char)n;
            }
        }
    }

    /* Select the size of dataset to use to identify optimal codons            */
    if (p = garg(0, NULL, "-coa_num", GARG_NEXT | GARG_EXACT))
    {
        strcpy(pm->junk, p);
        if ((p = strchr(pm->junk, '%')) != NULL)
        {
            p = '\0';
            pcoa->fop_gene = atoi(pm->junk) * -1;
        }
        else
        {
            pcoa->fop_gene = atoi(pm->junk);
        }
    }

    /* These option are mutually exclusive and are normally selected using the */
    /* the bulk output menu (menu 8)                                           */

    if (p = garg(0, NULL, "-raau", GARG_EXACT))
        pm->bulk = 'L';
    if (p = garg(0, NULL, "-cu", GARG_EXACT))
        pm->bulk = 'C';
    if (p = garg(0, NULL, "-cutab", GARG_THERE))
        pm->bulk = 'O';
    if (p = garg(0, NULL, "-cutot", GARG_THERE))
    {
        pm->bulk = 'C';
        pm->totals = TRUE;
    }
    if (p = garg(0, NULL, "-reader", GARG_EXACT))
        pm->bulk = 'R';
    if (p = garg(0, NULL, "-rscu", GARG_EXACT))
        pm->bulk = 'S';
    if (p = garg(0, NULL, "-tidy", GARG_EXACT))
        pm->bulk = 'T';
    if (p = garg(0, NULL, "-fasta", GARG_EXACT))
        pm->bulk = 'T';
    if (p = garg(0, NULL, "-aau", GARG_EXACT))
        pm->bulk = 'A';
    if (p = garg(0, NULL, "-transl", GARG_THERE))
        pm->bulk = 'N';
    if (p = garg(0, NULL, "-base", GARG_THERE))
        pm->bulk = 'B';
    if (p = garg(0, NULL, "-dinuc", GARG_THERE))
        pm->bulk = 'D';
    if (p = garg(0, NULL, "-noblk", GARG_EXACT))
        pm->bulk = 'X';

    /* -t is used to change the column separator used in the output files     */
    /* at present it must be a space, tab or comma                            */
    /* Must occur after -transl or it misreads transl as a seperator          */
    if (p = garg(0, NULL, "-t", GARG_NEXT | GARG_SUBSQ))
    {
        strcpy(pm->junk, p);
        n = 0;
        do
        {
            c = pm->junk[n++];
        } while (strchr("'\"\0", (int)c) != NULL);
        if (strchr("\t, ", (int)c) == NULL)
        {
            printf("WARNING: The chosen separator %s is unsuitable use"
                   "comma, tab or space\n",
                   pm->junk);
        }
        else
        {
            pm->seperator = c;
        }
    }

    pm->codonW = TRUE;
    if (pm->bulk == 'X')
        pm->bulk = 'C';

    /* By this point we should have processed all the command line arguments  */
    /* so now we test for any remaining, these are unrecognised               */

    while ((p = garg(0, NULL, "-", GARG_THERE)))
        if (strcmp(&p, "-"))
            continue;
        else {
            /* if we are running without a menu then abort this run         */
            sprintf(pm->junk, "Unrecognised argument %s", p);
            my_exit(99, pm->junk);
        }

    /* Anything remaining should be file names                                */
    /* The first name should be the input file name                           */

    if (p = garg(0, NULL, "", GARG_THERE))
    {
        if ((pm->inputfile = open_file("", p, "r", FALSE)) == NULL)
        {
            printf("Could not open input file - %s\n", p);
            my_exit(1, "failed to open file in proc_commline");
        }
        else
            strncpy(pm->curr_infilename, pm->junk, MAX_FILENAME_LEN - 1);
    }
    /* The second should be the output filename                               */
    if (p = garg(0, NULL, "", GARG_THERE))
    {
        if (strcmp(&p, "-"))
            pm->outputfile = stdout;
        else if ((pm->outputfile = open_file("", p, "w",
                                        (int)pm->verbose)) == NULL)
        {
            printf("Could not open output file - %s\n", p);
            my_exit(1, "commline out file");
        }
        else
            strncpy(pm->curr_outfilename, pm->junk, MAX_FILENAME_LEN - 1);
    }

    /* The third which only occurs if the programme is running as CodonW      */
    if (pm->codonW && (p = garg(0, NULL, "", GARG_THERE)))
    {
        if ((pm->tidyoutfile = open_file("", p, "w",
                                         (int)pm->verbose)) == NULL)
        {
            printf("Could not open blkoutput file - %s\n", p);
            my_exit(1, "commline blk outfile");
        }
        else
            strncpy(pm->curr_tidyoutname, pm->junk, MAX_FILENAME_LEN - 1);
    }

    /* Now check the command line is empty ... it should be at this point     */
    while (p = garg(0, NULL, "", GARG_THERE))
        printf("This command line parameter was not recognised %s\n", p);

    /* IF no file name was found on the command line and the programme is     */
    /* impersonating another programme or we decided not to use the menu      */
    /* we need to load an input file name                                     */

    if (!pm->inputfile)
        my_exit(1, "Please provide inputfile as argument");
    if (!pm->outputfile)
        my_exit(1, "Please provide outputfile as argument");
    if (!pm->tidyoutfile)
        my_exit(1, "Please provide bulk file as argument");

    return 1;
}
/****************** Garg     ***********************************************/
/* This subroutine strips of the commandline arguments and passes them back*/
/* to the calling function. Each time it is called with argc and argv non  */
/* null the commandline is refreshed. If called with these are null args   */
/* a commandline pre-stored is used, this commandline is striped arg by arg*/
/* as they are identified                                                  */
/* This subroutine was developed as a collaboration with Colin McFarlane   */
/* GARG_EXACT           The argument must match targ exactly               */
/* GARG_THERE           The targ may be sub-string of the argument         */
/* GARG_SUBSQ           The string immediate after targ is returned        */
/* GARG_NEXT            The next argument after targ is returned           */
/* else                 return NULL                                        */
/***************************************************************************/
char *garg(int argc, char *argv[], const char *targ, int mode)
{
    static char *argw[MAX_ARGS];
    static int done[MAX_ARGS];
    static int argn;

    int arg = 1, nc;

    if (argv)
    {
        if (--argc < 1)
            return NULL;
        for (argn = 0; argn < argc; argn++)
        {
            argw[argn] = argv[argn + 1];
            done[argn] = 0;
        }
    }
    nc = mode & GARG_EXACT ? BUFSIZ : strlen(targ);

    for (arg = 0; arg < argn; arg++)
        if ((0 == strncmp(targ, argw[arg], nc)) && !done[arg])
        {
            done[arg] = 1;
            if (mode & GARG_THERE)
                return argw[arg];
            if (mode & GARG_SUBSQ)
                return &argw[arg][nc];
            if (mode & GARG_NEXT)
            {
                done[++arg < argn ? arg : --arg] = 1;
                return argw[arg];
            }
            return argw[arg];
        }
    return NULL;
}
