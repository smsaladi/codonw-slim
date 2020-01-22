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

#include "../include/codonW.h"

/* define genetic codes   */
GENETIC_CODE_STRUCT cu_ref[] = {
    {
        "Universal Genetic code",
        "TGA=* TAA=* TAG=*",
        {
            0,
            1, 6, 10, 18, 1, 6, 10, 18, 2, 6, 11, 11, 2, 6, 11, 19,
            2, 7, 12, 20, 2, 7, 12, 20, 2, 7, 13, 20, 2, 7, 13, 20,
            3, 8, 14,  6, 3, 8, 14,  6, 3, 8, 15, 20, 4, 8, 15, 20,
            5, 9, 16, 21, 5, 9, 16, 21, 5, 9, 17, 21, 5, 9, 17, 21
        }
    },
    {
        "Vertebrate Mitochondrial code",
        "AGR=* ATA=M TGA=W",
        {
            0,
            1, 6, 10, 18, 1, 6, 10, 18, 2, 6, 11, 19, 2, 6, 11, 19,
            2, 7, 12, 20, 2, 7, 12, 20, 2, 7, 13, 20, 2, 7, 13, 20,
            3, 8, 14,  6, 3, 8, 14,  6, 4, 8, 15, 11, 4, 8, 15, 11,
            5, 9, 16, 21, 5, 9, 16, 21, 5, 9, 17, 21, 5, 9, 17, 21
        }
    },
    {
        "Yeast Mitochondrial code",
        "CTN=* ATA=M TGA=W",
        {
            0,
            1, 6, 10, 18, 1, 6, 10, 18, 2, 6, 11, 19, 2, 6, 11, 19,
            8, 7, 12, 20, 8, 7, 12, 20, 8, 7, 13, 20, 8, 7, 13, 20,
            3, 8, 14,  6, 3, 8, 14,  6, 4, 8, 15, 20, 4, 8, 15, 20,
            5, 9, 16, 21, 5, 9, 16, 21, 5, 9, 17, 21, 5, 9, 17, 21
        }
    },
    {
        "Filamentous fungi Mitochondrial code",
        "TGA=W",
        {
            0,
            1, 6, 10, 18, 1, 6, 10, 18, 2, 6, 11, 19, 2, 6, 11, 19,
            2, 7, 12, 20, 2, 7, 12, 20, 2, 7, 13, 20, 2, 7, 13, 20,
            3, 8, 14,  6, 3, 8, 14,  6, 3, 8, 15, 20, 4, 8, 15, 20,
            5, 9, 16, 21, 5, 9, 16, 21, 5, 9, 17, 21, 5, 9, 17, 21
        }
    },
    {
        "Insects and Plathyhelminthes Mitochondrial code",
        "ATA=M TGA=W AGR=S",
        {
            0,
            1, 6, 10, 18, 1, 6, 10, 18, 2, 6, 11, 19, 2, 6, 11, 19,
            2, 7, 12, 20, 2, 7, 12, 20, 2, 7, 13, 20, 2, 7, 13, 20,
            3, 8, 14,  6, 3, 8, 14,  6, 4, 8, 15,  6, 4, 8, 15,  6,
            5, 9, 16, 21, 5, 9, 16, 21, 5, 9, 17, 21, 5, 9, 17, 21
        }
    },
    {
        "Nuclear code of Cilitia",
        "UAA=Q=Gln  UAG=Q",
        {
            0,
            1, 6, 10, 18, 1, 6, 10, 18, 2, 6, 13, 11, 2, 6, 13, 19,
            2, 7, 12, 20, 2, 7, 12, 20, 2, 7, 13, 20, 2, 7, 13, 20,
            3, 8, 14,  6, 3, 8, 14,  6, 3, 8, 15, 20, 4, 8, 15, 20,
            5, 9, 16, 21, 5, 9, 16, 21, 5, 9, 17, 21, 5, 9, 17, 21
        }
    },
        {"Nuclear code of Euplotes",
        "UGA=C",
        {
            0,
            1, 6, 10, 18, 1, 6, 10, 18, 2, 6, 11, 18, 2, 6, 11, 19,
            2, 7, 12, 20, 2, 7, 12, 20, 2, 7, 13, 20, 2, 7, 13, 20,
            3, 8, 14,  6, 3, 8, 14,  6, 3, 8, 15, 20, 4, 8, 15, 20,
            5, 9, 16, 21, 5, 9, 16, 21, 5, 9, 17, 21, 5, 9, 17, 21
        }
    },
    {
        "Mitochondrial code of Echinoderms",
        "UGA=W AGR=S AAA=N",
        {
            0,
            1, 6, 10, 18, 1, 6, 10, 18, 2, 6, 11, 19, 2, 6, 11, 19,
            2, 7, 12, 20, 2, 7, 12, 20, 2, 7, 13, 20, 2, 7, 13, 20,
            3, 8, 14,  6, 3, 8, 14,  6, 3, 8, 14,  6, 4, 8, 15,  6,
            5, 9, 16, 21, 5, 9, 16, 21, 5, 9, 17, 21, 5, 9, 17, 21
        }
    }
};

/* define amino acid info     */
AMINO_STRUCT amino_acids = {
    {
        "X",
        "F", "L", "I", "M", "V",
        "S", "P", "T", "A", "Y",
        "*", "H", "Q", "N", "K",
        "D", "E", "C", "W", "R", "G"
    },
    {
        "UNK",
        "Phe", "Leu", "Ile", "Met", "Val",
        "Ser", "Pro", "Thr", "Ala", "Tyr",
        "TER", "His", "Gln", "Asn", "Lys",
        "Asp", "Glu", "Cys", "Trp", "Arg", "Gly"
    },
    {
        "BAD",
        "UUU", "UCU", "UAU", "UGU", "UUC", "UCC", "UAC", "UGC", "UUA", "UCA", "UAA", "UGA", "UUG", "UCG", "UAG", "UGG",
        "CUU", "CCU", "CAU", "CGU", "CUC", "CCC", "CAC", "CGC", "CUA", "CCA", "CAA", "CGA", "CUG", "CCG", "CAG", "CGG",
        "AUU", "ACU", "AAU", "AGU", "AUC", "ACC", "AAC", "AGC", "AUA", "ACA", "AAA", "AGA", "AUG", "ACG", "AAG", "AGG",
        "GUU", "GCU", "GAU", "GGU", "GUC", "GCC", "GAC", "GGC", "GUA", "GCA", "GAA", "GGA", "GUG", "GCG", "GAG", "GGG"
    }
};

FOP_STRUCT fop_ref[] = {
    {
        "Escherichia coli",
        "Ikemura (1985) Mol. Biol. Evol. 2:13-34 (updated by INCBI 1991)",
        {
            0,
            2, 3, 2, 2, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 3, 2, 2, 3, 3, 2, 2, 2, 2, 3, 3, 3, 2,
            2, 3, 2, 2, 3, 3, 3, 3, 2, 2, 3, 2, 2, 2, 2, 2,
            3, 3, 2, 3, 2, 2, 3, 3, 2, 2, 3, 2, 2, 3, 2, 2
        }
    },
    {
        "Bacillus subtilis",
        "Sharp et al (1990) Genetics & Biotech of Bacilli vol3 pp89-98",
        {
            0,
            2, 3, 2, 2, 3, 1, 3, 2, 2, 2, 2, 2, 2, 1, 2, 2,
            3, 3, 2, 3, 2, 1, 2, 3, 2, 3, 3, 1, 2, 2, 2, 1,
            2, 3, 2, 2, 3, 1, 3, 2, 1, 2, 3, 2, 2, 2, 2, 1,
            3, 3, 2, 3, 2, 1, 3, 2, 3, 2, 3, 2, 2, 2, 2, 1
        }
    },
    {
        "Dictyostelium discoideum",
        "Sharp and Devine (1989) Nucl. Acids Res 17:5029-5039)",
        {
            0,
            2, 2, 2, 2, 3, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 3, 3, 2, 3, 2, 2, 3, 3, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 3, 2,
            2, 2, 2, 3, 3, 3, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2
        }
    },
    {
        "Aspergillus nidulans",
        "Lloyd and Sharp (1991) Mol. Gen. Genet 230: 288-294",
        {
            0,
            2, 2, 2, 2, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 3, 2,
            2, 2, 2, 2, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 3, 2,
            2, 3, 2, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 3, 2
        }
    },
    {
        "Saccharomyces cerevisiae",
        "Sharp and Cowe (1991) Yeast 7:657-678",
        {
            0,
            2, 3, 2, 3, 3, 3, 3, 2, 2, 2, 2, 2, 3, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 3, 2, 2, 3, 3, 2, 2, 2, 2, 2,
            3, 3, 2, 2, 3, 3, 3, 2, 2, 2, 2, 3, 2, 2, 3, 2,
            3, 3, 2, 3, 3, 2, 3, 2, 2, 2, 3, 2, 2, 2, 2, 2
        }
    },
    {
        "Drosophila melanogaster",
        "Shields et al. (1988) Mol Biol Evol 5: 704-716",
        {
            0,
            2, 2, 2, 2, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 3, 2, 3, 3, 3, 2, 2, 2, 2, 3, 2, 3, 2,
            2, 2, 2, 2, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 3, 2,
            2, 2, 2, 2, 3, 3, 3, 3, 2, 2, 2, 2, 3, 2, 3, 2
        }
    },
    {
        "Caenorhabditis elegans",
        "Stenico, Lloyd and Sharp Nuc. Acids Res. 22: 2437-2446(1994)",
        {
            0,
            2, 2, 2, 2, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2,
            3, 2, 2, 3, 3, 2, 3, 3, 2, 3, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 3, 2,
            2, 3, 2, 2, 3, 3, 3, 2, 2, 2, 2, 3, 2, 2, 3, 2
        }
    },
    {
        "Neurospora crassa",
        "Lloyd and Sharp (1993)",
        {
            0,
            2, 3, 2, 2, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 3, 2,
            2, 3, 2, 2, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 3, 2,
            2, 2, 2, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 3, 2
        }
    }
};

CAI_STRUCT cai_ref[] = {  
    {
        "Escherichia coli",
        "No reference",
        {
            0.000F,
            0.296F, 1.000F, 0.239F, 0.500F, 1.000F, 0.744F, 1.000F, 1.000F,
            0.020F, 0.077F, 0.000F, 0.000F, 0.020F, 0.017F, 0.000F, 1.000F,
            0.042F, 0.070F, 0.291F, 1.000F, 0.037F, 0.012F, 1.000F, 0.356F,
            0.007F, 0.135F, 0.124F, 0.004F, 1.000F, 1.000F, 1.000F, 0.004F,
            0.185F, 0.965F, 0.051F, 0.085F, 1.000F, 1.000F, 1.000F, 0.410F,
            0.003F, 0.076F, 1.000F, 0.004F, 1.000F, 0.099F, 0.253F, 0.002F,
            1.000F, 1.000F, 0.434F, 1.000F, 0.066F, 0.122F, 1.000F, 0.724F,
            0.495F, 0.586F, 1.000F, 0.010F, 0.221F, 0.424F, 0.259F, 0.019F
        }
    },
    {
        "Bacillus subtilis",
        "No reference",
        {
            0.00F,
            0.571F, 1.000F, 0.500F, 1.000F, 1.000F, 0.021F, 1.000F, 1.000F,
            1.000F, 0.458F, 0.000F, 0.000F, 0.036F, 0.021F, 0.000F, 1.000F,
            0.857F, 1.000F, 1.000F, 1.000F, 0.143F, 0.071F, 0.083F, 0.609F,
            0.500F, 0.714F, 1.000F, 0.022F, 0.071F, 0.143F, 0.214F, 0.043F,
            0.500F, 1.000F, 0.417F, 0.125F, 1.000F, 0.033F, 1.000F, 0.208F,
            0.071F, 0.867F, 1.000F, 0.435F, 1.000F, 0.200F, 0.097F, 0.022F,
            1.000F, 1.000F, 0.417F, 0.955F, 0.188F, 0.025F, 1.000F, 0.773F,
            0.750F, 0.275F, 1.000F, 1.000F, 0.438F, 0.125F, 0.412F, 0.045F
        }
    },
    {
        "Saccharomyces cerevisiae",
        "Sharp and Cowe (1991) Yeast 7:657-678",
        {
            0.00F,
            0.113F, 1.000F, 0.071F, 1.000F, 1.000F, 0.693F, 1.000F, 0.077F,
            0.117F, 0.036F, 0.000F, 0.000F, 1.000F, 0.005F, 0.000F, 1.000F,
            0.006F, 0.047F, 0.245F, 0.137F, 0.003F, 0.009F, 1.000F, 0.002F,
            0.039F, 1.000F, 1.000F, 0.002F, 0.003F, 0.002F, 0.007F, 0.002F,
            0.823F, 0.921F, 0.053F, 0.021F, 1.000F, 1.000F, 1.000F, 0.031F,
            0.003F, 0.012F, 0.135F, 1.000F, 1.000F, 0.006F, 1.000F, 0.003F,
            1.000F, 1.000F, 0.554F, 1.000F, 0.831F, 0.316F, 1.000F, 0.020F,
            0.002F, 0.015F, 1.000F, 0.002F, 0.018F, 0.001F, 0.016F, 0.004F
        }
    }
};

AMINO_PROP_STRUCT amino_prop = {
    { /* hydropathicity values */
        0.00F,
        2.80F, 3.80F, 4.50F, 1.90F, 4.20F, 
        -0.8F, -1.6F, -0.7F, 1.80F, -1.3F,
        1.00F, -3.2F, -3.5F, -3.5F, -3.9F,
        -3.5F, -3.5F, 2.50F, -0.9F, -4.5F, -0.4F
    },
    { /* am i aromatic ?       */
        0,
        1, 0, 0, 0, 0, 
        0, 0, 0, 0, 1,
        0, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0
    }
};

MENU_STRUCT Z_menu = {
    'X',   /*This default is set in proc_commline to CU        */
    false, /*totals                                            */
    true,  /*warnings about sequence data are to be displayed  */

    false, /*fop                                               */
    false, /*cai                                               */
    false, /*cbi                                               */
    false, /*bases                                             */
    false, /*gc3s                                              */
    false, /*gc                                                */
    false, /*enc                                               */
    false, /* silent base                                      */
    false, /* Length silent codons                             */
    false, /* length in codons                                 */
    false, /* hydrophobicity                                   */
    false, /* aromaticity                                      */

    ',', /* default seperator                                */

    0, /* genetic code                                     */
    0, /* type of fop_species                              */
    0, /* type of cai_species                              */

    NULL, /* Null pointer input file                          */
    NULL, /* Null pointer outputfile                          */
    NULL, /* Null pointer tidyout file                        */
    NULL, /* Null codon usage file                            */
    NULL, /* Null pointer fopfile                             */
    NULL, /* Null pointer caifile                             */
    NULL, /* Null pointer cbifile                             */
    NULL, /* Null pointer the logfile name                    */
    NULL, /* error to stderr                                */

    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL
};

REF_STRUCT Z_ref = {
    cu_ref,
    fop_ref,
    cai_ref,

    &amino_acids,
    &amino_prop
};
