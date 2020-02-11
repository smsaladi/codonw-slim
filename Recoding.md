Data Recoding
=============

CodonW converts sequences into a numerical format for computation
under the hood. The python interface uses pandas Series.


Nucleotides are recoded as T/U=1, C=2, A=3, G=4. The 20 standard
amino acids and the termination codons are recoded as integer values
in the range 1 to 21, note that stop codons are assigned the amino
acid value 11. The decision about whether a codon is
synonymous, or how many members are in a particular amino acid
synonymous family are dependent on the genetic code chosen.


Numerical values used to recode amino acids.
Code    AA3    AA1     Code    AA3   AA1
1       Phe    F       2       Leu    L
3       Ile    I       4       Met    M
5       Val    V       6       Ser    S
7       Pro    P       8       Thr    T
9       Ala    A       10      Tyr    Y
11      STOP   *       12      His    H
13      Gln    Q       14      Asn    N
15      Lys    K       16      Asp    D
17      Glu    E       18      Cys    C
19      Trp    W       20      Arg    R
21      Gly    G




Each codon is recoded into an integer value in the range 1 to 64,
see Table 1. The formulae used to recode the codons is:
        
    code = ((p1-1)*16)+P2+((p3-1)*4)    1 <= code <= 64

Where each of the three codon positions is represented by P1, P2 and
P3. Using this recoding convention, the codon ATG has the value 45.
 
    code = ((3-1)*16)+1+((4-1)*4) = 45

Unrecognised or non-translatable bases, codons or amino acids are
all assigned the value `0`.


Code  Codon   AA   Code   Codon    AA    Code   Codon   AA      Code  Codon  AA
1      UUU    Phe   2      UCU     Ser   3       UAU    Tyr     4     UGU    Cys
5      UUC          6      UCC           7       UAC            8     UGC
9      UUA    Leu   10     UCA           11      UAA    STOP    12    UGA    STOP
13     UUG          14     UCG           15      UAG            16    UGG    Trp
17     CUU          18     CCU     Pro   19      CAU    His     20    CGU    Arg
21     CUC          22     CCC           23      CAC            24    CGC
25     CUA          26     CCA           27      CAA    Gln     28    CGA
33     AUU     Ile  34     ACU     Thr   35      AAU    Asn     36    AGU    Ser
37     AUC          38     ACC           39      AAC            40    AGC
29     CUG          30     CCG           31      CAG            32    CGG
41     AUA          42     ACA           43      AAA    Lys     44    AGA    Arg
45     AUG     Met  46     ACG           47      AAG            48    AGG
49     GUU     Val  50     GCU     Ala   51      GAU    Asp     52    GGU    Gly
53     GUC          54     GCC           55      GAC            56    GGC
57     GUA          58     GCA           59      GAA    Glu     60    GGA
61     GUG          62     GCG           63      GAG            64    GGG
