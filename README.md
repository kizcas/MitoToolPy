MitoToolPy

Copyright © 2014 MitoTool & DomeTree Team 

Author: Long Fan (mitotool@gmail.com)


===================================
##Brief introduction

MitoToolPy written in python and running on various operation systems (e.g. Windows, Linux, Mac), is designed for data analysis of mitochondrial DNA of eight domestic animals, including cattle, chicken, dog, horse, goat, pig, sheep and yak. It utilizes dometree (www.dometree.org) as the reference, and its analysis for mtDNA sequences include (1) outputting variants scored relative to the reference sequences; (2) classifying haplogroup status (haplogrouping); and (3) checking potential errors. MitoToolPy can also handle partial mtDNA genome (e.g. control region) sequences and can serve as a convenient "barcoding tool" for strains or breeds of domestic animals.

MitoToolPy is licensed under GPLv3, and consists of two python scripts:

(1) MitoToolPy-Seq

    mitotoolpy-seq.py - classify haplogroup according to mtDNA sequence.

(2) MitoToolPy-Var

    mitotoolpy-var.py - classify haplogroup according to the variants of mtDNA sequence.
    mitotoolpy-var.py - scan possible errors with the help of speculated haplogroups.


===================================
##Installation

MitoToolPy requires the preinstallation of Python and Biopython. Meanwhile, ClulstalW 2 will be downloaded together for pairwise alignment between reference sequence and query sequence. Therefore, after the download, please make sure that ClustalW (in 'bin' directory) has right to execute on your system. Supposed your system is not mentioned below, please download one of below versions and replace original ClustalW (in 'bin' directory) with compiled ClustalW of your system.


===================================
##Manual of MitoToolPy

(1) MitoToolPy-Seq (i.e., mitotoolpy-seq.py)

This script is used for classifying haplogroup according to input sequence(s), outputting variants scored relative to the reference sequences, listing potential missing variants and private variants.

The brief version of the following usage instructions can be obtained by entering `python mitotoolpy-seq.py -h` at

Usage: `python mitotoolpy-seq.py -s [species] -r [region] -i [input] -o [output]`

where

**[species]** = an animal selected from 'cattle', 'chicken', 'dog', 'horse', 'goat', 'pig', 'sheep' and 'yak'.

**[region]** = a region selected from 'dloop', 'nondloop', 'whole' and '[begin_position]:[end_position]'; Default: 'whole'.

**[input]** = input fasta file.

**[output]** = output file or 'stdout'; Default: 'stdout' (i.e., the interface of command line prompt).

Usage: `python mitotoolpy-seq.py -v` , which displays current version number.

Usage: `python mitotoolpy-seq.py -h` , which displays help document.

Notes:

(1) Input file would be a single fasta file, which could store single record or multiple records.

(2) Output file would be a tab-delimited text file with five columns:

*Column 1)  Sample name

*Column 2)  Recommended Haplogroup

*Column 3)  Missing variants of recommended haplogroup

*Column 4)  Private variants of recommended haplogroup

*Column 5)  Variants of selected region                   


(2) MitoToolPy-Var (i.e., mitotoolpy-var.py)

This script is used for classifying haplogroup according to input variants, listing potential missing variants and private variants. Given candidate haplogroup(s) is/are provided by users, it could also detect missing variants of claimed haplogroup.

The brief version of the following usage instructions can be obtained by entering "python mitotoolpy-var.py -h" at

Usage: `python mitotoolpy-var.py -s [species] -r [region] -i [input] -o [output]`

where
**[species]** = an animal selected from 'cattle', 'chicken', 'dog', 'horse', 'goat', 'pig', 'sheep' and 'yak'.

**[region]** = a region selected from 'dloop', 'nondloop', 'whole' and '[begin_position]:[end_position]'; Default: 'whole'.

**[input]** = input text file, see following notes for details.

**[output]** = output file or 'stdout'; Default: 'stdout' (i.e., the interface of command line prompt).

`python mitotoolpy-var.py -v` , which displays current version number.

`python mitotoolpy-var.py -h` , which displays help document.

Notes:

(1) Input file should be a tab-delimited text file with two or three columns:

*Column 1)  Sample name

*Column 2)  Variants

*Column 3 (optional))  Candidate haplogroup(s) claimed by user. Multiple haplogroups seperated by comma, colon or semicolon could be given together.

(2) Output file would be a tab-delimited text file with five columns:

*Column 1)  Sample name

*Column 2)  Recommended Haplogroup

*Column 3)  Missing variants of recommended haplogroup

*Column 4)  Private variants of recommended haplogroup

*Column 5)  Missing variants of haplogroup(s) claimed by user

*Column 6)  Variants of selected region                   


===================================
##History

Oct 26, 2014: Version 1.0 released

