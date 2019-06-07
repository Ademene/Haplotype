#!/usr/bin/env python
# -*- coding: utf-8 -*-

### usage in prompt : $ python Fasta_to_GenomeFile.py Genome.fasta


import sys
filein= sys.argv[1]

### The fasta must not be wrapped to 80 characters, and have the whole sequence on the same line for each ">".

INPUT=open(filein,'r')
OUT=open(filein.replace('.fasta','.genome'),'w')

### Print the header
OUT.write("Contig\tLength\n")

for lines in INPUT:
    if lines.startswith('>'):
        Scaf=(lines[1:].strip())
        OUT.write(Scaf + '\t')
    else:
        OUT.write(str(len(lines)) + '\n')
        
INPUT.close()
OUT.close()
