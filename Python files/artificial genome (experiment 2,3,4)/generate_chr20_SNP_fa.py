"""
Generate artificial genome with SNP only for chr20
_experiment_= "2"
__author__ = “Jiajia Xu”
__copyright__ = “Copyright 2019, Individual Computing Project”
__email__ = “jiajia.xu@anu.edu.au”
__status__ = “Production”
"""

import collections
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import os
import sys

##########################
# generate another haplotype from reference fasata file
##########################



inputFile="../../../raw_data/chr20_hs37d5.fa"
outputFile= "../artificial_chr20_hs37d5.fa"

#read reference haplotype in
nucleotide = ['A', 'T', 'C', 'G']
record=SeqIO.read(open(inputFile),"fasta")
print(record.id)
H2list=list(record.seq.upper())
length=len(H2list)

#randomly generate a list of 0.1% positions from all positions in the refernece haplotype
import random
position_set=set()
for i in range(int(length/1000)):
    position = random.randint(1, length)
    position_set.add(position)
position_list = sorted(list(position_set))

#randomly change the nucleotide type from the positions list to generate a different haplotype
count = 0
for i in position_list:
    if H2list[i] == 'N':
        continue
    increment = random.randint(1, 3)
    index = (nucleotide.index(H2list[i]) + increment) % 4
    H2list[i] = nucleotide[index]
    count += 1

#store two haplotypes in one fasta file
rec2 = SeqRecord(Seq(''.join(e for e in H2list)), id=record.id,description=record.description)
SeqIO.write(rec2, outputFile, 'fasta')




