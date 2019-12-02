"""
Gnerate ground truth VCF for chr20 with reference and artificial genome
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

inputFile1="../../../raw_data/chr20_hs37d5.fa"
inputFile2= "../../../raw_data/artifical_hs37d5.fa"
outputFile="../../../raw_data/artificial_chr20_hs37d5.vcf"

foutVCF = open(outputFile,"w")
foutVCF.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG666\n")

records1=list(SeqIO.parse(inputFile1,"fasta"))
records2=list(SeqIO.parse(inputFile2,"fasta"))

record1=list(records1[19].seq)
print(records1[19].id)
record2=list(records2[19].seq)
print(records2[19].id)

length=len(record1)

x_list=[]
for i in range(length):
    if record1[i]!=record2[i]:
        x_list.append(i)

print(len(x_list))
count = 1
for i in x_list:

    foutVCF.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (20, i+1, count, record1[i], record2[i], '50', 'PASS', '', 'GT:DP:ADALL:AD:GQ:IGT:IPS:PS', '1|1:0:0,0:0,0:99:1/1:.:HOMVAR') )
    count += 1

foutVCF.close()


