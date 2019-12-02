"""
extract any chromosomes from reference genome
_experiment_= "2,3,4"
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
#chr1 to chr Y from reference fasta file
##########################

inputFile= "../../../raw_data/hs37d5.fa"
outputFile = "chr1_22_hs37d5.fa"

records=[]
count=0
for record in SeqIO.parse(inputFile, "fasta"):
    if len(record.id)<=2 and record.id!="MT" and record.id!="X" and record.id!="Y":
        print(record.id)
        print(record.description)
        records.append(record)
        count+=1
SeqIO.write(records, outputFile, "fasta")
print(count)
