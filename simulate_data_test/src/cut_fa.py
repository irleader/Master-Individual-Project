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