import collections
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import os
import sys

##########################
#chr20 from reference fasta file
##########################

refFile="../../../raw_data/references/SNP_INDEL_chr1_22_hs37d5.fa"
outFile = "SNP_INDEL_chr201.fa"

for record in SeqIO.parse(refFile, "fasta"):
    if record.id == "20":
        SeqIO.write(record, outFile, "fasta")
        
    
    
    