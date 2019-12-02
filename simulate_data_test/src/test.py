from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#
# inputFile="../../../raw_data/2H_chr20_hs37d5.fa"
# outputFile="../../../raw_data/artificial_chr20_hs37d5.fa"
#
#
# records=list(SeqIO.parse(inputFile,"fasta"))
# print(records[0].id)
# print(records[1].id)
# record=SeqRecord(records[1].seq,id="20",description="dna:chromosome chromosome:GRCh37:20:1:63025520:1")
# SeqIO.write(record,outputFile,"fasta")
import random
# nucleotide = ['A', 'T', 'C', 'G']
#
# o=["A","C","G","T"]
# a=["A","C","G","T"]
# p=[1,2,3]
# count = 0
# for k in p:
#     a.insert(k + count, "O")
#     print(k,o[k-1]," ",a[k-1+count]+a[k+count])
#     count = count + 1
# print(a)

# p1=[3]
# for k in p1:
#     a[k-1]="_"
#     print(a)
#     print(k-1,o[k-2]+o[k-1],a[k-2])

# inputfile="../../../raw_data/chr20_hs37d5.fa"
# sequence=list(SeqIO.read(inputfile,"fasta").seq)
# print(sequence[64584],sequence[64585],sequence[64586])
import csv


file2=open("pacbio.simulate.output.csv",'r')
reader2=csv.reader(file2)
P=list(reader2)
#P.pop(0)