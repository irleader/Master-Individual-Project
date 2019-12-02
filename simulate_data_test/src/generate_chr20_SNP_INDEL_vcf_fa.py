from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

inputFile="../../../raw_data/chr20_hs37d5.fa"
outputFile= "../../../raw_data/SNP_INDEL_chr20_hs37d5.fa"
outputFile2="../../../raw_data/SNP_INDEL_chr20_hs37d5.vcf"

foutVCF = open(outputFile2,"w")
foutVCF.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG666\n")

#read reference haplotype in
nucleotide = ['A', 'T', 'C', 'G']

record=SeqIO.read(inputFile,"fasta")



#randomly generate a list of 0.167% (1/600)positions from all positions in the refernece haplotype
# insertion : deletion : SNP = 1:1:8
import random

originseq=list(record.seq.upper())
H1seq=list(record.seq.upper())
H2seq=list(record.seq.upper())
length=len(originseq)
print("chromosome length: ",length)
SNP_position_set=set()
INSERT_position_set=set()
DELETE_position_set=set()
for j in range(int(length/600)):
    position=random.randint(1,length)
    if originseq[position-1]!="N":
        if position % 10 == 3:
            INSERT_position_set.add(position)
            continue
        if position % 10 == 7:
            DELETE_position_set.add(position)
            continue
        else:
            SNP_position_set.add(position)
print("empty set",INSERT_position_set.intersection(DELETE_position_set))
print("empty set", INSERT_position_set.intersection(SNP_position_set))
print("empty set", SNP_position_set.intersection(DELETE_position_set))
INSERT_position_set=sorted(list(INSERT_position_set))
print("insert SIZE", len(INSERT_position_set))
DELETE_position_set=sorted(list(DELETE_position_set))
print("delete SIZE", len(DELETE_position_set))
SNP_position_set=sorted(list(SNP_position_set))
print("SNP SIZE", len(SNP_position_set))

for k in SNP_position_set:
    if originseq[k-1]== "N":
        print("wrong position!")
        continue
    increment=random.randint(1,3)
    index= (nucleotide.index(H2seq[k-1]) + increment) % 4
    H2seq[k-1]=nucleotide[index]
    foutVCF.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (record.id, k , 0, originseq[k-1], H2seq[k-1], '50', 'PASS', '', 'GT:DP:ADALL:AD:GQ:IGT:IPS:PS','1|1:0:0,0:0,0:99:1/1:.:HOMVAR'))
print("size after SNP",len(H2seq))

for k in DELETE_position_set:
    if originseq[k-1]== "N":
        print("wrong position!")
        continue
    H2seq[k-1]="_"
    foutVCF.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (record.id, k-1, 0, originseq[k-2]+originseq[k-1], H2seq[k-2], '50', 'PASS', '', 'GT:DP:ADALL:AD:GQ:IGT:IPS:PS', '1|1:0:0,0:0,0:99:1/1:.:HOMVAR'))
print("size after delete", len(H2seq))

count=0
for k in INSERT_position_set:
    if originseq[k-1]== "N":
        print("wrong position!")
        continue
    H2seq.insert(k+count,nucleotide[random.randint(0,3)])
    foutVCF.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (record.id, k, 0, originseq[k-1], H2seq[k-1+count]+H2seq[k+count], '50', 'PASS', '','GT:DP:ADALL:AD:GQ:IGT:IPS:PS', '1|1:0:0,0:0,0:99:1/1:.:HOMVAR'))
    count=count+1
print("size after insert", len(H2seq))

new_H2seq=[]
for i in H2seq:
    if i!="_":
        new_H2seq.append(i)
print(length,"and",len(INSERT_position_set)-len(DELETE_position_set),"same as",len(new_H2seq))



new_H2seq=SeqRecord(Seq(''.join(e for e in new_H2seq)),id=record.id)


foutVCF.close()
SeqIO.write(new_H2seq, outputFile,'fasta')







