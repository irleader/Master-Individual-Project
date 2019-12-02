from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

inputFile="../../raw_data/2H_chr20_hs37d5.fa"
outputFile="../copy_chr20_hs37d5.vcf"

foutVCF = open(outputFile,"w")
foutVCF.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG666\n")

records=[]
for record in SeqIO.parse(inputFile, "fasta"):
    records.append(record)

H1list=list(records[0].seq)
H2list=list(records[1].seq)

length=len(H1list)

x_list=[]
for i in range(length):
    if H1list[i]!=H2list[i]:
        x_list.append(i)

print(x_list)
count = 1
for i in x_list:

    foutVCF.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (20, i+1, count, H1list[i], H2list[i], '50', 'PASS', '', 'GT:DP:ADALL:AD:GQ:IGT:IPS:PS', '0|1:0:0,0:0,0:99:0/1:.:HOMVAR') )
    count += 1

foutVCF.close()