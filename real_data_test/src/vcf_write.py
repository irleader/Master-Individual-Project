# import vcf
# vcf_reader=vcf.Reader(open("../wgs.chr20.happy.output.vcf",'r'))
# for record in vcf_reader:
#     print(record)

import csv

with open('wgs.chr20.csv','r') as file:
    reader=csv.reader(file)
    modified_wgs_variants=list(reader)
    modified_wgs_variants.pop(0)

with open('pacbio.chr20.csv','r') as file:
    reader=csv.reader(file)
    modified_pacbio_variants=list(reader)
    modified_pacbio_variants.pop(0)

dict_wgs_variants={}
for i in range(len(modified_wgs_variants)):
    dict_wgs_variants[modified_wgs_variants[i][1]]=modified_wgs_variants[i][0:1]+modified_wgs_variants[i][2:]

dict_pacbio_variants={}
for i in range(len(modified_pacbio_variants)):
    dict_pacbio_variants[modified_pacbio_variants[i][1]] = modified_pacbio_variants[i][0:1] + modified_pacbio_variants[i][2:]




