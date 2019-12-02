"""
Clean and extract all important info from Hap.py VCF file, remove duplicated variant calls.
Convert Hap.py VCF files into csv files for further processing.

_experiment_= "1"
__author__ = “Jiajia Xu”
__copyright__ = “Copyright 2019, Individual Computing Project”
__email__ = “jiajia.xu@anu.edu.au”
__status__ = “Production”
"""


# This function cleans and extracts important info from VCF file of Hap.py
def read_happy_vcf(filename):
    file = open(filename, "r")
    variants = []
    variant=[]
    count = 0
    pos_list=[]
    len_6_format_row=[]
    for line in file:

        if line.startswith('#'):
            continue
        rows = line.split()
        chrom=int(rows[0])
        pos = int(rows[1])
        ref = rows[3].strip()
        alt = rows[4].strip()
        info=rows[7].strip()
        format=rows[8].strip().split(':')
        truth=rows[9].strip().split(':')
        query=rows[10].strip().split(':')

        assert len(format) == 6 or len(format) == 7

        if len(format)==7:
            t_genotype = truth[format.index('GT')]
            t_decision=truth[format.index('BD')]
            t_decision_sub=truth[format.index('BK')]
            t_quality=truth[format.index('QQ')]
            t_additional_info=truth[format.index('BI')]
            t_variant_type=truth[format.index('BVT')]
            t_location_type=truth[format.index('BLT')]

            q_genotype = query[format.index('GT')]
            q_decision=query[format.index('BD')]
            q_decision_sub=query[format.index('BK')]
            q_quality=query[format.index('QQ')]
            q_additional_info=query[format.index('BI')]
            q_variant_type=query[format.index('BVT')]
            q_location_type=query[format.index('BLT')]

        if len(format)==6:
            len_6_format_row.append(pos)
            t_genotype = truth[format.index('GT')]
            t_decision=truth[format.index('BD')]
            t_decision_sub='.'
            t_quality=truth[format.index('QQ')]
            t_additional_info=truth[format.index('BI')]
            t_variant_type=truth[format.index('BVT')]
            t_location_type=truth[format.index('BLT')]

            q_genotype = query[format.index('GT')]
            q_decision=query[format.index('BD')]
            q_decision_sub='.'
            q_quality=query[format.index('QQ')]
            q_additional_info=query[format.index('BI')]
            q_variant_type=query[format.index('BVT')]
            q_location_type=query[format.index('BLT')]

        variant= [chrom,pos,ref,alt,info,t_genotype,t_decision,t_decision_sub,t_quality,t_additional_info,t_variant_type,t_location_type,
                  q_genotype,q_decision,q_decision_sub,q_quality,q_additional_info,q_variant_type,q_location_type]
        #print(variant)
        variants.append(variant)
        pos_list.append(pos)
        count=count+1
    file.close()
    return pos_list,count,variants,len_6_format_row




# Use the above function to read Hap.py output VCF (from both Illumina and Pacbio reads)
p1,count1,wgs_variants,len_6_format_row_wgs=read_happy_vcf("../../raw_Data/wgs.chr20.happy.output.vcf")
p2,count2,pacbio_variants,len_6_format_row_pacbio=read_happy_vcf("../../raw_Data/pacbio.chr20.happy.output.vcf")

print('wgs original position size: ',len(p1))
print('pacbio original position size: ',len(p2))




# Check the duplicated variant calls and remove them from Illumina, PacBio
duplicate1=p1.copy()
duplicate2=p2.copy()

for x in set(duplicate1):
    duplicate1.remove(x)
print('wgs duplicate position size: ',len(duplicate1))
for x in set(duplicate2):
    duplicate2.remove(x)
print('pacbio duplicate position size: ',len(duplicate2))


duplicate=set(duplicate1).union(set(duplicate2))
print('total duplicate position size:',len(duplicate))

wgs_indices=set(p1).difference(duplicate)
print('wgs position size after removing mutual position: ', len(wgs_indices))
pacbio_indices=set(p2).difference(duplicate)
print('pacbio position size after removing mutual position: ', len(pacbio_indices))

wgs_single_indices=set(p1).difference(set(duplicate1))
print('wgs single position size: ',len(wgs_single_indices))
pacbio_single_indices=set(p2).difference(set(duplicate2))
print('pacbio single position size: ',len(pacbio_single_indices))





# store Illumina and PacBio variant calls into dictionary
modified_wgs_variants=[]
for i in range(len(wgs_variants)):
    if wgs_variants[i][1] in wgs_indices:
        modified_wgs_variants.append(wgs_variants[i])

modified_pacbio_variants=[]
for i in range(len(pacbio_variants)):
    if pacbio_variants[i][1] in pacbio_indices:
        modified_pacbio_variants.append(pacbio_variants[i])


wgs_single_variants=[]
for i in range(len(wgs_variants)):
    if wgs_variants[i][1] in wgs_single_indices:
        wgs_single_variants.append(wgs_variants[i])

pacbio_single_variants=[]
for i in range(len(pacbio_variants)):
    if pacbio_variants[i][1] in pacbio_single_indices:
        pacbio_single_variants.append(pacbio_variants[i])


wgs_duplicate_variants=[]
for i in range(len(wgs_variants)):
    if wgs_variants[i][1] in duplicate1:
        wgs_duplicate_variants.append(wgs_variants[i])

pacbio_duplicate_variants=[]
for i in range(len(pacbio_variants)):
    if pacbio_variants[i][1] in duplicate2:
        pacbio_duplicate_variants.append(pacbio_variants[i])




# save variant calls into CSV file
import csv

with open("wgs.chr20.csv","w") as file:
    writer=csv.writer(file)
    writer.writerow(['chrom','pos','ref','alt','info','t_genotype','t_decision','t_decision_sub','t_quality','t_additional_info','t_variant_type','t_location_type',
                  'q_genotype','q_decision','q_decision_sub','q_quality','q_additional_info','q_variant_type','q_location_type'])
    writer.writerows(modified_wgs_variants)

with open("pacbio.chr20.csv","w") as file:
    writer=csv.writer(file)
    writer.writerow(['chrom','pos','ref','alt','info','t_genotype','t_decision','t_decision_sub','t_quality','t_additional_info','t_variant_type','t_location_type',
                  'q_genotype','q_decision','q_decision_sub','q_quality','q_additional_info','q_variant_type','q_location_type'])
    writer.writerows(modified_pacbio_variants)





#clean truth variant calls and store them in CSV
with open("truth.chr20.csv","r") as file:
    reader=csv.reader(file)
    truth_variants=list(reader)
    truth_variants.pop(0)

total_count=len(truth_variants)
print('truth original position size',total_count)

import numpy as np
truth_variants=np.array(truth_variants)
truth_position=list(map(int,list(truth_variants[:,1])))

truth_indices=list(set(truth_position).difference(duplicate))
print('truth position size after removing mutual position',len(truth_indices))

truth_variants=truth_variants.tolist()

modified_truth_variants=[]
for i in range(len(truth_variants)):
    if int(truth_variants[i][1]) in truth_indices:
        modified_truth_variants.append(truth_variants[i])

with open("cleaned_truth.chr20.csv","w") as file:
    writer=csv.writer(file)
    writer.writerow(['chrom','pos','ref','alt','info','qual', 'filter', 'info', 'genotype', 'depth','depth_all_allele',
                   'depth_each_allele', 'quality', 'origin_genotype', 'phase_IGT','phase_GT'])
    writer.writerows(modified_truth_variants)



