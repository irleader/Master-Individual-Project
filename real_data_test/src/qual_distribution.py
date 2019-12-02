#Phred quality scores Q, base-calling error probabilities P
# Q=-10 * log10 (P)

#QUAL: output vcf
#Phred-scaled quality score for the assertion made in ALT. i.e. −10log10 prob(call in ALT is wrong).
# If ALT is ‘.’ (no variant) then this is −10log10 prob(variant), and if ALT is not ‘.’ this is −10log10 prob(no variant).
# If unknown, the missing value should be specified. (Numeric)


# output vcf: (column 9)
#GQ: conditional genotype quality

# truth vcf: (column 12) Variant quality for ROC creation
#GQ: Net Genotype quality across all datasets, calculated from GQ scores of callsets supporting the consensus GT,
# using only one callset from each dataset


import csv
import numpy as np

with open("wgs.chr20.output.csv","r") as file:
    wgs_reader=csv.reader(file)
    wgs_variants=list(wgs_reader)
    wgs_variants.pop(0)
    wgs_variants_array = np.array(wgs_variants)

with open("pacbio.chr20.output.csv","r") as file:
    pacbio_reader=csv.reader(file)
    pacbio_variants=list(pacbio_reader)
    pacbio_variants.pop(0)
    pacbio_variants_array = np.array(pacbio_variants)

with open("wgs.chr20.csv","r") as file:
    wgs_hap_reader=csv.reader(file)
    wgs_hap_variants=list(wgs_hap_reader)
    wgs_hap_variants.pop(0)
    wgs_hap_variants_array = np.array(wgs_hap_variants)

with open("pacbio.chr20.csv","r") as file:
    pacbio_hap_reader=csv.reader(file)
    pacbio_hap_variants=list(pacbio_hap_reader)
    pacbio_hap_variants.pop(0)
    pacbio_hap_variants_array = np.array(pacbio_hap_variants)

with open("cleaned_truth.chr20.csv","r") as file:
    truth_reader=csv.reader(file)
    truth_variants=list(truth_reader)
    truth_variants.pop(0)

truth_variants_array=np.array(truth_variants)

#output qual and GQ
wgs_qual=wgs_variants_array[:,[1,5]]
pacbio_qual=pacbio_variants_array[:,[1,5]]
truth_qual=truth_variants_array[:,[1,5]]

wgs_GQ=wgs_variants_array[:,[1,9]]
pacbio_GQ=pacbio_variants_array[:,[1,9]]
truth_GQ=truth_variants_array[:,[1,12]]

#hap query QQ,same as qual
wgs_hap_qQQ=wgs_hap_variants_array[:,[1,15]]
pacbio_hap_qQQ=pacbio_hap_variants_array[:,[1,15]]

#check if query QQ is same as qual
wgs_index=list(wgs_qual[:,0])
wgs_hap_index=list(wgs_hap_qQQ[:,0])
wgs_hap_qQQ_value=list(map(float,list(wgs_hap_qQQ[:,1])))
wgs_qual_value=list(map(float,list(wgs_qual[:,1])))

index=[]
index2=[]
for i in range(len(wgs_hap_index)):
    if wgs_hap_index[i] in wgs_index:
        index.append(wgs_index.index(wgs_hap_index[i]))
        index2.append(i)

new_wgs_qual_value=[]
for i in range(len(index)):
    new_wgs_qual_value.append(wgs_qual_value[index[i]])

new_wgs_hap_qQQ_value=[]
for i in range(len(index2)):
    new_wgs_hap_qQQ_value.append(wgs_hap_qQQ_value[index2[i]])


# with open("wgs.chr20.c3.csv","r") as file:
#     wgs_hap_reader=csv.reader(file)
#     wgs_hap_c3_variants=list(wgs_hap_reader)
#     wgs_hap_c3_variants.pop(0)
#     wgs_hap_c3_variants=np.array(wgs_hap_c3_variants)
#
# with open("pacbio.chr20.c3.csv","r") as file:
#     pacbio_hap_reader=csv.reader(file)
#     pacbio_hap_c3_variants=list(pacbio_hap_reader)
#     pacbio_hap_c3_variants.pop(0)
#     pacbio_hap_c3_variants=np.array(pacbio_hap_c3_variants)
#
# with open("wgs.chr20.c6.csv","r") as file:
#     wgs_hap_reader=csv.reader(file)
#     wgs_hap_c6_variants=list(wgs_hap_reader)
#     wgs_hap_c6_variants.pop(0)
#     wgs_hap_c6_variants=np.array(wgs_hap_c6_variants)
#
# with open("pacbio.chr20.c6.csv","r") as file:
#     pacbio_hap_reader=csv.reader(file)
#     pacbio_hap_c6_variants=list(pacbio_hap_reader)
#     pacbio_hap_c6_variants.pop(0)
#     pacbio_hap_c6_variants=np.array(pacbio_hap_c6_variants)
#
# wgs_hap_c3_qQQ=list(map(float,list(wgs_hap_c3_variants[:,15])))
# wgs_hap_c6_qQQ=list(map(float,list(wgs_hap_c6_variants[:,15])))
# pacbio_hap_c3_qQQ=list(map(float,list(pacbio_hap_c3_variants[:,15])))
# pacbio_hap_c6_qQQ=list(map(float,list(pacbio_hap_c6_variants[:,15])))
#
#
# import matplotlib.pyplot as plt
#
# plt.hist(wgs_hap_c3_qQQ,bins='auto',alpha=0.7)
# plt.hist(wgs_hap_c6_qQQ,bins='auto',alpha=0.7)
# plt.show()
#
# plt.hist(pacbio_hap_c3_qQQ,bins='auto',alpha=0.7)
# plt.hist(pacbio_hap_c6_qQQ,bins='auto',alpha=0.7)
# plt.show()
