# import vcf
# vcf_reader=vcf.Reader(open("../wgs.chr20.happy.output.vcf",'r'))
# for record in vcf_reader:
#     print(record)
#
# import allel
# callset=allel.read_vcf('../pacbio.chr20.happy.output.vcf')
# print(callset.keys())

# import vcf
# import csv
#
# header_reader = vcf.Reader(filename='../header.vcf.gz')
# header_writer = vcf.Writer(open('../../test.vcf', 'w'), header_reader)
# for record in header_reader:
#     header_writer.write_record(record)


# def read_happy_vcf(filename):
#     file = open(filename, "r")
#     variants = []
#     variant=[]
#     count = 0
#     pos_list=[]
#     len_6_format_row=[]
#     for line in file:
#
#         if line.startswith('#'):
#             continue
#         print(line)
#         print(type(line))
#         rows = line.split()
#         chrom=int(rows[0])
#         pos = int(rows[1])
#         ref = rows[3].strip()
#         alt = rows[4].strip()
#         info=rows[7].strip()
#         format=rows[8].strip().split(':')
#         truth=rows[9].strip().split(':')
#         query=rows[10].strip().split(':')
#
#
#
#         assert len(format) == 6 or len(format) == 7
#
#         if len(format)==7:
#             t_genotype = truth[format.index('GT')]
#             t_decision=truth[format.index('BD')]
#             t_decision_sub=truth[format.index('BK')]
#             t_quality=truth[format.index('QQ')]
#             t_additional_info=truth[format.index('BI')]
#             t_variant_type=truth[format.index('BVT')]
#             t_location_type=truth[format.index('BLT')]
#
#             q_genotype = query[format.index('GT')]
#             q_decision=query[format.index('BD')]
#             q_decision_sub=query[format.index('BK')]
#             q_quality=query[format.index('QQ')]
#             q_additional_info=query[format.index('BI')]
#             q_variant_type=query[format.index('BVT')]
#             q_location_type=query[format.index('BLT')]
#
#         if len(format)==6:
#             len_6_format_row.append(pos)
#             t_genotype = truth[format.index('GT')]
#             t_decision=truth[format.index('BD')]
#             t_decision_sub='.'
#             t_quality=truth[format.index('QQ')]
#             t_additional_info=truth[format.index('BI')]
#             t_variant_type=truth[format.index('BVT')]
#             t_location_type=truth[format.index('BLT')]
#
#             q_genotype = query[format.index('GT')]
#             q_decision=query[format.index('BD')]
#             q_decision_sub='.'
#             q_quality=query[format.index('QQ')]
#             q_additional_info=query[format.index('BI')]
#             q_variant_type=query[format.index('BVT')]
#             q_location_type=query[format.index('BLT')]
#             continue
#
#
#
#         variant= [chrom,pos,ref,alt,info,t_genotype,t_decision,t_decision_sub,t_quality,t_additional_info,t_variant_type,t_location_type,
#                   q_genotype,q_decision,q_decision_sub,q_quality,q_additional_info,q_variant_type,q_location_type]
#
#         #print(variant)
#         variants.append(variant)
#         pos_list.append(pos)
#         count=count+1
#
#     file.close()
#     return pos_list,count,variants,len_6_format_row
#
# p1,count1,wgs_variants,len_6_format_row_wgs=read_happy_vcf("../wgs.chr20.happy.output.vcf")
# p2,count2,pacbio_variants,len_6_format_row_pacbio=read_happy_vcf("../pacbio.chr20.happy.output.vcf")

#clean truth csv
# import csv
#
# with open("truth.chr20.csv","r") as file:
#     reader=csv.reader(file)
#     truth_variants=list(reader)
#     truth_variants.pop(0)
#
# total_count=len(truth_variants)
# print('truth original position size',total_count)
#
# import numpy as np
# truth_variants=np.array(truth_variants)
# truth_position=list(map(int,list(truth_variants[:,1])))

# import numpy as np
# SD=np.loadtxt("../../raw_data/MosaicSDs_Human_hg19.txt",delimiter="/t",dtype=str)

file=open("../../raw_data/MosaicSDs_Human_hg19.txt", 'r')

SD=[]
for line in file:
    tmp=[]
    for element in line.split():
        tmp.append(element)
    SD.append(tmp)

chr20_SD=[]

for i in range(len(SD)):
    if SD[i][0]=='chr20':
        chr20_SD.append(SD[i][1:3])

for i in range(len(chr20_SD)):
    chr20_SD[i]=list(map(int,chr20_SD[i]))


