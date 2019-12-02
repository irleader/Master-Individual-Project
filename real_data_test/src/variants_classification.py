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


wgs_keys=dict_wgs_variants.keys()
wgs_pos=set(wgs_keys)
pacbio_keys=dict_pacbio_variants.keys()
pacbio_pos=set(pacbio_keys)


common_pos=wgs_pos.intersection(pacbio_pos)

wgs_unique_pos=list(wgs_pos-common_pos)
pacbio_unique_pos=list(pacbio_pos-common_pos)

common_pos=list(common_pos)

def filter_variants(input,pos):
    output={}
    for i in range(len(pos)):
        output[pos[i]]=input[pos[i]]
    return output


dict_wgs_unique_variants=filter_variants(dict_wgs_variants,wgs_unique_pos)


dict_pacbio_unique_variants=filter_variants(dict_pacbio_variants,pacbio_unique_pos)


dict_wgs_common_variants=filter_variants(dict_wgs_variants,common_pos)


dict_pacbio_common_variants=filter_variants(dict_pacbio_variants,common_pos)


def dict_to_list(input):
    output=[]
    for key in input.keys():
        temp=[]
        temp.append(input[key][0])
        temp.append(key)
        temp=temp+input[key][1:]
        output.append(temp)
    return output

wgs_unique_variants=dict_to_list(dict_wgs_unique_variants)

wgs_common_variants=dict_to_list(dict_wgs_common_variants)

pacbio_unique_variants=dict_to_list(dict_pacbio_unique_variants)

pacbio_common_variants=dict_to_list(dict_pacbio_common_variants)


def write_csv(input1,input2):
    with open(input1,"w") as file:
        writer = csv.writer(file)
        writer.writerow(
            ['chrom', 'pos', 'ref', 'alt', 'info', 't_genotype', 't_decision', 't_decision_sub', 't_quality',
             't_additional_info', 't_variant_type', 't_location_type',
             'q_genotype', 'q_decision', 'q_decision_sub', 'q_quality', 'q_additional_info', 'q_variant_type',
             'q_location_type'])
        writer.writerows(input2)

write_csv("wgs.unique.chr20.csv",wgs_unique_variants)

write_csv("wgs.common.chr20.csv",wgs_common_variants)

write_csv("pacbio.unique.chr20.csv",pacbio_unique_variants)

write_csv("pacbio.common.chr20.csv",pacbio_common_variants)








