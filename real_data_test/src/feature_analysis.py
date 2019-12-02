import csv

with open("pacbio.chr20.csv","r") as file:
    reader=csv.reader(file)
    variants=list(reader)
    variants.pop(0)

total_count=len(variants)

# info
confident_region_count=0
for i in range(len(variants)):
    #print(type(wgs_unique_variants[i][4]))
    if "Regions=CONF,TS_contained" in variants[i][4]:
        confident_region_count+=1
#print(confident_region_count)


import numpy as np
variants_array=np.array(variants)

info=variants_array[:,4]

t_location=variants_array[:,11]
#t_genotype=variants_array[:,5]
t_deicision=variants_array[:,6]

q_location=variants_array[:,18]
#q_genotype=variants_array[:,12]
q_deicision=variants_array[:,13]


# category 1: called by truth, not called by query, FN by truth, confident area
c1=[]
for i in range(len(variants)):
    if t_location[i]!='nocall' and q_location[i]=='nocall' and t_deicision[i]=='FN' and q_deicision[i]=='.' \
            and ('Regions=CONF,TS_contained' or 'Regions=CONF,TS_boundary' in info[i]):
        c1.append(variants[i])

# category 2: called by truth, not called by query, UNK by truth, unconfident area
c2=[]
for i in range(len(variants)):
    if t_location[i]!='nocall' and q_location[i]=='nocall' and t_deicision[i]=='UNK' and q_deicision[i]=='UNK' \
            and (('Regions=CONF,TS_contained' and "Regions=CONF,TS_boundary") not in info[i]):
        c2.append(variants[i])

# category 3: called by truth, called by query, TP by truth, confident area
c3=[]
for i in range(len(variants)):
    if t_location[i]!='nocall' and q_location[i]!='nocall' and t_deicision[i]=='TP' and q_deicision[i]=='TP' \
            and ('Regions=CONF,TS_boundary' or 'Regions=CONF,TS_contained' in info[i]):
        # if t_location[i]!=q_location[i]:
        #     print(t_location[i],q_location[i])
        c3.append(variants[i])

# category 4: called by truth, called by query, UNK by truth, unconfident area
c4=[]
for i in range(len(variants)):
    if t_location[i]!='nocall' and q_location[i]!='nocall' and t_deicision[i]=='UNK' and q_deicision[i]=='UNK' \
            and (("Regions=CONF,TS_contained" and "Regions=CONF,TS_boundary") not in info[i]):
        # if t_location[i]!=q_location[i]:
        #     print(t_location[i],q_location[i])
        c4.append(variants[i])

# category 5: not called by truth, called by query, FP by query, confident area
c5=[]
for i in range(len(variants)):
    if t_location[i]=='nocall' and q_location[i]!='nocall' and t_deicision[i]=='.' and q_deicision[i]=='FP' \
            and ("Regions=CONF,TS_contained" or "Regions=CONF,TS_boundary" in info[i]):
        c5.append(variants[i])

# category 6: not called by truth, called by query, UNK by truth, unconfident area
c6=[]
for i in range(len(variants)):
    if t_location[i]=='nocall' and q_location[i]!='nocall' and t_deicision[i]=='UNK' and q_deicision[i]=='UNK' \
            and (("Regions=CONF,TS_contained" and "Regions=CONF,TS_boundary") not in info[i]):
        c6.append(variants[i])

# category 7: called by truth, called by query, FN by truth, FP by query, confident area
c7=[]
for i in range(len(variants)):
    if t_location[i]!='nocall' and q_location[i]!='nocall' and t_deicision[i]=='FN' and q_deicision[i]=='FP' \
            and ("Regions=CONF,TS_contained" or "Regions=CONF,TS_boundary" in info[i]):
        # if t_location[i]!=q_location[i]:
        #     print(t_location[i],q_location[i])
        c7.append(variants[i])

# category 8: not called by truth, called by query, TP by query, confident area
c8=[]
for i in range(len(variants)):
    if t_location[i]=='nocall' and q_location[i]!='nocall' and t_deicision[i]=='.' and q_deicision[i]=='TP' \
            and ("Regions=CONF,TS_contained" or "Regions=CONF,TS_boundary" in info[i]):
        # if t_location[i]!=q_location[i]:
        #     print(t_location[i],q_location[i])
        c8.append(variants[i])

# category 9: called by truth, not called by query, TP by truth, confident area
c9=[]
for i in range(len(variants)):
    if t_location[i]!='nocall' and q_location[i]=='nocall' and t_deicision[i]=='TP' and q_deicision[i]=='.' \
            and ("Regions=CONF,TS_contained" or "Regions=CONF,TS_boundary" in info[i]):
        c9.append(variants[i])


# import csv
#
# with open("wgs.chr20.c3.csv","w") as file:
#     writer=csv.writer(file)
#     writer.writerow(['chrom','pos','ref','alt','info','t_genotype','t_decision','t_decision_sub','t_quality','t_additional_info','t_variant_type','t_location_type',
#                   'q_genotype','q_decision','q_decision_sub','q_quality','q_additional_info','q_variant_type','q_location_type'])
#     writer.writerows(c3)
#
# with open("wgs.chr20.c6.csv","w") as file:
#     writer=csv.writer(file)
#     writer.writerow(['chrom','pos','ref','alt','info','t_genotype','t_decision','t_decision_sub','t_quality','t_additional_info','t_variant_type','t_location_type',
#                   'q_genotype','q_decision','q_decision_sub','q_quality','q_additional_info','q_variant_type','q_location_type'])
#     writer.writerows(c6)












