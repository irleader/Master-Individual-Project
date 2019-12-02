"""
Read CSV variant calls from Illumina, PacBio, and Truth
Divide all three files by intersections into 7 regions
Draw distribution graphs for variant information
Gengerate CSV for region C and E SNP
_experiment_= "3"
__author__ = “Jiajia Xu”
__copyright__ = “Copyright 2019, Individual Computing Project”
__email__ = “jiajia.xu@anu.edu.au”
__status__ = “Production”
"""

# read CSV variant calls files
import csv
import numpy as np

with open("wgs.simulate.csv","r") as file:
    reader=csv.reader(file)
    I_and_T=list(reader)
    I_and_T.pop(0)
I_and_T_array=np.array(I_and_T)
I_and_T_pos= I_and_T_array[:, 1]
I_and_T_pos=list(map(int, I_and_T_pos.tolist()))

with open("pacbio.simulate.csv","r") as file:
    reader=csv.reader(file)
    P_and_T=list(reader)
    P_and_T.pop(0)
P_and_T_array=np.array(P_and_T)
P_and_T_pos= P_and_T_array[:, 1]
P_and_T_pos=list(map(int, P_and_T_pos.tolist()))

with open("cleaned_truth.simulate.csv","r") as file:
    reader=csv.reader(file)
    T=list(reader)
    T.pop(0)
T_array=np.array(T)
T_pos=T_array[:,1]
T_pos=list(map(int,T_pos.tolist()))

total_pos= list(set(P_and_T_pos).union(set(I_and_T_pos)).union(set(T_pos)))


print('total size: ',len(total_pos))
print('truth set size: ', len(T_pos))
P_pos=[]
for i in range(len(P_and_T_pos)):
    #if P_and_T[i][18]!='nocall':
        P_pos.append(int(P_and_T[i][1]))
print('pacbio set size: ',len(P_pos))
I_pos=[]
for i in range(len(I_and_T_pos)):
    #if I_and_T[i][18]!='nocall':
        I_pos.append(int(I_and_T[i][1]))
print('illumina set size: ',len(I_pos))

I_and_T_dict={}
for i in range(len(I_and_T_pos)):
    I_and_T_dict[I_and_T_pos[i]]=I_and_T[i]
P_and_T_dict={}
for i in range(len(P_and_T_pos)):
    P_and_T_dict[P_and_T_pos[i]]=P_and_T[i]

# A: (Union of FN from illumina and pacbio),  FN by truth decsision, nocall or wrong call by query
A_pos=set()
for i in I_and_T_dict.keys():
    if I_and_T_dict[i][6]=="FN":
        A_pos.add(i)
for i in P_and_T_dict.keys():
    if P_and_T_dict[i][6]=="FN":
        A_pos.add(i)
A_pos=list(A_pos)
print('A size',len(A_pos))

A=[]
for i in range(len(A_pos)):
    A_element= I_and_T_dict[A_pos[i]] + P_and_T_dict[A_pos[i]]
    A.append(A_element)


# B and F: (FP by illumina), FP by query decision
B_and_F_pos=set()
for i in I_and_T_dict.keys():
    if I_and_T_dict[i][13]=="FP":
        B_and_F_pos.add(i)



# C and F: (FP by pacbio), FP by query decision
C_and_F_pos=set()
for i in P_and_T_dict.keys():
    if P_and_T_dict[i][13]=="FP":
        C_and_F_pos.add(i)


# F: (FP by both illumina and pacbio), intersection between FP by illumina and pacbio
F_pos=B_and_F_pos.intersection(C_and_F_pos)
F_pos=list(F_pos)
print('F size',len(F_pos))

F=[]
for i in range(len(F_pos)):
    F_element= I_and_T_dict[F_pos[i]] + P_and_T_dict[F_pos[i]]
    F.append(F_element)

#B
B_pos= B_and_F_pos.difference(F_pos)
B_pos=list(B_pos)
print('B size',len(B_pos))

B=[]
for i in range(len(B_pos)):
    B_element= I_and_T_dict[B_pos[i]]
    B.append(B_element)

#C
C_pos= C_and_F_pos.difference(F_pos)
C_pos=list(C_pos)
print('C size',len(C_pos))

C=[]
for i in range(len(C_pos)):
    C_element= P_and_T_dict[C_pos[i]]
    C.append(C_element)

# D and G: (TP by illumina), TP by truth or query decision
D_and_G_pos=set()
for i in I_and_T_dict.keys():
    if I_and_T_dict[i][13]=="TP":
        D_and_G_pos.add(i)
    # if I_and_T_dict[i][6]=="TP":
    #     D_and_G_pos.add(i)

# E and G: (TP by pacbio), TP by truth or query decision
E_and_G_pos=set()
for i in P_and_T_dict.keys():
    if P_and_T_dict[i][13]=="TP":
        E_and_G_pos.add(i)
    # if P_and_T_dict[i][6]=="TP":
    #     E_and_G_pos.add(i)


#G : (TP by both illumina and pacbio), intersection
G_pos=D_and_G_pos.intersection(E_and_G_pos)
G_pos=list(G_pos)
print('G size',len(G_pos))

G=[]
for i in range(len(G_pos)):
    G_element= I_and_T_dict[G_pos[i]] + P_and_T_dict[G_pos[i]]
    G.append(G_element)

#D
D_pos=D_and_G_pos.difference(G_pos)
D_pos=list(D_pos)
print('D size',len(D_pos))
D=[]
for i in range(len(D_pos)):
    D_element= I_and_T_dict[D_pos[i]]
    D.append(D_element)

#E
E_pos=E_and_G_pos.difference(G_pos)
E_pos=list(E_pos)
print('E size',len(E_pos))
E=[]
for i in range(len(E_pos)):
    E_element= P_and_T_dict[E_pos[i]]
    E.append(E_element)

# This function checks ratio of SNP and INDEL
def percentage(x):
    insertion_count=0
    deletion_count=0
    snp_count=0
    snp=[]
    for i in x:
        if i[17]=='INDEL':
            if len(i[2])>len(i[3]):
                deletion_count+=1
            if len(i[2])<len(i[3]):
                insertion_count+=1
        if i[17]=='SNP':
            snp_count+=1
            snp.append(i[1])
    insertion_percent=round(insertion_count/len(x),4)
    deletion_percent=round(deletion_count/len(x),4)
    snp_percent=round(snp_count/len(x),4)
    return snp,insertion_percent,deletion_percent,snp_percent,insertion_count,deletion_count,snp_count

# check percentage of SNP and INDEL in these regions
# aa,ba,ca,da,ea,f=percentage(A)
# ab,bb,cb,db,eb,fb=percentage(B)
snp_c,ac,bc,cc,dc,ec,fc=percentage(C)
# ad,bd,cd,dd,ed,fd=percentage(D)
snp_e,ae,be,ce,de,ee,fe=percentage(E)
# af,bf,cf,df,ef,ff=percentage(F)
# ag,bg,cg,dg,eg,fg=percentage(G)
print("percentage of INDEL in region C",ac+bc)
print("number of SNP in region C",fc)
print("percentage of SNP in region E",ce)
print("number of SNP in region E: ",fe)

# extact Region C and E SNPs
snp_c=list(map(int,snp_c))
snp_e=list(map(int,snp_e))
snp_c.sort()
snp_e.sort()

#find genotype quality value
gq_c=[]
snp_c_dict={}
for j in snp_c:
    for i in C:
        if j==int(i[1]):
            snp_c_dict[j]=[float(i[15])]
            gq_c.append(float(i[15]))

gq_e=[]
snp_e_dict={}
for j in snp_e:
    for i in E:
        if j==int(i[1]):
            snp_e_dict[j]=[float(i[15])]
            gq_e.append(float(i[15]))


# draw distribution for Genotype Quality
import matplotlib.pyplot as plt

plt.hist(gq_c,bins='auto',alpha=0.6,label='Region C SNP')
plt.hist(gq_e,bins='auto',alpha=0.6,label='Region E SNP')
plt.title('Genotype Quality distribution')
plt.xlabel('Genotype Quality')
plt.ylabel('Count')
plt.legend()

plt.show()


# identify whether in repeats or not
file=open("../../../raw_data/annotation/hg19RepeatMasker.txt", 'r')

RM=[]
for line in file:
    tmp=[]
    for element in line.split():
        tmp.append(element)
    RM.append(tmp)

chr20_RM=[]

for i in range(len(RM)):
    if RM[i][0]=='chr20':
        chr20_RM.append(RM[i][1:3])
for i in range(len(chr20_RM)):
    chr20_RM[i]=list(map(int,chr20_RM[i]))

for i in snp_c:
    for j in chr20_RM:
        if i>=j[0] and i<=j[1]:
            snp_c_dict[i].append(1)
            break
    if len(snp_c_dict[i])==1:
        snp_c_dict[i].append(0)

for i in snp_e:
    for j in chr20_RM:
        if i>=j[0] and i<=j[1]:
            snp_e_dict[i].append(1)
            break
    if len(snp_e_dict[i])==1:
        snp_e_dict[i].append(0)

#identify whether in SD or not
file=open("../../../raw_data/annotation/MosaicSDs_Human_hg19.txt", 'r')

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

for i in snp_c:
    for j in chr20_SD:
        if i>=j[0] and i<=j[1]:
            snp_c_dict[i].append(1)
            break
    if len(snp_c_dict[i])==2:
        snp_c_dict[i].append(0)

for i in snp_e:
    for j in chr20_SD:
        if i>=j[0] and i<=j[1]:
            snp_e_dict[i].append(1)
            break
    if len(snp_e_dict[i])==2:
        snp_e_dict[i].append(0)

# add DP value (Read depth) and variant allele fraction
file2=open("pacbio.simulate.output.csv",'r')
reader2=csv.reader(file2)
P=list(reader2)
P.pop(0)
rd_c=[]
vaf_c=[]
for i in snp_c:
    for j in P:
        if i==int(j[1]):
            snp_c_dict[i].append(int(j[10]))
            rd_c.append(int(j[10]))
            snp_c_dict[i].append(float(j[12]))
            vaf_c.append(float(j[12]))
rd_e=[]
vaf_e=[]
for i in snp_e:
    for j in P:
        if i==int(j[1]):
            snp_e_dict[i].append(int(j[10]))
            rd_e.append(int(j[10]))
            snp_e_dict[i].append(float(j[12]))
            vaf_e.append(float(j[12]))

# draw distribution for Read Depth
import matplotlib.pyplot as plt

plt.hist(rd_c,bins='auto',alpha=0.6,label='Region C SNP')
plt.hist(rd_e,bins='auto',alpha=0.6,label='Region E SNP')
plt.title('Read depth distribution')
plt.xlabel('Read depth')
plt.ylabel('Count')
plt.legend()

plt.show()

# draw distribution for allele fraction
import matplotlib.pyplot as plt

plt.hist(vaf_c,bins='auto',alpha=0.6,label='Region C SNP')
plt.hist(vaf_e,bins='auto',alpha=0.6,label='Region E SNP')
plt.title('Variant allele fraction distribution')
plt.xlabel('Variant allele fraction')
plt.ylabel('Count')
plt.legend()

plt.show()

# add label for prediction
for i in snp_c:
    snp_c_dict[i].append(0)
for i in snp_e:
    snp_e_dict[i].append(1)


#combine and add position
snp_dict=dict(list(snp_c_dict.items())+list(snp_e_dict.items()))
snp=[]
for i in snp_dict.keys():
    snp_dict[i].insert(0, i)
    snp.append(snp_dict[i])

#write to csv
with open("C_E_dataset.csv","w") as file:
    writer=csv.writer(file)
    writer.writerow(['pos','qual','repeats','sd','read depth','allele fraction','truth'])
    writer.writerows(snp)
