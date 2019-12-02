"""
Read CSV variant calls from Illumina, PacBio, and Truth
Divide all three files by intersections into 7 regions
_experiment_= "1"
__author__ = “Jiajia Xu”
__copyright__ = “Copyright 2019, Individual Computing Project”
__email__ = “jiajia.xu@anu.edu.au”
__status__ = “Production”
"""


# read CSV variant calls files
import csv
import numpy as np

with open("wgs.chr20.csv","r") as file:
    reader=csv.reader(file)
    I_and_T=list(reader)
    I_and_T.pop(0)
I_and_T_array=np.array(I_and_T)
I_and_T_pos= I_and_T_array[:, 1]
I_and_T_pos=list(map(int, I_and_T_pos.tolist()))

with open("pacbio.chr20.csv","r") as file:
    reader=csv.reader(file)
    P_and_T=list(reader)
    P_and_T.pop(0)
P_and_T_array=np.array(P_and_T)
P_and_T_pos= P_and_T_array[:, 1]
P_and_T_pos=list(map(int, P_and_T_pos.tolist()))

with open("cleaned_truth.chr20.csv","r") as file:
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
    if I_and_T_dict[i][13]=="TP" or I_and_T_dict[i][6]=="TP":
        D_and_G_pos.add(i)

# E and G: (TP by pacbio), TP by truth or query decision
E_and_G_pos=set()
for i in P_and_T_dict.keys():
    if P_and_T_dict[i][13]=="TP" or P_and_T_dict[i][6]=="TP":
        E_and_G_pos.add(i)

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
    for i in x:
        if i[17]=='INDEL':
            if len(i[2])>len(i[3]):
                deletion_count+=1
            if len(i[2])<len(i[3]):
                insertion_count+=1
        if i[17]=='SNP':
            snp_count+=1
    insertion_percent=round(insertion_count/len(x),4)
    deletion_percent=round(deletion_count/len(x),4)
    snp_percent=round(snp_count/len(x),4)
    return insertion_percent,deletion_percent,snp_percent


# check percentage of SNP and INDEL in these regions
a,b,c=percentage(C)
print('INSERTION percentage in region C:', a)
print('DELETION percentage in region C:', b)
print('SNP percentage in region C:', c)
a,b,c=percentage(E)
print('INSERTION percentage in region E:', a)
print('DELETION percentage in region E:', b)
print('SNP percentage in region E:', c)
a,b,c=percentage(F)
print('INSERTION percentage in region F:', a)
print('DELETION percentage in region F:', b)
print('SNP percentage in region F:', c)


# check percentage of confident region in these region
count=0
for i in range(len(T_pos)):
    if 'CONF' in I_and_T_dict[T_pos[i]][4]:
        count+=1
print('confident region in truth vcf: ',count/len(T_pos))

count1=0
for i in range(len(T_pos)):
    if 'CONF' in P_and_T_dict[T_pos[i]][4]:
        count1+=1

print('confident region in truth vcf: ',count1/len(T_pos))


count_A=0
for i in range(len(A_pos)):
    if 'CONF' in I_and_T_dict[A_pos[i]][4]:
        count_A+=1
print('confident region in A: ',count_A/len(A_pos))


count_B=0
for i in range(len(B_pos)):
    if 'CONF' in I_and_T_dict[B_pos[i]][4]:
        count_B+=1
print('confident region in B: ',count_B/len(B_pos))

count_C=0
for i in range(len(C_pos)):
        if 'CONF' in P_and_T_dict[C_pos[i]][4]:
            count_C+=1
print('confident region in C: ',count_C/len(C_pos))


count_D=0
for i in range(len(D_pos)):
        if 'CONF' in I_and_T_dict[D_pos[i]][4]:
            count_D+=1
print('confident region in D: ',count_D/len(D_pos))

count_E=0
for i in range(len(E_pos)):
        if 'CONF' in P_and_T_dict[E_pos[i]][4]:
            count_E+=1
print('confident region in E: ',count_E/len(E_pos))


count_F=0
for i in range(len(F_pos)):
        if 'CONF' in I_and_T_dict[F_pos[i]][4]:
            count_F+=1
print('confident region in F: ',count_F/len(F_pos))

count_G=0
for i in range(len(G_pos)):
        if 'CONF' in I_and_T_dict[G_pos[i]][4]:
            count_G+=1
print('confident region in G: ',count_G/len(G_pos))

