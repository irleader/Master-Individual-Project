"""
Check proportion of segmental duplications in each region
_experiment_= "2"
__author__ = “Jiajia Xu”
__copyright__ = “Copyright 2019, Individual Computing Project”
__email__ = “jiajia.xu@anu.edu.au”
__status__ = “Production”
"""

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
    if P_and_T[i][18]!='nocall':
        P_pos.append(int(P_and_T[i][1]))
print('pacbio set size: ',len(P_pos))
I_pos=[]
for i in range(len(I_and_T_pos)):
    if I_and_T[i][18]!='nocall':
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
    if I_and_T_dict[i][6]=="TP":
        D_and_G_pos.add(i)

# E and G: (TP by pacbio), TP by truth or query decision
E_and_G_pos=set()
for i in P_and_T_dict.keys():
    if P_and_T_dict[i][13]=="TP":
        E_and_G_pos.add(i)
    if P_and_T_dict[i][6]=="TP":
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


# find all segmental duplications in each region
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

# c_A=0
# sd_A_pos=[]
# for i in range(len(A_pos)):
#     for j in range(len(chr20_SD)):
#         if A_pos[i]>= chr20_SD[j][0] and A_pos[i]<= chr20_SD[j][1]:
#             sd_A_pos.append(A_pos[i])
#             c_A+=1


def find_sd(x):
    y=[]
    c=0
    for i in range(len(x)):
        for j in range(len(chr20_SD)):
            if x[i]>= chr20_SD[j][0] and x[i]<= chr20_SD[j][1]:
                y.append(x[i])
                c=c+1
                break
    return y,c

sd_A_pos,c_A=find_sd(A_pos)
sd_B_pos,c_B=find_sd(B_pos)
sd_C_pos,c_C=find_sd(C_pos)
sd_D_pos,c_D=find_sd(D_pos)
sd_E_pos,c_E=find_sd(E_pos)
sd_F_pos,c_F=find_sd(F_pos)
sd_G_pos,c_G=find_sd(G_pos)
sd_total_pos,c_total=find_sd(total_pos)
sd_T_pos,c_T=find_sd(T_pos)
sd_I_pos,c_I=find_sd(I_pos)
sd_P_pos,c_P=find_sd(P_pos)

print(c_A)
print(c_B)
print(c_C)
print(c_D)
print(c_E)
print(c_F)
print(c_G)
# print(c_total)
# print(c_T)
# print(c_I)
# print(c_P)