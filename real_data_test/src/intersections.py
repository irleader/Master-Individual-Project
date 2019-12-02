# 7 sections ABCDEFG
import csv
import numpy as np

with open("wgs.chr20.csv","r") as file:
    reader=csv.reader(file)
    I_and_T=list(reader)
    I_and_T.pop(0)

with open("pacbio.chr20.csv","r") as file:
    reader=csv.reader(file)
    P_and_T=list(reader)
    P_and_T.pop(0)

with open("cleaned_truth.chr20.csv","r") as file:
    reader=csv.reader(file)
    T=list(reader)
    T.pop(0)
T_array=np.array(T)


I_and_T_array=np.array(I_and_T)
I_and_T_pos= I_and_T_array[:, 1]
# I_info= I_array[:, 4]
# I_t_location= I_array[:, 11]
# #I_t_genotype=I_array[:,5]
# I_t_deicision= I_array[:, 6]
# I_q_location= I_array[:, 18]
# #I_q_genotype=I_array[:,12]
# I_q_deicision= I_array[:, 13]


P_and_T_array=np.array(P_and_T)
P_and_T_pos= P_and_T_array[:, 1]
# P_info= P_array[:, 4]
# P_t_location= I_array[:, 11]
# #I_t_genotype=I_array[:,5]
# P_t_deicision= I_array[:, 6]
# P_q_location= I_array[:, 18]
# #I_q_genotype=I_array[:,12]
# P_q_deicision= I_array[:, 13]

T_pos=T_array[:,1]
T_pos=list(map(int,T_pos.tolist()))
I_and_T_pos=list(map(int, I_and_T_pos.tolist()))
P_and_T_pos=list(map(int, P_and_T_pos.tolist()))
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


# F: called by both illumina and pacbio but not in truth
F_pos=list(set(P_and_T_pos).intersection(I_and_T_pos).difference(T_pos))

print('F size',len(F_pos))
F=[]
for i in range(len(F_pos)):
    F_element= I_and_T_dict[F_pos[i]] + P_and_T_dict[F_pos[i]]
    F.append(F_element)

with open("F.csv","w") as file:
    writer=csv.writer(file)
    writer.writerow(['chrom','pos','ref','alt','info','t_genotype','t_decision','t_decision_sub','t_quality','t_additional_info','t_variant_type','t_location_type',
                  'q_genotype','q_decision','q_decision_sub','q_quality','q_additional_info','q_variant_type','q_location_type',
                     'chrom','pos','ref','alt','info','t_genotype','t_decision','t_decision_sub','t_quality','t_additional_info','t_variant_type','t_location_type',
                  'q_genotype','q_decision','q_decision_sub','q_quality','q_additional_info','q_variant_type','q_location_type'])
    writer.writerows(F)

# B: called only by illumina
B_pos= list(set(I_and_T_pos).difference(T_pos).difference(P_and_T_pos))
print('B size',len(B_pos))
B=[]
for i in range(len(B_pos)):
    B_element= I_and_T_dict[B_pos[i]]
    B.append(B_element)

with open("B.csv","w") as file:
    writer=csv.writer(file)
    writer.writerow(['chrom','pos','ref','alt','info','t_genotype','t_decision','t_decision_sub','t_quality','t_additional_info','t_variant_type','t_location_type',
                  'q_genotype','q_decision','q_decision_sub','q_quality','q_additional_info','q_variant_type','q_location_type'])
    writer.writerows(B)


# C: called only by pacbio
C_pos=list(set(P_and_T_pos).difference(T_pos).difference(I_and_T_pos))
print('C size',len(C_pos))
C=[]
for i in range(len(C_pos)):
    C_element= P_and_T_dict[C_pos[i]]
    C.append(C_element)

with open("C.csv","w") as file:
    writer=csv.writer(file)
    writer.writerow(['chrom','pos','ref','alt','info','t_genotype','t_decision','t_decision_sub','t_quality','t_additional_info','t_variant_type','t_location_type',
                  'q_genotype','q_decision','q_decision_sub','q_quality','q_additional_info','q_variant_type','q_location_type'])
    writer.writerows(C)





# A: called only by truth
A_pos=list(set(T_pos).difference(I_pos).difference(P_pos))
print('A size',len(A_pos))
A=[]
for i in range(len(A_pos)):
    A_element= I_and_T_dict[A_pos[i]] + P_and_T_dict[A_pos[i]]
    A.append(A_element)

# D: called by illumina and truth but not pacbio
D_pos=list(set(I_pos).intersection(T_pos).difference(P_pos))
print('D size',len(D_pos))
D=[]
for i in range(len(D_pos)):
    D_element= I_and_T_dict[D_pos[i]]
    D.append(D_element)

# E: called by truth and pacbio but not illumina
E_pos=list(set(T_pos).intersection(P_pos).difference(I_pos))
print('E size',len(E_pos))
E=[]
for i in range(len(E_pos)):
    E_element= P_and_T_dict[E_pos[i]]
    E.append(E_element)

# G: called by all three
G_pos=list(set(T_pos).intersection(I_pos).intersection(P_pos))
print('G size',len(G_pos))
G=[]
for i in range(len(G_pos)):
    G_element= I_and_T_dict[G_pos[i]] + P_and_T_dict[G_pos[i]]
    G.append(G_element)



count=0
for i in range(len(T_pos)):
    if 'CONF' in I_and_T_dict[T_pos[i]][4]:
        count+=1

count1=0
for i in range(len(I_and_T)):
    if I_and_T[i][11]!='nocall':
        count1+=1

print('confident region in truth vcf: ',count/count1)

count_A=0
for i in range(len(A_pos)):
    if I_and_T_dict[A_pos[i]][11]!='nocall' and I_and_T_dict[A_pos[i]][18]=='nocall' and P_and_T_dict[A_pos[i]][18]=='nocall':
        if 'CONF' in I_and_T_dict[A_pos[i]][4]:
            count_A+=1
print('confident region in A: ',count_A/len(A_pos))


count_B=0
for i in range(len(B_pos)):
    if I_and_T_dict[B_pos[i]][11]=='nocall' and I_and_T_dict[B_pos[i]][18]!='nocall':
        if 'CONF' in I_and_T_dict[B_pos[i]][4]:
            count_B+=1
print('confident region in B: ',count_B/len(B_pos))

count_C=0
for i in range(len(C_pos)):
    if P_and_T_dict[C_pos[i]][11]=='nocall' and P_and_T_dict[C_pos[i]][18]!='nocall':
        if 'CONF' in P_and_T_dict[C_pos[i]][4]:
            count_C+=1
print('confident region in C: ',count_C/len(C_pos))


count_D=0
for i in range(len(D_pos)):
    if I_and_T_dict[D_pos[i]][11] != 'nocall' and I_and_T_dict[D_pos[i]][18] != 'nocall':
        if 'CONF' in I_and_T_dict[D_pos[i]][4]:
            count_D+=1
print('confident region in D: ',count_D/len(D_pos))

count_E=0
for i in range(len(E_pos)):
    if P_and_T_dict[E_pos[i]][11] != 'nocall' and P_and_T_dict[E_pos[i]][18] != 'nocall':
        if 'CONF' in P_and_T_dict[E_pos[i]][4]:
            count_E+=1
print('confident region in E: ',count_E/len(E_pos))


count_F=0
for i in range(len(F_pos)):
    if I_and_T_dict[F_pos[i]][11]=='nocall' and I_and_T_dict[F_pos[i]][18]!='nocall' and P_and_T_dict[F_pos[i]][18]!='nocall':
        if 'CONF' in I_and_T_dict[F_pos[i]][4]:
            count_F+=1
print('confident region in F: ',count_F/len(F_pos))

count_G=0
for i in range(len(G_pos)):
    if I_and_T_dict[G_pos[i]][11]!='nocall' and I_and_T_dict[G_pos[i]][18]!='nocall' and P_and_T_dict[G_pos[i]][18]!='nocall':
        if 'CONF' in I_and_T_dict[G_pos[i]][4]:
            count_G+=1
print('confident region in G: ',count_G/len(G_pos))