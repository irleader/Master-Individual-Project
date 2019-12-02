import csv

with open("wgs.chr20.csv","r") as file:
    reader=csv.reader(file)
    variants=list(reader)
    variants.pop(0)

total_count=len(variants)

import numpy as np

#decision mismatch
selected_variants1=[]
variants_array=np.array(variants)
t_decision= variants_array[:, 6]
q_decision= variants_array[:, 13]

for i in range(len(t_decision)):
    if t_decision[i]!=q_decision[i]:
        selected_variants1.append(variants[i])

import variants_classification as vc
vc.write_csv('decision_mismatch.csv',selected_variants1)

#q_decision is FP
