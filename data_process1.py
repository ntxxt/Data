from dataclasses import field
from genericpath import isdir
import os
from os import walk

"""
#summary of DUD-E
folder = './all'
sub_folders = [name for name in os.listdir(folder) if os.path.isdir(os.path.join(folder, name))]
with open('DUE-E.txt', 'w') as f:
    for i in sub_folders:
        f.write("%s\n" % i)
    print('done')

#summary of LIT-PCAB
dir = './'
def listfiles(dir):
    files = os.listdir(dir)
    allfile = list()
    for entry in files:
        fullpath = os.path.join(dir, entry)
        if os.path.isdir(fullpath):
            allfile = allfile + listfiles(fullpath)
        else:
            allfile.append(fullpath)
    return allfile

allfile = listfiles(dir)
proteins = []
for i in allfile:
    try:
        protein = i.split('/')[2].split('_')[0]
        if len(protein) == 4 :
            proteins.append(protein)
    except: pass
#print(proteins)

def write_list(list, filename):
    with open(filename, 'w') as f:
        for i in list:
            f.write("%s\n" % i)
    print('done')

write_list(proteins, 'LIT-PCAB.txt')


#summary of CASF-2013
folder = '/Users/xiaotongxu/structural_datasets/CASF/CASF-2013/coreset'
sub_folders = [name for name in os.listdir(folder) if os.path.isdir(os.path.join(folder, name))]
with open('CASF-2013.txt', 'w') as f:
    for i in sub_folders:
        f.write("%s\n" % i)
    print('done')
"""

#process Binding_MOAD data file
import pandas as pd
file_path = '/Users/xiaotongxu/structural_datasets/Binding_MOAD'
name = 'nr_bind.csv'
with open(os.path.join(file_path, name), 'r') as f:
    lines = f.readlines()
lines_ = [line for line in lines if "Representative" in line]
PDB_ID = [lines_[i].split()[-1].replace(",","") for i in range(len(lines_))]
id_line = []
for i in range(len(lines)):
    if "Representative" in lines[i-1]:
      #if lines[i].split(",")[-7] == "valid":
      id_line.append(lines[i])
ligand_smile = [id_line[i].split(",")[-2] for i in range(len(id_line))]
unit = [id_line[i].split(",")[-3] for i in range(len(id_line))]
activity = [id_line[i].split(",")[-4] for i in range(len(id_line))]
valid = [id_line[i].split(",")[-7] for i in range(len(id_line))]
#print(len(PDB_ID), len(valid), len(ligand_smile), len(unit), len(activity))

combine = pd.DataFrame(
    {'PDB_ID': PDB_ID, 'ligand_smile': ligand_smile, 'unit': unit, 'Activity': activity, "Valid": valid}
)
combine.to_csv('Binding_MOAD.csv')


#not working dur to strange csv file formatting 
"""
index = []
for i in range(len(first)):
    try:
        if first[i-1] != 0 and first[i] == 0:
            index.append(i)
    except: pass
PDB_ID = file.iloc[:,2][index]

id = index[0]
id2 = index[1]
valid = file.iloc[:,3]
#print(valid[:20])
#print(id,id2)
#for i in range(id,id2):
    #print(i)
    
    #if file.iloc[:3][i] != 0 and file.iloc[:4][i] == 'valid':
        #print(i)
pdb = file.iloc[:,2][id]
#print(pdb)
"""
