import pandas as pd
import pickle
import os
import json
from rdkit import Chem
from rdkit.Chem import Draw

my_file = open("./mols/data/vocabulary.txt", "r")
data = my_file.read()
data_into_list = data.replace('\n', ' ').split(" ")

block_names = {"block_name":{}}
n = 0
for i in data_into_list:
    block_names["block_name"][str(n)] = i
    print(n)
    n += 1


my_json_str = json.dumps(block_names)
with open('./mols/data/new_blocks_PDB_105.json', 'w') as file:
   file.write(my_json_str)

my_file = open("./mols/data/vocabulary.txt", "r")
data = my_file.read()
data_into_list = data.replace('\n', ' ').split(" ")
block_rs = {"block_r":{}}
n = 0
for i in data_into_list:
    smi = i
    mol = Chem.MolFromSmiles(smi)
    r_block = []
    for i, atom in enumerate(mol.GetAtoms()):
        idx = atom.GetIdx()
        r_block.append(idx)
    block_rs["block_r"][str(n)] = r_block
    n += 1


my_json_str = json.dumps(block_rs)

with open('./mols/data/new_new_blocks_PDB_105.json', 'w') as file:
   file.write(my_json_str)