import pandas as pd
import pickle
import os
import json
from rdkit import Chem
from rdkit.Chem import Draw
import numpy as np
import utils.chem as chem
from rdkit import Chem
import csv
from rdkit.Chem import BRICS

 


smi_list = []
with open("./mols/data/training_set.csv", 'r') as file:
  csvreader = csv.reader(file)
  for row in csvreader:
    smi_list.append(row)
# print(smi_list[0:10])
blocks = [] # molecule fragements after decomposition
# bonds = [] # broken bonds on the molecule
# block_r = [] # list: atom indicies
time = 0 # limit the inputs for test purpose


block_names = {"block_name":{}} # create a json file
num = 0

for smi in smi_list:
    mol = Chem.MolFromSmiles(smi[0])
    res_block = list(BRICS.BRICSDecompose(mol))
    # res_bond = list(BRICS.FindBRICSBonds(mol))
    res = [] # blocks
    # res_block_r = []
    
    for temp in sorted(res_block):
        new_mol=Chem.MolFromSmiles(temp)
        new_smi=Chem.MolToSmiles(new_mol,isomericSmiles=False)
        # bond_indicies = []
        # n = 0
        # get block_r
        # for i, atom in enumerate(new_mol.GetAtoms()):
        #     if (atom.GetSymbol() == "*"):
        #         if (atom.GetIdx() == 0):
        #             bond_indicies.append(0) 
        #             n = n + 1
        #         else: 
        #             bond_indicies.append(atom.GetIdx() - n - 1)
        #             n = n + 1
        # res_block_r.append(bond_indicies)
        # remove the free bonds
        # new_smi=new_smi.replace("(*)", "")
        # new_smi=new_smi.replace("*", "")
        res.append(new_smi)

    # block_r.append(res_block_r)
    blocks.append(res)
    block_names["block_name"][str(num)] = res
    # blocks.append(sorted(res_block))
    # bonds.append(sorted(res_bond))
    num = num + 1
    # time = time + 1
    # if time == 10:
    #     break
    
file.close()

# Create a dataframe for our bloc.csv

# dict = {'molecule': smi_list, 'block': blocks, 'bond': bonds, 'block_r': block_r}
# dict = {'molecule': smi_list, 'block': blocks}
# block_df = pd.DataFrame(dict)
# block_df.to_csv('./mols/data/Decomposition_blocks.csv', index=False)


my_json_str = json.dumps(block_names)

with open('./mols/data/new_new_new.json', 'w') as file:
   file.write(my_json_str)