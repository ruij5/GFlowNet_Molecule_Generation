import pandas as pd
import pickle
import os
import json
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import BRICS
import csv
df = pd.read_json('./mols/data/new_new_new.json')

block = df["block_name"].to_list()
# print(block)
# print(block)
# 
# 
# smi_block_r_dict = {}
# block_name = []
# block_smi = []
# block_r = []
# for temp in block:
#     block_r_temp = []
#     # print(temp)
#     for smi in temp:
#         print(smi)
#         mol = Chem.MolFromSmiles(smi)
#         smi = Chem.MolToSmiles(mol)
#         mol_temp=mol
#         bool_n = False
#         for i, atom in enumerate(mol.GetAtoms()):
#             smi_temp=Chem.MolToSmiles(mol, rootedAtAtom = atom.GetIdx())
#             smi_temp=smi_temp.replace("(*)", "")
#             smi_temp=smi_temp.replace("*", "")
#             # print(smi_temp)
#             if (smi_temp in block_smi):
#                 new_smi_temp=Chem.MolToSmiles(mol, rootedAtAtom = atom.GetIdx())
#                 mol_temp=Chem.MolFromSmiles(new_smi_temp)
#                 bool_n = True
#                 break
# 
#         if (bool_n == False):
#             smi_temp=smi.replace("(*)","")
#             smi_temp=smi_temp.replace("*", "")
#             block_smi.append(smi_temp)
# 
#         n = 0 # count how many free bonds
#         for i, atom in enumerate(mol_temp.GetAtoms()):
#             if (atom.GetSymbol() == "*"):
#                 if (atom.GetIdx() == 0):
#                     block_r_temp.append(0)
#                     n = n + 1
#                 else: 
#                     block_r_temp.append(atom.GetIdx() - n - 1)
#                     n = n + 1
#         block_r.append(block_r_temp)
# print(block_r)



smi_block_r_dict = {}
block_name = []
block_smi = []
block_r = []
for block_list in block:
    block_r_temp = []
    # print(temp)
    for smi in block_list:
        mol = Chem.MolFromSmiles(smi)
        smi = Chem.MolToSmiles(mol,isomericSmiles=False)
        mol_temp=mol
        bool_n = False
        for i, atom in enumerate(mol.GetAtoms()):
            smi_temp=Chem.MolToSmiles(mol, rootedAtAtom = atom.GetIdx())
            smi_temp=smi_temp.replace("(*)", "")
            smi_temp=smi_temp.replace("*", "")
            # print(smi_temp)
            if (smi_temp in block_smi):
                new_smi_temp=Chem.MolToSmiles(mol, rootedAtAtom = atom.GetIdx())
                mol_temp=Chem.MolFromSmiles(new_smi_temp)
                bool_n = True
                break

        if (bool_n == False):
            smi_temp=smi.replace("(*)","")
            smi_temp=smi_temp.replace("*", "")
            block_smi.append(smi_temp)

        n = 0 # count how many free bonds
        for i, atom in enumerate(mol_temp.GetAtoms()):
            if (atom.GetSymbol() == "*"):
                if (atom.GetIdx() == 0):
                    block_r_temp.append(0)
                    n = n + 1
                else: 
                    block_r_temp.append(atom.GetIdx() - n - 1)
                    n = n + 1
        block_r.append(block_r_temp)
print(block_r)







# create the json file


# block_json = {"block_name":{},"block_smi":{},"block_r":{}}

# for index in range(len(r_block)):
#         block_json["block_name"][str(n + index)] = smi
#         block_json["block_smi"][str(n + index)] = smi
#         new_r_block = []
#         new_r_block.append(r_block[index])
#         for idx in range(len(r_block)):
#             if (idx != index):
#                 new_r_block.append(r_block[idx])
#         block_json["block_r"][str(n + index)] = new_r_block



# my_json_str = json.dumps(block_json)

# with open('./mols/data/final_blocks_PDB_105.json', 'w') as file:
#    file.write(my_json_str)