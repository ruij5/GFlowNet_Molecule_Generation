import pandas as pd
import mol_mdp_ext
import pickle
import os
import json
from rdkit import Chem
from rdkit.Chem import Draw
import numpy as np
from gflownet import Dataset as dataset
from mol_mdp_ext import MolMDPExtended, BlockMoleculeDataExtended
from utils.molMDP import BlockMoleculeData, MolMDP
import utils.chem as chem
from rdkit import Chem

my_file = open("./mols/data/vocabulary.txt", "r")
data = my_file.read()
data_into_list = data.replace('\n', ' ').split(" ")

# block_names = {"block_name":{}}
# n = 0
# for i in data_into_list:
#     block_names["block_name"][str(n)] = i
#     print(n)
#     n += 1

gold = Chem.MolFromSmiles('[Au]')
block_json = {"block_name":{},"block_smi":{},"block_r":{}}
n = 0
for smi in data_into_list:

    mol = Chem.MolFromSmiles(smi)
    r_block = []
    for i, atom in enumerate(mol.GetAtoms()):
        try:
            idx = atom.GetIdx()
            molA, _ = chem.mol_from_frag(jun_bonds=[[0,1,0,idx]],frags=[gold, mol])
            r_block.append(idx)
        except:
            print(idx)
    for index in range(len(r_block)):
        block_json["block_name"][str(n + index)] = smi
        block_json["block_smi"][str(n + index)] = smi
        new_r_block = []
        new_r_block.append(r_block[index])
        for idx in range(len(r_block)):
            if (idx != index):
                new_r_block.append(r_block[idx])
        block_json["block_r"][str(n + index)] = new_r_block

    n += len(r_block)










# blocks = pd.read_json("mols/data/new_blocks_PDB_105.json")
# # blocks = pd.read_json("mols/data/new_blocks_PDB_105.json")
# block_smi = blocks["block_smi"].to_list()
# block_rs = blocks["block_r"].to_list()
# block_nrs = np.asarray([len(r) for r in block_rs])
# block_mols = [Chem.MolFromSmiles(smi) for smi in blocks["block_smi"]]
# block_natm = np.asarray([b.GetNumAtoms() for b in block_mols])


# translation_table = {}
# for blockidx in range(len(block_mols)):
#     # print(blockidx)
#     atom_map = {}
#     for j in range(len(block_mols)):
#         if block_smi[blockidx] == block_smi[j]:
#             atom_map[block_rs[j][0]] = j
#     # print(atom_map)
#     translation_table[blockidx] = atom_map
# gold = Chem.MolFromSmiles('[Au]')

# for blockidx in range(len(block_mols)):
#     for j in block_rs[blockidx]:
#         # print("start")
#         # print(blockidx)
#         # print(j)
#         # print(translation_table[blockidx])
#             try:
#                 molA, _ = chem.mol_from_frag(jun_bonds=[[0,1,0,j]],frags=[gold, block_mols[blockidx]])
#             except:
#                 print(j)
#             # print(block_duplicate)
#             # print(atom)



my_json_str = json.dumps(block_json)

with open('./mols/data/new_new_blocks_PDB_105.json', 'w') as file:
   file.write(my_json_str)