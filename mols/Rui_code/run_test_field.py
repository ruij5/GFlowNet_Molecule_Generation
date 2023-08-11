import pandas as pd
import pickle
import os
import json
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import BRICS
from rdkit.Chem import BRICS
from rdkit.Chem.Scaffolds import MurckoScaffold

# m = Chem.MolFromSmiles('CSCCC(NCc1ccc2c(c1)CN(C)C2)c1nnc2ccccn12')
# smi = Chem.MolToSmiles(m, 1)
# res = list(BRICS.FindBRICSBonds(m))
# for tuple in res:
#     print(type(tuple[0][0]))
#     for i in tuple:
#         print(i[0])

# a = "*c1cn[nH]c1*"
# print(a.replace("*", ""))

# mol = Chem.MolFromSmiles('CC(=CC=N)C1CC=C(CCc2cnn(-c3ccccn3)c2C(=N)N)CC1')
# fragment = '[14*]c1ccccn1'
# mol.GetSubstructMatch(fragment)






smi="c1c[nH]c2ccccc12"
# smi_2= MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles("n1cncn1"))
# smi_2 = "*C1CCN(*)CC1"
# m = Chem.MolFromSmiles(smi)
mol = Chem.MolFromSmiles(smi)
n=0
for atom in mol.GetAtoms():
    n = n + 1
print(n)
smi1 = Chem.MolToSmiles(mol,isomericSmiles=False)
# m_2 = Chem.MolFromSmiles(smi_2)
# smi=Chem.MolToSmiles(m, kekuleSmiles = True)
# smi_2=Chem.MolToSmiles(m_2, rootedAtAtom = 0)
# for i, atom in enumerate(m.GetAtoms()):
#     if (atom.GetIdx() == 1):
#         em1 = Chem.EditableMol(m)
#         print(atom.GetIdx())
#         em1 = em1.RemoveAtom(atom.GetIdx())
#         mol = em1.GetMol()
#         print(Chem.MolToSmiles(mol))

print(smi1)


smi2="c1c[nH]c2ccccc12"
mol2 = Chem.MolFromSmiles(smi2)
smi2 = Chem.MolToSmiles(mol2,isomericSmiles=False)
print(smi2)