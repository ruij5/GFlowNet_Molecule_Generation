import rdkit.Chem as Chem
# for bond in mol.GetBonds():
#     a1 = bond.GetBeginAtom().GetIdx()
#     a2 = bond.GetEndAtom().GetIdx()
#     t = bond.GetBondType()
#     print(t)



import pandas as pd
mol = Chem.MolFromSmiles('c1ccccc1')


# m1 = Chem.MolToSmiles(mol, kekuleSmiles=True)
# mol2 = Chem.MolFromSmiles(m1)
# m3 = Chem.MolToSmiles(mol2, kekuleSmiles=True)
# m2 = Chem.MolToSmiles(mol, kekuleSmiles=False)
# print(m1)
# print(m2)
# print(mol2)
# print(m3)

smile = Chem.MolFragmentToSmiles(mol, [0,1,2,3,4,5], kekuleSmiles=False)
print(smile)
smile = Chem.MolFragmentToSmiles(mol, [0,1,2,3,4,5], kekuleSmiles=True)
print(smile)